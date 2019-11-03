#!/usr/bin/env python3
# Convert Scotch output to VCF format
# Called by doPredict.sh as
# 	python makeVCF.py [tsv results] [fasta ref] [vcf results stub]
# Writes the full results in VCF format to ${stub}.vcf,

import csv
import os
from pathlib import Path
import pysam
import textwrap
from typing import Any, Dict, List
import typing
import subprocess
import sys

# constants
CHROMS = list(str(c) for c in range(1, 23)) + ["X", "Y"]

# types of Scotch calls
PRED_TYPES = ["<DEL_L>", "<DEL_R>", "<INS>"]
# types of calls for which we subtract one from the position, to adjust indexing

OUTPUT_DELIMITER = "\t"
# vcf fields
ID = "."
QUAL = "100"
FILTER = "PASS"
INFO_TEMPLATE = "CALLED_BY={}" 
FORMAT = "GT"
GT = "./."
ENCODE_GT = "0/1"

# run a script
def run_script(script_name, *args):
	scotch_dir: Path = Path(__file__).absolute().parent
	script: Path = scotch_dir / script_name
	str_args = [str(a) for a in args]
	print(f"Running {script} with {str_args}")

	if script.suffix != ".py":
		subprocess.call([script] + str_args)
	else:
		subprocess.call(["python", script] + str_args)

# get chromosome lengths from fasta reference for ##contig headers
def get_chrom_lengths(fasta: Any) -> Dict[str, int]:
	return {chrom: fasta.get_reference_length(chrom) for chrom in CHROMS}

# write_header to provided csv writer
# returns number of header lines
def write_header(writer: Any, chrom_lengths: Dict[str, int]) -> None:
	
	# generic vcf headers
	headers: [str] = ["##fileformat=VCFv4.1"]
	headers.append("##phasing=none")
	headers.append("##ALT=<ID=DEL_L, Description=\"Deletion Start\">")
	headers.append("##ALT=<ID=DEL_R, Description=\"Deletion End\">")
	headers.append("##ALT=<ID=INS,Description=\"Insertion\">")
	headers.append("##INFO=<ID=CALLED_BY,Number=1,Type=String,Description=\"Indel callers that made a correlating call\">")
	headers.append("##FORMAT=<ID=GT,Number=1, Type=String,Description=\"Genotype\">")

	for header in headers:
		writer.writerow([header])

	# contig headers
	for chrom in CHROMS:
		chrom_length: int = chrom_lengths[chrom]
		writer.writerow([f"##contig=<ID={chrom},length={chrom_length}>"])

	# column headers
	writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])

# get nucleotides from FASTA at position
def get_nucs(ref: Any, chrom: str, start: int, end: int = None) -> str:
	if not end: 
		end = start + 1
	
	# pysam is 0-based (https://pysam.readthedocs.io/en/latest/glossary.html#term-region) so subtract 1
	zero_based_start: int = start - 1
	zero_based_end: int = end - 1
	return ref.fetch(chrom, zero_based_start, zero_based_end).upper()

# write a variant with given fields to a list of writers
def write_variant(writer: Any, chrom: str, pos: int, ref: str, alt: str, info: str, gt: str) -> None:
	variant_row = [chrom, pos, ID, ref, alt, QUAL, FILTER, info, FORMAT, gt]
	writer.writerow(variant_row)

# process variant, writing to VCFs
def process_variant(variant: List[str], writer: Any, fasta: Any) -> None:

	# unpack fields
	[chrom, raw_pos, pred_type, length, called_by] = variant
	pos: int = int(raw_pos)
	info: str = INFO_TEMPLATE.format(called_by)
	assert pred_type in PRED_TYPES, f"Variant at {chrom}:{pos} has unexpected type {pred_type}"

	# get ref
	ref: str = get_nucs(fasta, chrom, pos)

	# write results to standard vcf
	write_variant(writer, chrom, pos, ref, pred_type, info, GT)

# sort numerically by position
def sort_output(output_tsv, sorted_output_tsv) -> None:

	# skip header
	skip_header_cmd = f"""awk '/^#/' {output_tsv} > {sorted_output_tsv}"""
	skip_header_output = subprocess.check_output(skip_header_cmd, shell=True)

	# sort the rest numerically by position
	sort_cmd = f"""awk '!/^#/' {output_tsv} | sort -t$'\t' -k2,2n >> {sorted_output_tsv}"""
	sort_output = subprocess.check_output(sort_cmd, shell=True)
	print(f"sort output: {sort_output}")

	# rm unsorted and check error
	rm_cmd = f"""[ "$(wc -l < {output_tsv})" -eq "$(wc -l < {sorted_output_tsv})" ] && rm {output_tsv} || echo sort error"""
	rm_output = subprocess.check_output(rm_cmd, shell=True)
	print(f"rm output: {rm_output}")
	
if __name__ == "__main__":

	# parse args
	tsv_results_path = sys.argv[1]
	fasta_path = sys.argv[2]
	vcf_results_stub = sys.argv[3]

	# read in FASTA reference
	fasta = pysam.FastaFile(fasta_path)

	# set up output
	def get_unsorted_and_sorted_paths(stub):
		return (f"{stub}.unsorted.vcf", f"{stub}.vcf")

	(unsorted_results_path, results_path) = get_unsorted_and_sorted_paths(vcf_results_stub)
	results_vcf = open(unsorted_results_path, "w")
	writer = csv.writer(results_vcf, delimiter=OUTPUT_DELIMITER, quoting=csv.QUOTE_NONE, quotechar="")

	# write VCF headers to output files
	chrom_lengths: Dict[str, int] = get_chrom_lengths(fasta)
	write_header(writer, chrom_lengths)

	# process variants
	with open(tsv_results_path, "r") as t: 
		for variant in csv.reader(t, delimiter="\t"):
			process_variant(variant, writer, fasta)

	# close output files
	results_vcf.close()

	# sort output
	sort_output(unsorted_results_path, results_path)

	# produce encoded VCFs
	run_script("encode.py", results_path, vcf_results_stub, fasta_path)
