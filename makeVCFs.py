#!/usr/bin/env python3
# Convert Scotch output to VCF format
# Called by doPredict.sh as
# 	python makeVCF.py [tsv results] [fasta ref] [vcf results stub]
# Writes the full results in VCF format to ${output-stub}.vcf,

# Shifts the positions of del_L, dOne, and ins down by 1

import csv
import os
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

# get chromosome lengths from fasta reference for ##contig headers
def get_chrom_lengths(fasta: Any) -> Dict[str, int]:
	return {chrom: fasta.get_reference_length(chrom) for chrom in CHROMS}

# write_header to provided csv writer
# returns number of header lines
def write_header(writer: Any, chrom_lengths: Dict[str, int]) -> int:
	
	n_header_lines = 0

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
		n_header_lines += 1

	# contig headers
	for chrom in CHROMS:
		chrom_length: int = chrom_lengths[chrom]
		writer.writerow([f"##contig=<ID={chrom},length={chrom_length}>"])
		n_header_lines += 1

	# column headers
	writer.writerow(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
	n_header_lines += 1

	return n_header_lines

# get nucleotides from FASTA at position
def get_nucs(ref: Any, chrom: str, start: int, end: int = None) -> str:
	if not end: 
		end = start + 1
	
	# pysam is 0-based (https://pysam.readthedocs.io/en/latest/glossary.html#term-region) so subtract 1
	zero_based_start: int = start - 1
	zero_based_end: int = end - 1
	return ref.fetch(chrom, zero_based_start, zero_based_end).upper()

# for encoded vcfs, pick an arbitrary distinct allele for alt
def get_alt_for_ref(ref: str) -> str:
	ref_alt_map: Dict[str, str] = {
		"A": "T",
		"T": "C", 
		"C": "G", 
		"G": "A"
	}
	if ref in ref_alt_map: 
		return ref_alt_map[ref]
	else:
		# probably an N, ambiguous base
		return "A"

# write a variant with given fields to a list of writers
def write_variant(writers: List[Any], chrom: str, pos: int, ref: str, alt: str, info: str, gt: str) -> None:
	variant_row = [chrom, pos, ID, ref, alt, QUAL, FILTER, info, FORMAT, gt]
	for writer in writers:
		writer.writerow(variant_row)

# process variant, writing to VCFs
def process_variant(variant: List[str], writers: Dict[str, Any], fasta: Any) -> None:

	# unpack fields
	[chrom, raw_pos, pred_type, length, called_by] = variant
	pos: int = int(raw_pos)
	info: str = INFO_TEMPLATE.format(called_by)
	assert pred_type in PRED_TYPES, f"Variant at {chrom}:{pos} has unexpected type {pred_type}"

	if pred_type == "<DEL_L>":
		# this is what $SCOTCH/others-22/encode/encode.sh does
		pos += 1
		# TODO: if do this should sort, because otherwise this could mess up even sorted input
	
	# get ref
	ref: str = get_nucs(fasta, chrom, pos)

	# write results to standard vcf
	results_writers = [writers["results"]]
	write_variant(results_writers, chrom, pos, ref, pred_type, info, GT)

	# write encoded results to encoded vcfs
	encode_all_writer = [writers["encode_all"]]
	encoded_alt: str = get_alt_for_ref(ref)
	writers = encode_all_writer + [writers[f"encode_{pred_type}"]]
	write_variant(writers, chrom, pos, ref, encoded_alt, info, ENCODE_GT)

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
	(unsorted_del_L_path, del_L_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_del_L")
	(unsorted_del_R_path, del_R_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_del_R")
	(unsorted_ins_path, ins_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_ins")
	(unsorted_all_path, all_path) = get_unsorted_and_sorted_paths(f"{vcf_results_stub}.encode_all")

	results_vcf = open(unsorted_results_path, "w")
	encoded_del_L_results_vcf = open(unsorted_del_L_path, "w")
	encoded_del_R_results_vcf = open(unsorted_del_R_path, "w")
	encoded_ins_results_vcf = open(unsorted_ins_path, "w")
	encoded_all_results_vcf = open(unsorted_all_path, "w")

	def writer_for_vcf(vcf: Any) -> Any:
		return csv.writer(vcf, delimiter=OUTPUT_DELIMITER, quoting=csv.QUOTE_NONE, quotechar="")
	
	writers = {
		"results": writer_for_vcf(results_vcf),
		"encode_<DEL_L>": writer_for_vcf(encoded_del_L_results_vcf),
		"encode_<DEL_R>": writer_for_vcf(encoded_del_R_results_vcf),
		"encode_<INS>": writer_for_vcf(encoded_ins_results_vcf),
		"encode_all": writer_for_vcf(encoded_all_results_vcf),
	}
	
	# write VCF headers to output files
	n_header_lines = 0
	chrom_lengths: Dict[str, int] = get_chrom_lengths(fasta)
	for _, writer in writers.items():
		n_header_lines = write_header(writer, chrom_lengths)

	# process variants
	with open(tsv_results_path, "r") as t: 
		for variant in csv.reader(t, delimiter="\t"):
			process_variant(variant, writers, fasta)

	# close output files
	for vcf in [results_vcf, encoded_del_L_results_vcf, encoded_del_R_results_vcf, encoded_ins_results_vcf, encoded_all_results_vcf]:	
		vcf.close()

	# sort output
	sort_output(unsorted_results_path, results_path)
	sort_output(unsorted_del_L_path, del_L_path)
	sort_output(unsorted_del_R_path, del_R_path)
	sort_output(unsorted_ins_path, ins_path)
	sort_output(unsorted_all_path, all_path)
	

