#!/usr/bin/env python3

import argparse
import csv
from enum import Enum
from pathlib import Path
import subprocess
from types import SimpleNamespace
import typing
from typing import Any, List

# Constants
CHROM_INDEX = 0
POS_INDEX = 1
INDEL_TYPE_INDEX = 2
DIST_THRESHOLD = 3
DELIMITER = "\t"

class Caller(Enum):
	SCOTCH = "Scotch"
	DEEPVARIANT = "DeepVariant"
	GATKHC = "GATK-HC"
	VARSCAN = "VarScan"
	PINDELL = "Pindel-L"

class IndelType(Enum): 
	DEL_L = "<DEL_L>"
	DEL_R = "<DEL_R>"
	INS = "<INS>"

# A class that wraps around a generator (reader) that yields
# indel breakpoints from a breakpoints file, and includes
class VariantReader(SimpleNamespace):
	# name of caller that called variants
	caller_name: Caller
	
	# variant the generator is currently on
	current: List

	# generator yielding indel breakpoints
	reader: Any

	# the indel type the reader is currently considering
	# important: because we advance a group of readers by just incrementing the one with the lowest position
	# we can skip over correlates because they need to be not just near in position
	# but also have the same indel type
	# so we iterate over the variants multiple times, each time pulling out one kind of indel type
	filter_for: IndelType

	# whether we've hit the end of the generator
	has_next: bool

	# list of callers who have calls that correlate with current
	correlates: [Caller]

	# whether we've printed current (don't want to print twice)
	have_printed_current: bool

	# for when this reader is in a list of readers, whether this one was the most recently advanced
	is_newest: bool

# from a variant record, return enum item representing indel type
def get_indel_type(variant: [str]) -> IndelType:
	indel_type_map = {
		"<DEL_L>": IndelType.DEL_L,
		"<DEL_R>": IndelType.DEL_R,
		"<INS>": IndelType.INS
	}
	assert variant[INDEL_TYPE_INDEX] in indel_type_map.keys(), \
		f"Variant at {variant[CHROM_INDEX]}:{variant[POS_INDEX]} has unexpected type {variant[INDEL_TYPE_INDEX]}"
	return indel_type_map[variant[INDEL_TYPE_INDEX]]

# return a VariantReader object that wraps around a list of breakpoints
def get_reader(tsv: Any, filter_for: IndelType, caller_name: Caller) -> VariantReader:

	print(tsv)
	reader = (variant for variant in csv.reader(tsv, delimiter=DELIMITER)
		if variant[INDEL_TYPE_INDEX] == filter_for.value)
	print(f"{caller_name}: {next(reader)}")

	return VariantReader(
		caller_name=caller_name,
		current=next(reader),
		reader=reader,
		filter_for=filter_for,
		has_next=True,
		correlates=[],
		have_printed_current=False,
		is_newest=False
	)

# removing duplicates (same chrom, pos, indel type) from output
# and sort numerically by position
def sort_output(output_tsv: Path, sorted_output_tsv: Path) -> None:

	# first sort considers everything up to ">" in allele, and sorts on that uniquely
	# second sort sorts all calls by position in (asc.) numeerical order
	sort_cmd = f"""sort -t">" -k1,1 -u {output_tsv} | sort -t$'\t' -k2,2n > {sorted_output_tsv}"""
	print(f"Sorting: {sort_cmd}")
	output = subprocess.check_output(sort_cmd, shell=True)
	print(f"Output: {output}")

# write current variant if has correlates
def check_current(reader: VariantReader) -> None:

	if reader.correlates and not reader.have_printed_current:
		called_in = ",".join([reader.caller_name.value] + [c.value for c in reader.correlates])
		output_writer.writerow(reader.current + [called_in])
		reader.have_printed_current = True

# run a script with given args
def run_script(script_name, *args) -> None:
	metal_dir: Path = Path(__file__).absolute().parent
	script: Path = metal_dir / script_name
	str_args = [str(a) for a in args]
	print(f"Running {script} with {str_args}")
	
	if script.suffix != ".py": 
		subprocess.call([script] + str_args)
	else:
		subprocess.call(["python", script] + str_args)

# given a list of VariantReaders, advance the one that's at the lowest position
# and still has more variants to yield
def advance_readers(readers: [VariantReader]):

	readers_with_next = [reader for reader in readers if reader.has_next]
	
	# if no readers have next, we're done
	if not readers_with_next: 
		for reader in readers:
			check_current(reader)
		return False
	
	def get_chrom_idx(c):
		if c.isdigit(): 
			return int(c)
		elif c == "X": 
			return 23
		return 24
	
	sorted_readers = sorted(readers_with_next,
		key = lambda r: (get_chrom_idx(r.current[CHROM_INDEX]), int(r.current[POS_INDEX])))
	lowest = sorted_readers[0]

	lowest_next = next(lowest.reader, None)

	if lowest_next:

		# we successfully incremented this reader, print the current variant if it had a correlate
		check_current(lowest)
		lowest.current = lowest_next
		lowest.correlates = []
		lowest.have_printed_current = False

		# set only this reader to have is_newest = True
		for reader in readers:
			reader.is_newest = False
		lowest.is_newest = True

	else:
		# this reader has no more variants to yield,
		# advance the one with the next lowest position
		lowest.has_next = False
		advance_readers(readers)

	return True

# compare current variants in each reader, 
# looking for calls that have correlates from other callers
def compare_readers(readers: [VariantReader]): 

	# compare the most recently advanced reader to the rest
	for reader in readers:
		if reader.is_newest:
			query_reader = reader
			break

	query_chrom = query_reader.current[CHROM_INDEX]
	query_pos = int(query_reader.current[POS_INDEX])
	query_filter_for = query_reader.filter_for
	query_caller_name = query_reader.caller_name

	# compare to results from other callers
	for other_reader in readers:
		if other_reader is query_reader:
			continue

		other_chrom = other_reader.current[CHROM_INDEX]
		if query_chrom != other_chrom: 
			continue

		other_pos = int(other_reader.current[POS_INDEX])
		other_caller_name = other_reader.caller_name

		# Scotch and Pindel insertions must have correlates in DeepVariant, GATK HC, or VarScan
		if (query_filter_for is IndelType.INS and
			(query_caller_name is Caller.SCOTCH or query_caller_name is Caller.PINDELL) and
			(other_caller_name is Caller.SCOTCH or other_caller_name is Caller.PINDELL)):
				continue
		

		# check if calls are within distance threshold
		# (don't need to check types match because all iterators only do one type at a time)
		if abs(query_pos - other_pos) < DIST_THRESHOLD:
			query_reader.correlates.append(other_caller_name)
			other_reader.correlates.append(query_caller_name)

# start the comparison process, looking for correlating variants from different callesr
def start_compare(readers: [VariantReader]) -> None:

	# at start, need to compare each list of variants to each other
	for reader in readers:
		for r in readers:
			r.is_newest = False
			r.correlates = []

		reader.is_newest = True
		compare_readers(readers)

		for r in readers:
			check_current(r)

	# advance through variant lists
	while advance_readers(readers):
		positions = [r.current[POS_INDEX] for r in readers]
		all_positions = " - ".join(positions)
		
		compare_readers(readers)

	for r in readers:
		check_current(r)

if __name__ ==  "__main__":

	parser = argparse.ArgumentParser(description="Process args")
	parser.add_argument("-s", "--scotch_vcf", required=True, type=str, help="Path to Scotch VCF")
	parser.add_argument("-d", "--deepvariant_vcf", required=True, type=str, help="Path to DeepVariant VCF")
	parser.add_argument("-g", "--gatkhc_vcf", required=True, type=str, help="Path to GATK HC VCF")
	parser.add_argument("-v", "--varscan_vcf", required=True, type=str, help="Path to Varscan VCF")
	parser.add_argument("-p", "--pindell_vcf", required=True, type=str, help="Path to Pindel-L VCF")
	parser.add_argument("-r", "--ref_fasta", required=True, type=str, help="Path to reference FASTA")
	parser.add_argument("-o", "--output_dir", required=True, type=str, help="Path to output directory")
	args = parser.parse_args()
	print(args)

	# get breakpoints
	ref_fasta: Path = Path(args.ref_fasta)
	output_dir: Path = Path(args.output_dir)
	assert output_dir.is_dir(), f"output_dir {args.output_dir} must be a directory that exists"

	get_breakpoints = "getBreakpoints.sh"
	breakpoints_dir: Path = output_dir / "breakpoints/"
	breakpoints_dir.mkdir(exist_ok=True) #exist ok ? 

	vcfs = {
		"scotch": args.scotch_vcf,
		"deepvariant": args.deepvariant_vcf,
		"gatkhc": args.gatkhc_vcf,
		"varscan": args.varscan_vcf,
		"pindell": args.pindell_vcf,
	}
	breakpoints_tsvs = {}

	for caller_name, vcf in vcfs.items():
		assert Path(vcf).is_file(), f"--{caller_name} must be a VCF file that exists"
		breakpoints_tsv: Path = breakpoints_dir / f"{caller_name}.breakpoints.tsv"
		
		assert not breakpoints_tsv.exists(), f"Metal writes to {breakpoints_tsv} but that already exists: please delete or move"
		run_script(get_breakpoints, vcf, breakpoints_tsv)

		breakpoints_tsvs[caller_name] = breakpoints_tsv
	
	# TODO: check metal.tsv dne also?
	output_tsv: Path = output_dir / "metal.unsorted.tsv"
	assert not output_tsv.exists(), f"Metal writes to {output_tsv} but that alredy exists: please delete or move"	
	output = open(output_tsv, "w")
	output_writer = csv.writer(output, delimiter=DELIMITER)

	# compare breakpoints
	with open(breakpoints_tsvs["scotch"]) as scotch_variants, \
		open(breakpoints_tsvs["deepvariant"]) as deepvariant_variants, \
		open(breakpoints_tsvs["gatkhc"]) as gatkhc_variants, \
		open(breakpoints_tsvs["varscan"]) as varscan_variants, \
		open(breakpoints_tsvs["pindell"]) as pindell_variants:

		for indel_type in IndelType:

			print(f"Comparing calls of type {indel_type}...")
			readers = [
				get_reader(scotch_variants, indel_type, Caller.SCOTCH),
				get_reader(deepvariant_variants, indel_type, Caller.DEEPVARIANT),
				get_reader(gatkhc_variants, indel_type, Caller.GATKHC),
				get_reader(varscan_variants, indel_type, Caller.VARSCAN), 
				get_reader(pindell_variants, indel_type, Caller.PINDELL)
			]

			start_compare(readers)
			# go back to start of file so can read again from start later
			for variants_list in [scotch_variants, deepvariant_variants, gatkhc_variants,
				varscan_variants, pindell_variants]:
				variants_list.seek(0)
	
	output.close()

	sorted_output_tsv: Path = output_dir / "metal.tsv"
	assert not sorted_output_tsv.exists(), f"Metal writes to {sorted_output_tsv} but that already exists: please delete or move"
	sort_output(output_tsv, sorted_output_tsv)

	# make VCFs
	print("Running makeVCFs.py...")
	make_vcfs = "makeVCFs.py"

	metal_output_stub: str = str((output_dir / "metal").resolve())
	run_script(make_vcfs, sorted_output_tsv, ref_fasta, metal_output_stub)

	print("Done.")
