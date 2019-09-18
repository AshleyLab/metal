#!/usr/bin/env python3

import argparse
import csv
from enum import Enum
from pathlib import Path
import subprocess
import typing
from typing import Any, List, NamedTuple

# Constants
DIST_THRESHOLD = 3
DELIMITER = "\t"

class VariantReader(NamedTuple):
	current: List
	reader: Any
	has_next: bool
	should_print_current: bool

class Callers(Enum):
	SCOTCH = 1
	DEEPVARIANT = 2
	GATKHC = 3
	VARSCAN = 4
	PINDELL = 5

deepvariant_vcf_path = ""
gatkhc_vcf_path = ""
pindell_vcf_path = "" # should be pindel-l
varscan_vcf_path = ""
scotch_vcf_path = ""

def print_variant(variant: List) -> None:
	print(variant)

def get_reader(vcf: Any) -> VariantReader:
	reader = csv.reader(vcf, delimiter=DELIMITER)
	return VariantReader(
		current=next(reader),
		reader=reader,
		has_next=True,
		should_print_current=False
	)

def run_script(script_name, *args):
	metal_dir: Path = Path(__file__).absolute().parent
	script: Path = metal_dir / script_name
	str_args = [str(a) for a in args]
	print(f"Running {script} with {str_args}")
	
	if script.suffix != ".py": 
		subprocess.call([script] + str_args)
	else:
		subprocess.call(["python", script] + str_args)

def advance_readers(readers: [VariantReader]):

	readers_with_next = [reader for reader in readers if reader.has_next]
	
	# if no readers have next, we're done
	if not readers_with_next: 
		for reader in readers:
			if reader.should_print_current:
				print_variant(reader.current)
		return False

	lowest = min(readers_with_next, key=(lambda r: r.current.pos))
	lowest_next = next(lowest.reader, None)

	if lowest_next:

		# print current if should
		if lowest.should_print_current: 
			print_variant(lowest.current)
			lowest.should_print_current = False
		lowest.current = lowest_next
		lowest.should_print_current = False

		# set all the other is_newest to false
		for reader in readers:
			reader.is_newest = False
		lowest.is_newest = True

	else:
		# at the end of this file, advance the next loewst
		lowest.has_next = False
		advance_readers(readers)

	return True

# VariantReader, List of VariantReader
def compare_readers(readers: [VariantReader]): 

	# actually just need to compare newest to rest
	for reader in readers:
		if reader.is_newest:
			newest_reader = reader

	newest_pos = newest_reader.current[POS_INDEX] #make sure is int

	# compare to others
	for other_reader in readers:
		if other_reader.is_newest: 
			continue

		other_pos = other_reader.current[POS_INDEX]
		if abs(newest_pos - other_pos) < DIST_THRESHOLD:
			newest_reader.should_print_current = True
			other_reader.should_print_current = True

def start_compare(readers: [VariantReader]) -> None:

	# at start, need to compare each list of variants to each other
	for i in range(len(readers)):
		for r in readers:
			r.is_newest = False
			r.should_print_current = False
		readers[i].is_newest = True
		compare_readers(readers)
		for r in readers:
			if r.should_print_current:
				print_variant(r.current)

	# leave the last one with is_newest = True
	for r in readers:
		r.should_print_current = False

	# advance through variant lists
	while True:
		compare_readers(readers)
		has_more = advance_readers(readers)
		if not has_more:
			break

	for r in readers:
		if r.should_print_current:
			print_variant(r.current)


if __name__ ==  "__main__":

	parser = argparse.ArgumentParser(description="Process args")
	parser.add_argument("-s", "--scotch_vcf", required=True, type=str, help="Path to Scotch VCF")
	parser.add_argument("-d", "--deepvariant_vcf", required=True, type=str, help="Path to DeepVariant VCF")
	parser.add_argument("-g", "--gatkhc_vcf", required=True, type=str, help="Path to GATK HC VCF")
	parser.add_argument("-v", "--varscan_vcf", required=True, type=str, help="Path to Varscan VCF")
	parser.add_argument("-p", "--pindell_vcf", required=True, type=str, help="Path to Pindel-L VCF")
	parser.add_argument("-o", "--output_dir", required=True, type=str, help="Path to output directory")
	args = parser.parse_args()
	print(args)

	# get breakpoints
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

	for caller_name, vcf in vcfs.items():
		assert Path(vcf).is_file(), f"--{caller_name} must be a VCF file that exists"
		breakpoints_tsv: Path = breakpoints_dir / f"{caller_name}.breakpoints.tsv"
		# TODO: check dne so not overwriting?
		run_script(get_breakpoints, vcf, breakpoints_tsv)


	with open(args.scotch_vcf) as scotch_variants, \
		open(args.deepvariant_vcf) as deepvariant_variants, \
		open(args.gatkhc_vcf) as gatkhc_variants, \
		open(args.varscan_vcf) as varscan_variants, \
		open(args.pindell_vcf) as pindell_variants:

		variants = (scotch_variants, deepvariant_variants, gatkhc_variants, varscan_variants, pindell_variants)
		#readers = [get_reader(v) for v in variants]      

		print("starting compare")
		#start_compare(readers)


	print("Done.")
