#!/usr/bin/env python3

import csv
import typing
from typing import Any, List, NamedTuple

# Constants
DIST_THRESHOLD = 3
DELIMITER = "\t"

class VariantReader(NamedTuple):
	current: List
	reader: Any
	should_print_current: bool

deepvariant_vcf_path = ""
gatkhc_vcf_path = ""
pindell_vcf_path = "" # should be pindel-l
varscan_vcf_path = ""
scotch_vcf_path = ""

def get_reader(vcf: Any) -> Any:
	return csv.reader(vcf, delimiter=DELIMITER)

# VariantReader, List of VariantReader
def process_variants(readers: [VariantReader]): 

	# actually just need to compare newest to rest
	for reader in readers
		if reader.is_newest:
			query = reader

	query_pos = query.current[POS_INDEX] #make sure is int

	# also want to keep track of reader with lowest pos
	# so we know which to advance later
	lowest_pos = query_pos
	reader_with_lowest_pos: VariantReader = query

	# compare
	for reader in readers:
		if reader.is_newest: # it is the query one
			continue

		other_pos = other.current[POS_INDEX]

		if abs(query_pos - other_pos) < DIST_THRESHOLD:
			query.should_print_current = True
			other.should_print_current = True

		if other_pos < lowest_pos:
			lowest_pos = other_pos
			reader_with_lowest_pos = other

	# set all the newest to false so only lowest will have is_newest true
	for reader in readers:
		reader.is_newest = False

	# iterate lowest
	if reader_with_lowest_pos.should_print_current:
		print_variant(reader_with_lowest_pos.current)
	
	reader_with_lowest_pos.should_print_current = False
	reader_with_lowest_pos.current = next(reader_with_lowest_pos.reader)
	reader_with_lowest_pos.is_newest = True

	# what if lowest post doesn't have next?	


with open(deepvariant_vcf_path) as deepvariant_vcf, \
	open(gatkhc_vcf_path) as gatkhc_vcf, \
	open(pindell_vcf_path) as pindell_vcf, \
	open(varscan_vcf_path) as varscan_vcf, \
	open(scotch_vcf_path) as scotch_vcf:

	
	deepvariant_reader = get_reader(deepvariant_vcf)
	gatkhc_reader = get_reader(gatkhc_vcf)
	pindell_reader = get_reader(pindell_vcf)
	varscan_reader = get_reader(varscan_vcf)
	scotch_reader = get_reader(scotch_vcf)

	variant_readers = [
		VariantReader(current=next(deepvariant_reader), reader=deepvariant_reader, should_print_current=False)
	]

	variants = { 	
		"deepvariant": next(deepvariant_reader),
		"gatkhc": next(gatkhc_reader),
		"pindell": next(pindell_reader),
		"varscan": next(varscan_reader),
		"scotch": next(scotch_reader)
	}
	
	while update_variants()
		process_variants(variants)		


