#!/bin/bash

# vcf calls from indel callers
vcf="$1"
# where to write output of extracted breakpoints

egrep -v ^# "$vcf" | awk -F'\t' '{

	OFS = FS
	chrom = $1
	pos = $2
	ref = $4
	alt = $5

	if (ref ~ /,/ || alt ~ /,/) {
		# skip multiallelic records
		next
	}

	if (alt == "<DEL_L>" || alt == "<DEL_R>" || alt == "<INS>") { 

		# already a breakpoint
		# matches Scotch (except for 1-bp dels), Pindel insertions
		print chrom, pos, alt, "NA"

	} else if (alt == "<DEL>") { 

		# Pindel deletion without allele
		# get endpoint from END tag in INFO
		
		info = $8
		split(info, infoItems, ";")
		endTag = infoItems[1]
		split(endTag, endValues, "=")
		end = endValues[2]
		delLength = end - pos

		print chrom, pos, "<DEL_L>", delLength
		print chrom, pos + delLength, "<DEL_R>", delLength

	} else if (length(ref) > length(alt)) { 

		# deletion
		delLength = length(ref) - length(alt)
		print chrom, pos, "<DEl_L>", delLength 
		print chrom, pos + delLength, "<DEL_R>", delLength

	} else if (length(ref) < length(alt)) { 

		# insertion
		insLength = length(alt) - length(ref) 
		print chrom, pos, "<INS>", insLength

	}
}'
