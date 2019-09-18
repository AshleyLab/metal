#!/bin/bash
# For Pindel(-L): extracts breakpoints from VCF

# vcf calls from indel callers
vcf="$1"
# where to write output of extracted breakpoints
output="$2"

egrep -v ^# "$vcf" | awk -F'\t' '{

	OFS = FS
	chrom = $1
	pos = $2
	ref = $4
	alt = $5
	
	if (ref ~ /,/ || alt ~ /,/) { 
		# skip multiallelic records
		# pindel shouldn"t have anyway
		next
	}

	if (alt == "<DEL>") {
		
		# deletion, but allele not reported
		# get length from END tag in INFO

		info = $8
		split(info, infoItems, ";")
		endTag = infoItems[1]
		split(endTag, endValues, "=")
		end = endValues[2]
		delLength = end - pos

		print chrom, pos, "<DEL_L>", delLength
		print chrom, pos + delLength, "<DEL_R>", delLength
	
	} else if (alt == "<INS>") {

		# insertion, but allele not reported
		print chrom, pos, "<INS>", "NA"

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

}' > "$output"
