#!/bin/bash
# For DeepVariant, GATK HC, VarScan2: extracts breakpoints from VCF

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
		next
	}

	if (length(ref) > length(alt)) { 
	
		# deletion
		delLength = length(ref) - length(alt)
		print chrom, pos, "<DEL_L>", delLength
		print chrom, pos + delLength, "<DEL_R>", delLength

	} else if (length(ref) < length(alt)) { 
	
		# insertion
		insLength = length(alt) - length(ref)
		print chrom, pos, "<INS>", insLength

	}

}' > "$output"
