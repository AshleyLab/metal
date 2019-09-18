#!/bin/bash
# For Scotch: extracts breakpoints from VCF

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
		# scotch shouldn"t have anyway
		next
	} 

	if (alt == "<DEL_L>" || alt == "<DEL_R>" || alt == "<INS>") {
		
		# already a breakpoint
		print chrom, pos, alt, "NA"

	} else if (length(ref) > length(alt)) { 
		
		# deletion
		delLength = length(ref) - length(alt)
		print chrom, pos, "<DEl_L>", delLength
		print chrom, pos + delLength, "<DEL_R>", delLength 
	
	}

}' > "$output"
