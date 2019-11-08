# metal

```
$ python metal.py -h
usage: metal.py [-h] -s SCOTCH_VCF -d DEEPVARIANT_VCF -g GATKHC_VCF -v
                VARSCAN_VCF -p PINDELL_VCF -r REF_FASTA -o OUTPUT_DIR

Process args

optional arguments:
  -h, --help            show this help message and exit
  -s SCOTCH_VCF, --scotch_vcf SCOTCH_VCF
                        Path to Scotch VCF
  -d DEEPVARIANT_VCF, --deepvariant_vcf DEEPVARIANT_VCF
                        Path to DeepVariant VCF
  -g GATKHC_VCF, --gatkhc_vcf GATKHC_VCF
                        Path to GATK HC VCF
  -v VARSCAN_VCF, --varscan_vcf VARSCAN_VCF
                        Path to Varscan VCF
  -p PINDELL_VCF, --pindell_vcf PINDELL_VCF
                        Path to Pindel-L VCF
  -r REF_FASTA, --ref_fasta REF_FASTA
                        Path to reference FASTA
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory
```

## Dependencies

Metal requires the following Python packages:
```
pysam
```
