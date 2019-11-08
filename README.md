# Metal

Metal is a meta indel caller that composes calls from [Scotch](https://github.com/AshleyLab/scotch), DeepVariant, GATK HaplotypeCaller, Pindel (run with the `-l` option for best results), and VarScan2. 

## Installation 

Clone this repository. 

```
git clone https://github.com/AshleyLab/metal.git
```

## Dependencies

Metal requires the following Python packages:
```
pysam
```

## Run

### Overview

After running Scotch, DeepVariant, GATK HaploypeCaller, VarScan2, and Pindel, run Metal like this. 

```

python metal.py -s $scotch_vcf -d $deepvariant_vcf -g $gatkhc_vcf \
	-v $varscan_vcf -p $pindell_vcf -r $ref_fasta -o $output_dir

```
### Input

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
### Output

Metal produces several output files. `metal.vcf` includes all the results in VCF format, with alternate alleles represented by `<DEL_L>`,`<DEL_R>` or `<INS>` representing a deletion start, deletion end or insertion breakpoint, respectively. 

#### Encoding

Metal also produced several VCF files where indel breakpoints are _encoded_ as regular variants. The motivation is that some tools (including Scotch and Pindel) do not report the nucleotide sequence of alternate alleles for all variants. (They may, for example, report just `<INS>` instead.) As a result, output VCFs that include their calls, like `metal.vcf`, may not be recognized as conforming to valid VCFs. 

`encode.py` translates each indel breakpoint into an SNV at the same locus with an arbitary alternate allele, producing strictly valid VCF output. `[stub].encode_del_L.vcf` includes deletion start breakpoints represented this way, `[stub].encode_del_R.vcf` includes deletion end breakpoints, `[stub].encode_ins.vcf` includes insertion breakpoints, and `[stub].encode_all.vcf` includes all breakpoints. 

Since this process preserves breakpoint position, these files can be input to benchmarking tools like GA4GH Benchmarking that execute a distance-based comparison to evlauate tools' performance. Truth VCFs and the VCFs output by other callers to be benchmarked should also be encoded by encode.py. The script is called as

```
python encode.py input.vcf output_stub reference.fa
```
