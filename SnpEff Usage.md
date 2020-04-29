## SnpEff Usage

##### Annotate the VCF file:

I created a new database with genome annotations using the r6.28 genome release from FlyBase - I used the .fasta files and .gtf files and the following [guide](http://snpeff.sourceforge.net/SnpEff_manual.html#databases). The other database (most recent one I could find from the available ones) was r6.12 from 2016. 

<u>Running SnpEff:</u>

You have to run from the snpEff folder...and VCFs need an INFO field which to annotate (--recode-INFO-all to keep original INFO field when making a new filtered VCF with vcftools). 

`DIR_PATH='/n/debivort_lab/Jamilla_seq/final_vcfs'`

Annotate all 2L variants (filtered for quality only):

```
java -jar snpEff.jar -v dmel_r6.28 $DIR_PATH/wildFlies_2L_filtered.vcf.gz > wildFlies_2L.ann.vcf
```

Annotate just the good-quality SNPs called for all individuals:

```
java -jar snpEff.jar -v dmel_r6.28 $DIR_PATH/allChr_fullGeno.recode.vcf > allChr_fullGeno.recode.ann.vcf
```

Note about variant counts in snpEff_summary.html: https://github.com/pcingola/SnpEff/wiki/Number-of-variants-in-VCF-and-HTML-summary-do-not-match

<u>Filter using SnpSift:</u>

Filter for just variants with warnings:

```
java -jar SnpSift.jar filter "(ANN[*].ERRORS =~ 'WARNING') || (ANN[*].ERRORS =~ 'INFO')" allChr_fullGeno.recode.ann.vcf > warning_snps.vcf

grep -v '#' warning_snps.vcf | wc -l #count SNPs retained
>2271

grep -v '#' warning_snps.vcf | grep -E 'WARNING|INFO|ERROR' | wc -l #any warnings remain?
>2271
```

SNPs that hit a gene:

```
java -jar SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has '5_UTR_variant') || (ANN[*].EFFECT has '3_UTR_variant')" allChr_fullGeno.recode.ann.vcf > genic_snps.vcf

grep -v '#' genic_snps.vcf | wc -l #count SNPs retained
>37329
```

SNPs that hit an exon (coding):

```
java -jar SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'synonymous_variant')" allChr_fullGeno.recode.ann.vcf > exonic_snps.vcf

grep -v '#' exonic_snps.vcf | wc -l #count SNPs retained
>25577
```

Missense only SNPs:

```
java -jar SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant')" allChr_fullGeno.recode.ann.vcf > missense_snps.vcf

grep -v '#' missense_snps.vcf | wc -l #count SNPs retained
>11221
```

Filter for just variants with NO warnings:

```
java -jar SnpSift.jar filter -n "(ANN[*].ERRORS =~ 'WARNING') || (ANN[*].ERRORS =~ 'INFO')" allChr_fullGeno.recode.ann.vcf > anno_PASS_snps.vcf

grep -v '#' anno_PASS_snps.vcf | wc -l #count SNPs retained
>44717

grep -v '#' anno_PASS_snps.vcf | grep -E 'WARNING|INFO|ERROR' | wc -l #any warnings remain?
>0
```

SNPs that hit a gene (no warnings):

```
java -jar SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has '5_UTR_variant') || (ANN[*].EFFECT has '3_UTR_variant')" anno_PASS_snps.vcf > genic_PASS_snps.vcf

grep -v '#' genic_PASS_snps.vcf | wc -l #count SNPs retained
>35904
```

SNPs that hit an exon (coding; no warnings):

```
java -jar SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'synonymous_variant')" anno_PASS_snps.vcf > exonic_PASS_snps.vcf

grep -v '#' exonic_PASS_snps.vcf | wc -l #count SNPs retained
>24683
```

Missense only SNPs (no warnings):

```
java -jar SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant')" anno_PASS_snps.vcf > missense_PASS_snps.vcf

grep -v '#' missense_PASS_snps.vcf | wc -l #count SNPs retained
>10809
```

##### Bootstrapping with different SNP types:

```
cd /n/debivort_lab/Jamilla_seq/final_vcfs
module load Anaconda3/5.0.1-fasrc02
source activate py_env

#genic snps
NUM_SNPS=$(grep -v '#' genic_PASS_snps.vcf | wc -l)
python ~/Seq-Data/variant_processing/vcf_to_mat.py $NUM_SNPS genic_PASS_snps.vcf genic

#exonic snps
NUM_SNPS=$(grep -v '#' exonic_PASS_snps.vcf | wc -l)
python ~/Seq-Data/variant_processing/vcf_to_mat.py $NUM_SNPS exonic_PASS_snps.vcf exon

#nonsyn snps
NUM_SNPS=$(grep -v '#' missense_PASS_snps.vcf | wc -l)
python ~/Seq-Data/variant_processing/vcf_to_mat.py $NUM_SNPS missense_PASS_snps.vcf nonsyn
```

```
cd /n/scratchlfs/debivort_lab/Jamilla/05_vcf

#genic snps bootstrap
POP_FILES=($(ls -1 /n/debivort_lab/Jamilla_seq/final_vcfs/sample_names/herit_pops))
NUMFILES=${#POP_FILES[@]}
ZBNUMFILES=$(($NUMFILES - 1))
sbatch --array=0-$ZBNUMFILES ~/Seq-Data/variant_processing/calcSegSites.sbatch
```

```
#exonic snps bootstrap
POP_FILES=($(ls -1 /n/debivort_lab/Jamilla_seq/final_vcfs/sample_names/herit_pops))
NUMFILES=${#POP_FILES[@]}
ZBNUMFILES=$(($NUMFILES - 1))
sbatch --array=0-$ZBNUMFILES ~/Seq-Data/variant_processing/calcSegSites.sbatch
```

```
#nonsyn snps bootstrap
POP_FILES=($(ls -1 /n/debivort_lab/Jamilla_seq/final_vcfs/sample_names/herit_pops))
NUMFILES=${#POP_FILES[@]}
ZBNUMFILES=$(($NUMFILES - 1))
sbatch --array=0-$ZBNUMFILES ~/Seq-Data/variant_processing/calcSegSites.sbatch
```

