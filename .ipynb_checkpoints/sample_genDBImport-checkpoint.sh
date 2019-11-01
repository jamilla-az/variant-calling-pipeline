#load some modules
module load jdk/1.8.0_45-fasrc01 #java

cd 04_testGenDB #make sure there is a tmp dir there and no existing my_database

~/gatk-4.1.3.0/gatk GenomicsDBImport \
--java-options "-Xmx12g" \
-V ../03_gvcfs/BERK_1_1217_S218.raw.g.vcf \
-V ../03_gvcfs/BERK_1_1218_S219.raw.g.vcf \
-V ../03_gvcfs/BERK_1_1269_S270.raw.g.vcf \
-V ../03_gvcfs/BERK_1_1270_S271.raw.g.vcf \
--genomicsdb-workspace-path my_database_totalGenome \
-L 2L \
-L 2R \
-L 3L \
-L 3R \
-L X \
-L 4 \
--tmp-dir=tmp 
