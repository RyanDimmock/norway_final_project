#!/usr/bin/bash
#SBATCH --job-name=NORWAY_north_all_bi
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32g
#SBATCH --time=4-00:00:00
#SBATCH --parsable
#SBATCH --output=/gpfs01/home/mbxrd1/err_out/%x%j.out
#SBATCH --error=/gpfs01/home/mbxrd1/err_out/%x%j.err

source $HOME/.bash_profile

ml gatk-uoneasy/4.1.5.0-GCCcore-9.3.0-Java-11

vcf=/share/Yant_Group_Shared/Cochlearia_MS_datasets/400_vcf/all_bi_best_snps.vcf.gz
ref=/gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ref.gen/C_excelsa_V5.fasta

#filter into population VCFs
gatk SelectVariants -R $ref -V $vcf --restrict-alleles-to BIALLELIC -sn GUL_1 -sn GUL_10 -sn GUL_2 -sn GUL_4 -sn GUL_5 -sn GUL_6 -sn GUL_7 -sn GUL_9 \
-sn LOD_1 -sn LOD_10 -sn LOD_2 -sn LOD_4 -sn LOD_5 -sn LOD_7 -sn LOD_8 -sn LOD_9 \
-sn ORS_1 -sn ORS_2 -sn ORS_3 -sn ORS_4 -sn ORS_5 -sn ORS_6 -sn ORS_7 -sn ORS_9 \
-sn TJE_1 -sn TJE_10 -sn TJE_2 -sn TJE_5 -sn TJE_6 -sn TJE_7 -sn TJE_8 -sn TJE_9 \
-sn BEA_10 -sn BEA_2 -sn BEA_4 \
-sn VAG_102 -sn VAG_3 -sn VAG_4 -sn VAG_7 \
-sn MEL_1 -sn MEL_2 \
-sn SKI_5 -sn SKI_4 \
-sn BRI_2 -sn BRI_5 -sn BRI_6 -sn BRI_9 \
-sn VES_10 -sn VES_2 -sn VES_3 -sn VES_4 -sn VES_5 -sn VES_6 -sn VES_8 -sn VES_9 \
-sn SOR_1 -sn SOR_10 -sn SOR_2 -sn SOR_3 -sn SOR_4 -sn SOR_5 -sn SOR_6 -sn SOR_7 \
-O /gpfs01/home/mbxrd1/scan_final/north_all_bi.vcf

vcf2=/gpfs01/home/mbxrd1/scan_final/north_all_bi.vcf
ref=/gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ref.gen/C_excelsa_V5.fasta

#filter into population VCFs
gatk SelectVariants -R $ref -V $vcf2 -sn GUL_1 -sn GUL_10 -sn GUL_2 -sn GUL_4 -sn GUL_5 -sn GUL_6 -sn GUL_7 -sn GUL_9 \
-sn LOD_1 -sn LOD_10 -sn LOD_2 -sn LOD_4 -sn LOD_5 -sn LOD_7 -sn LOD_8 -sn LOD_9 \
-sn ORS_1 -sn ORS_2 -sn ORS_3 -sn ORS_4 -sn ORS_5 -sn ORS_6 -sn ORS_7 -sn ORS_9 \
-sn TJE_1 -sn TJE_10 -sn TJE_2 -sn TJE_5 -sn TJE_6 -sn TJE_7 -sn TJE_8 -sn TJE_9 \
-sn BEA_10 -sn BEA_2 -sn BEA_4 \
-sn VAG_102 -sn VAG_3 -sn VAG_4 -sn VAG_7 \
-sn MEL_1 -sn MEL_2 \
-sn SKI_5 -sn SKI_4 \
-sn BRI_2 -sn BRI_5 -sn BRI_6 -sn BRI_9 \
-O /gpfs01/home/mbxrd1/scan_final/north_nonspring_nf.vcf

gatk SelectVariants -R $ref -V $vcf2 -sn VES_10 -sn VES_2 -sn VES_3 -sn VES_4 -sn VES_5 -sn VES_6 -sn VES_8 -sn VES_9 -O /gpfs01/home/mbxrd1/scan_final/SL1_VES_nf.vcf
gatk SelectVariants -R $ref -V $vcf2 -sn SOR_1 -sn SOR_10 -sn SOR_2 -sn SOR_3 -sn SOR_4 -sn SOR_5 -sn SOR_6 -sn SOR_7 -O /gpfs01/home/mbxrd1/scan_final/SL2_SOR_nf.vcf

# do variants to table
gatk VariantsToTable -V /gpfs01/home/mbxrd1/scan_final/north_nonspring_nf.vcf -R $ref -F CHROM -F POS -F AC -F AN --output /gpfs01/home/mbxrd1/scan_final/north_nonspring_raw.table
gatk VariantsToTable -V /gpfs01/home/mbxrd1/scan_final/SL1_VES_nf.vcf -R $ref -F CHROM -F POS -F AC -F AN --output /gpfs01/home/mbxrd1/scan_final/SL1_VES_raw.table
gatk VariantsToTable -V /gpfs01/home/mbxrd1/scan_final/SL2_SOR_nf.vcf -R $ref -F CHROM -F POS -F AC -F AN --output /gpfs01/home/mbxrd1/scan_final/SL2_SOR_raw.table