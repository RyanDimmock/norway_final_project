# norway_final_project

## extract_DP.R

Identifies and prints out the individuals who have a median sequencing depth below a specified threshold. The script first extracts the depth from a given VCF file and converts the resulting matrix into a data frame, replacing missing values with zeros. It then calculates the median depth for each individual and sets a depth threshold (set at 4 in the above script). Lastly, the script identifies individuals whose median depth falls below this threshold and prints their names. The resulting list represents the individuals with relatively low sequencing depth in the dataset.

## Ryan_ldp_PCA_estcov3.R
PCA plot:
- Takes the LD pruned vcf data as an input file. 
- Transforms the data into numeric genotypes for handling of polyploid data. 
- Sets threshold for missing data and minor allele frequency
- Estimates a covariance matrix
- Assigns ecotype and population data to the dataframe such that these data can be highlighted in the plot
- Plots the specified principal components using ggplot2

## final_sv_v2t_4scan.sh
Preparing cohorts for selection scan using GATK:
- Selects all individuals that will make up the three cohorts used in the scan restricting to biallelic.
- Then splits this vcf into the three separate cohorts using SelectVariants. This ensures that the resulting vcf files will have same number of rows and can be compared correctly in the scan.
- Uses variants to table to convert the vcf file into table format.
