library(vcfR)

# read in vcf
vcf <- read.vcfR("4dg_NORWAY_north__s1to6_bi_best_snps.vcf")

# extract depth
dp <- extract.gt(vcf, element='DP', as.numeric=TRUE)

# Convert matrix to data frame and replace NA with 0
dp_df <- as.data.frame(dp)
dp_df[is.na(dp_df)] <- 0

# Calculate median depth for each individual and set min depth threshold
medians <- apply(dp_df, 2, median)
threshold <- 8

# Find individuals with median depth below threshold and print
low_depth_individuals <- names(medians[medians < threshold])
print(low_depth_individuals)
