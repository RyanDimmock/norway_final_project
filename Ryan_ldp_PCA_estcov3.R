library(vcfR)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(StAMPP)
library(pegas)

setwd("/Users/ryandimmock/Documents/final_project/norway_data/")
vcf <- read.vcfR("4dg_NORWAY_north__s1to6_bi_best_snps_ns4_ldp.vcf") #Filtered and LD-pruned VCF file

#Transform VCF to numeric genotypes
df <- extract.gt(vcf)
df[df == "0|0"] <- 0
df[df == "0|1"] <- 1
df[df == "1|0"] <- 1
df[df == "1|1"] <- 2
df[df == "0/0"] <- 0
df[df == "0/1"] <- 1
df[df == "1/0"] <- 1
df[df == "1/1"] <- 2
df[df == "0/0/0/0"] <- 0
df[df == "0/0/0/1"] <- 1
df[df == "0/0/1/1"] <- 2
df[df == "0/1/1/1"] <- 3
df[df == "1/1/1/1"] <- 4
df[df == "0/0/0/0/0/0"] <- 0
df[df == "0/0/0/0/0/1"] <- 1
df[df == "0/0/0/0/1/1"] <- 2
df[df == "0/0/0/1/1/1"] <- 3
df[df == "0/0/1/1/1/1"] <- 4
df[df == "0/1/1/1/1/1"] <- 5
df[df == "1/1/1/1/1/1"] <- 6
df <- data.frame(apply(df,2,function(x)as.numeric(as.character(x))))

#Remove samples with > 50% missing data
mis <- apply(df,2,function(x)sum(is.na(x))/length(x))
df <- df[,mis <= 0.75]

#Calculate allele frequencies
ploidy <- apply(df,2,max,na.rm=T)
p <- apply(df,1,function(x)sum(x,na.rm=T)/sum(ploidy[!is.na(x)]))

#Removing individuals can change allele frequencies, so we make sure that maf > 0.05
df <- df[p >= 0.05 & p <= 0.95,]
p <- p[p >= 0.05 & p <= 0.95]

#Estimate a covariance matrix
n <- ncol(df)
cov <- matrix(nrow=n,ncol=n)
for(i in 1:n){
  for(j in 1:i){
    x <- mean(c(ploidy[i],ploidy[j]))
    cov[i,j] <- mean((df[,i]-x*p)*(df[,j]-x*p)/(x*p*(1-p)),na.rm=T)
    cov[j,i] <- cov[i,j]
  }	
}

# Perform PCA on the matrix
pc <- prcomp(cov, scale=T)
# xlab <- paste0("PC1 (", round(summary(pc)$importance[2]*100), "%)")
xlab <- paste0("PC2 (", round(summary(pc)$importance[5]*100), "%)")
ylab <- paste0("PC3 (", round(summary(pc)$importance[8]*100), "%)")

# Assigning ecotype according to the sample names
ecotype <- ifelse(grepl("^GUL_|^LOD_|^MEL_|^SKI_", colnames(df)), "Estuary",
                  ifelse(grepl("^ORS_|^TJE_|^BEA_|^VAG_|^BRI_", colnames(df)), "Beach",
                         ifelse(grepl("^SOR_|^VES_|^KVA_|^TRO_", colnames(df)), "Spring", NA)))

# Define pcs data frame
pcs <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],PC3=pc$x[,3],id=colnames(df),ploidy=ploidy, ecotype=ecotype)
pcs$population <- sapply(strsplit(as.character(pcs$id), "_"), "[", 1)

# Add location and unique code information
location <- c(rep("Lofoten", 6), rep("Troms", 7))
population <- c("TJE", "ORS", "GUL", "LOD", "VES", "SOR", "BEA", "VAG", "BRI", "MEL", "SKI", "TRO", "KVA")
unique_code <- c("BL1", "BL2", "EL1", "EL2", "SL1", "SL2", "BT1", "BT2", "BT3", "ET1", "ET2", "ST1", "ST2")
population_info <- data.frame(population, unique_code, location)

# Merge population information with pcs data frame
pcs <- merge(pcs, population_info, by = "population", all.x = TRUE)

# Create a new variable for the legend that includes the population and unique code
pcs$population_legend <- paste(pcs$population, "(", pcs$unique_code, ")")

# Order populations by ecotype and location
pcs <- pcs[order(pcs$ecotype, pcs$location), ]

# Update the levels of the population factor to reflect the new order
pcs$population_legend <- factor(pcs$population_legend, levels = unique(pcs$population_legend))

# Create separate color palettes for each ecotype
beach_palette <- colorRampPalette(c("salmon","gold2","red", "brown4", "darkorange3"))(length(unique(pcs$population[pcs$ecotype == "Beach"])))
estuary_palette <- colorRampPalette(c("green", "greenyellow", "forestgreen", "darkseagreen3"))(length(unique(pcs$population[pcs$ecotype == "Estuary"])))
spring_palette <- colorRampPalette(c("purple", "deeppink", "dodgerblue2", "slateblue1"))(length(unique(pcs$population[pcs$ecotype == "Spring"])))

# Combine the color palettes into a single vector
pop_colors <- c(beach_palette, estuary_palette, spring_palette)

# Assign the colors to the populations based on their order in the pcs data frame
pop_colors <- setNames(pop_colors, unique(pcs$population_legend))

# plot PCA
ggplot(pcs, aes(x = PC2, y = PC3, color = population_legend, shape = ecotype)) +
  geom_point(size = 3.5, alpha = 0.75) +
  scale_shape_manual(values = c(16, 17, 18)) +
  scale_color_manual(values = pop_colors) +
  guides(color = guide_legend(override.aes = list(shape = 15))) +
  labs(x = xlab, y = ylab, color = "Population", shape = "Ecotype") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 11, color = "black"),
        axis.ticks.length = unit(.15, "cm"),
        axis.ticks = element_line(color = "black", size = 0.5),
        axis.title = element_text(size = 12, color = "black"),
        plot.title = element_text(size = 14, color = "black", hjust = 0.5),
        legend.text = element_text(size = 11, color = "black"),
        legend.title = element_text(size = 12, color = "black"),
        legend.key = element_blank(),
        aspect.ratio = 1)

