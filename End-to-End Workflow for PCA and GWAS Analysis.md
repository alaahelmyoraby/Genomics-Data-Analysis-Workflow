### **1. Set Working Directory**
Ensure all required files are in the same directory and set it as the working directory in R:
```R
getwd()
setwd("") 
```

---

### **2. Perform PCA with PLINK**
Run PLINK to generate PCA eigenvalues and eigenvectors:
```bash
./plink.exe --bfile Qatari156_filtered_pruned --pca --out pc_values
```

Convert the `.eigenvec` file to a CSV format for easier manipulation:
```bash
awk '{$1=$1}1' OFS=',' pc_values.eigenvec > pca_values.csv
```

---

### **3. Process PCA Results in R**
Load the PCA results and label columns:
```R
pca_file <- read.table("pc_values.eigenvec", header = FALSE)
pca_col <- c("FID", "IID", paste0("PC", 1:(ncol(pca_file) - 2)))
colnames(pca_file) <- pca_col
```

---

### **4. Prepare Phenotype File**
Extract the first three columns (`FID`, `IID`, `Phenotype`) from the PCA file:
```R
phenotype <- pca_file[, 1:3]
write.table(phenotype, "phenotype.txt", row.names = FALSE, quote = FALSE)
```

Perform association analysis with the phenotype:
```bash
./plink.exe --bfile Qatari156_filtered_pruned --assoc --pheno phenotype.txt --out assoc_file
```

---

### **5. Identify Significant SNPs**
Set the Bonferroni correction threshold and filter significant SNPs:
```R
bonferroni_threshold <- 0.05 / 67735  #of SNP
qatar_assoc <- read.table("assoc_file.qassoc", header = TRUE)
significant_snps <- qatar_assoc[qatar_assoc$P < bonferroni_threshold, ]
significant_snps <- significant_snps[order(significant_snps$P), ]
significant_snps <- head(significant_snps)
write.csv(significant_snps, "significant_snps.txt", row.names = FALSE)
```

Extract significant SNPs using PLINK:
```bash
./plink.exe --bfile Qatari156_filtered_pruned --extract significant_snps.txt --make-bed --out significant_snps
```

---

### **6. Prepare SNP Data**
Recode the data for significant SNPs:
```bash
./plink.exe --bfile significant_snps --recode --out significant_snps
```

Read and label the `.ped` file:
```R
sig_ped <- read.table("significant_snps.ped", header = FALSE)
colnames(sig_ped) <- c("FID", "IID", "MID", "PID", "SEX", "Phenotype", paste0("SNP", 1:(ncol(sig_ped) - 6)))
```

Merge the PCA file and SNP data:
```R
sig_pe_pc <- merge(sig_ped, pca_file, by = "FID")
```

---

### **7. Generate SNP-Specific Labels**
Create concatenated SNP genotypes for visualization:
```R
sig_pe_pc$rs10466604 <- paste0(sig_pe_pc$SNP1, sig_pe_pc$SNP2)
sig_pe_pc$rs335339 <- paste0(sig_pe_pc$SNP3, sig_pe_pc$SNP4)
sig_pe_pc$rs7355960 <- paste0(sig_pe_pc$SNP5, sig_pe_pc$SNP6)
sig_pe_pc$rs16857866 <- paste0(sig_pe_pc$SNP7, sig_pe_pc$SNP8)
sig_pe_pc$rs1841575 <- paste0(sig_pe_pc$SNP9, sig_pe_pc$SNP10)
```

---

### **8. Visualize Results**
Visualize the association between a SNP and a principal component using `ggplot2`:
```R
library(ggplot2)

ggplot(sig_pe_pc, aes(x = rs1841575, y = PC1)) +
  geom_boxplot(fill = "brown") +
  theme_minimal() +
  labs(
    title = "rs1841575 vs PC1",
    y = "PC1"
  )
```

---

### **9. Compare Covariates**
Create a covariates file including the first five columns of the PCA data:
```R
covariates <- pca_file[, 1:5]
write.table(covariates, "covariates.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)
```

Run association analysis with `PC1` and `PC2` as covariates:
```bash
./plink.exe --bfile Qatari156_filtered_pruned --assoc \
--pheno phenotype.txt --covar covariates.txt --covar-name PC1,PC2 --out pc1_2_assoc_results
```

Load and filter results:
```R
covariant_pc1_2 <- read.table("pc1_2_assoc_results.qassoc", header = TRUE)
covariant_pc1_2 <- covariant_pc1_2[covariant_pc1_2$P < bonferroni_threshold, ]
covariant_pc1_2 <- covariant_pc1_2[order(covariant_pc1_2$P), ]
```

Repeat the process with `PC2` and `PC3` as covariates:
```bash
./plink.exe --bfile Qatari156_filtered_pruned --assoc \
--pheno phenotype.txt --covar covariates.txt --covar-name PC2,PC3 --out pc2_3_assoc_results
```

Filter results:
```R
covariant_pc2_3 <- read.table("pc2_3_assoc_results.qassoc", header = TRUE)
covariant_pc2_3 <- covariant_pc2_3[covariant_pc2_3$P < bonferroni_threshold, ]
covariant_pc2_3 <- covariant_pc2_3[order(covariant_pc2_3$P), ]
```

---

### **Notes**
- This workflow uses PCA components as covariates to adjust for population stratification.
- Significant SNPs are identified using a Bonferroni correction.
- Visualization steps require the `ggplot2` library in R.

### **Dependencies**
- [PLINK](https://www.cog-genomics.org/plink/2.0/)
- R with the following libraries:
  - `ggplot2`
