Here’s an updated version of the `README.md` file with your original filenames incorporated for clarity and consistency. 

```markdown
# GWAS Workflow with PLINK

This guide outlines the steps to perform Genome-Wide Association Studies (GWAS) using PLINK, with commands and outputs based on the given dataset.

---

## Prerequisites

- **PLINK**: Ensure the executable is installed.
- **Input Files**: 
  - `.vcf`: Variant Call Format file (e.g., `Qatari156_filtered_pruned.vcf`).
  - `.ped` and `.map`: Required for some conversion steps.

---

## Workflow Steps

### 1. Convert VCF to Binary Format
Convert the VCF file to PLINK binary files (`.bed`, `.bim`, `.fam`).

```bash
./plink --vcf Qatari156_filtered_pruned.vcf --make-bed --out Qatari156_filtered_pruned
```

---

### 2. Quality Control (QC)

#### a) Filter SNPs by Missing Data
Exclude SNPs with high missingness:
```bash
./plink --bfile Qatari156_filtered_pruned --geno 0.05 --make-bed --out Qatari156_filtered_pruned0.05geno
```

#### b) Hardy-Weinberg Equilibrium (HWE) Test
Remove SNPs deviating from Hardy-Weinberg Equilibrium:
```bash
./plink --bfile Qatari156_filtered_pruned0.05geno --hwe 0.000001 --make-bed --out Qatari156_filtered_pruned_hwe
```

#### c) Minor Allele Frequency (MAF) Filter
Retain SNPs with sufficient allele frequency:
```bash
./plink --bfile Qatari156_filtered_pruned_hwe --maf 0.05 --make-bed --out Qatari156_filtered_pruned_maf
```

---

### 3. Principal Component Analysis (PCA)
Perform PCA to account for population structure.

```bash
./plink --bfile Qatari156_filtered_pruned_maf --pca --out Qatari156_filtered_pruned_pca
```

Key Outputs:
- `Qatari156_filtered_pruned_pca.eigenval`: Eigenvalues.
- `Qatari156_filtered_pruned_pca.eigenvec`: Eigenvectors (principal components).

Convert `.eigenvec` to a CSV for downstream analysis:
```bash
awk '{$1=$1}1' OFS=',' Qatari156_filtered_pruned_pca.eigenvec > Qatari156_filtered_pruned_pca.csv
```

---

### 4. GWAS Analysis

#### Logistic Regression
Perform association analysis using phenotypes and covariates.

```bash
./plink --bfile Qatari156_filtered_pruned_maf --assoc --pheno phenotype_new --covar covariats.txt --covar-name PC1,PC2 --out assoc_results
```

#### In R
To perform logistic regression manually:
```R
glm(phenotype ~ SNP + PC1 + PC2, family = binomial, data = dataset)
```

---

### 5. Advanced Operations

#### Extract Chromosome-Specific SNPs
For example, extract all SNPs from chromosome 1:
```bash
./plink --bfile Qatari156_filtered_pruned_maf --chr 1 --make-bed --out chr1_data
```

#### Calculate Allele Frequencies
```bash
./plink --bfile Qatari156_filtered_pruned_maf --freq --out allele_freqs
```

---

### 6. Checking Outputs and Debugging

#### File Listings
Example:
```bash
ls Qatari156_filtered_pruned_pca*
```
Output:
- `Qatari156_filtered_pruned_pca.eigenval`
- `Qatari156_filtered_pruned_pca.eigenvec`
- `Qatari156_filtered_pruned_pca.log`

#### Count SNPs and Individuals
```bash
wc -l Qatari156_filtered_pruned.bim  # Number of SNPs
wc -l Qatari156_filtered_pruned.fam  # Number of individuals
```

---

## Example Directory Structure

| File Name                               | Description                                   |
|-----------------------------------------|-----------------------------------------------|
| `Qatari156_filtered_pruned.bed`         | Binary genotype file.                         |
| `Qatari156_filtered_pruned.bim`         | SNP map file.                                 |
| `Qatari156_filtered_pruned.fam`         | Sample information file.                      |
| `Qatari156_filtered_pruned_pca.eigenval`| Eigenvalues from PCA.                         |
| `Qatari156_filtered_pruned_pca.eigenvec`| Eigenvectors (principal components).          |
| `assoc_results.assoc`                   | GWAS results.                                 |

---

## Notes and Tips

- **File Naming**: Use consistent, descriptive file names for easy tracking.
- **Adjust Filters**: Customize thresholds for `--geno`, `--hwe`, and `--maf` to suit your dataset and goals.
- **Check Logs**: Review `.log` files for potential warnings or errors.

---

## References

- [PLINK Documentation](https://www.cog-genomics.org/plink/)
- Standard guidelines for GWAS QC and population stratification correction.
```

### Key Improvements:
1. **File Names**: Your exact filenames are used, such as `Qatari156_filtered_pruned`, for familiarity.
2. **Commands**: Organized logically with real examples and their purposes clearly stated.
3. **Output Listings**: Real output files are highlighted to provide clear expectations.
4. **Clarity**: Technical steps are presented in a way that balances depth with readability.

This version is optimized for someone following along with your dataset. Let me know if you’d like further tweaks!
