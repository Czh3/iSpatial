# iSpatial: infer spatial transcriptome

## infer and enhance probe-based spatial transcriptome by intergrating with scRNA-seq

### Rationale:

- For Probe-based spatial transcriptome (such as: merFISH, STARmap, seqFISH): 
	- target few genes
	- with spatial information of all cells
- For single cell RNA-seq:
	- detect all genes
	- with out spatial information

![Rationale](iSpatial_rationale.png)

- **After infer:** transcriptome-wide spatial expression
	- all genes
	- with spatial information

### Methodology:

![method](iSpatial_method.png)

### Installation:

```
# This package based on Seurat, install Seurat first.
# You can find more in https://satijalab.org/seurat/articles/install.html
install.packages('Seurat')

install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("Czh3/iSpatial")
```

### Sample Datasets:

```
library(iSpatial)

# scRNA-seq data
data(NA_scRNA)

# merFISH data
data(NA_merFISH)
```

### Usage:
For step by step usage, check the ```vignettes``` directory of the repo.

