## Installation

The repository contains R, Python, shell, and Snakemake workflows.

### 1. Clone the repository

```bash
git clone https://github.com/joweihsieh/Multiseriate-Bifacial-Vascular-Cambium-Organization.git
cd Multiseriate-Bifacial-Vascular-Cambium-Organization
```

Place the following two dependency files in the repository root if they are not already present:

```text
install_R_packages.R
requirements.txt
```

### 2. Install Conda

Install Miniconda, Anaconda, or Miniforge using the official Conda documentation:

- Conda installation: https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html

After installation, verify that Conda is available:

```bash
conda --version
```

### 3. Create and activate the analysis environment

The recorded `grid` package version was 4.3.3, which is distributed with R 4.3.3. The environment below therefore uses R 4.3.3. 

```bash
conda create \
  --name xenium-analysis \
  --channel conda-forge \
  r-base=4.3.3 \
  python=3.10 \
  pip \
  --yes

conda activate xenium-analysis
```

Confirm that the active executables come from the Conda environment:

```bash
which Rscript
Rscript --version

which python
python --version

python -m pip --version
```

On Windows, use `where Rscript` and `where python` instead of `which`.

### 4. Install the R packages

Run the supplied installer from the repository root:

```bash
Rscript install_R_packages.R
```

The script installs or verifies the following recorded versions:

```text
data.table 1.18.4
readxl 1.5.0
Matrix 1.6.5
ggplot2 4.0.3
pheatmap 1.0.13
optparse 1.8.2
scales 1.4.0
grid 4.3.3
irlba 2.3.7
RANN 2.6.2
Seurat 5.5.0
RColorBrewer 1.1.3
magrittr 2.0.5
dplyr 1.2.1
igraph 2.3.2
MASS 7.3-60.0.1
mclust 6.1.2
openxlsx 4.2.8.1
```

### 5. Install the Python packages

Install the pinned Python dependencies into the active `xenium-analysis` environment:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The supplied `requirements.txt` contains:

```text
numpy==1.24.4
pandas==2.2.1
scipy==1.9.3
matplotlib==3.9.2
scikit-learn==1.3.0
```

### 6. Install Snakemake v7.20.0

Snakemake is best kept in a separate Conda environment so that its dependencies do not alter the analysis environment:

```bash
conda create \
  --name xenium-snakemake \
  --channel conda-forge \
  --channel bioconda \
  python=3.10 \
  snakemake=7.20.0 \
  --yes

conda activate xenium-snakemake
snakemake --version
```

Official Snakemake installation documentation:

- https://snakemake.readthedocs.io/en/v7.20.0/getting_started/installation.html

Reactivate the analysis environment before running R or Python scripts directly:

```bash
conda activate xenium-analysis
```

### 7. Install Cell Ranger v7.1.0

Cell Ranger is distributed separately by 10x Genomics and is not installed by `install_R_packages.R`, `pip`, or the Conda commands above. Download Cell Ranger v7.1.0 from the official 10x Genomics previous-version page and follow the accompanying installation instructions:

- Previous Cell Ranger releases: https://www.10xgenomics.com/support/software/cell-ranger/downloads/previous-versions
- Cell Ranger installation guide: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in

After downloading the Linux archive, a typical installation is:

```bash
tar -xzf cellranger-7.1.0.tar.gz
export PATH="/path/to/cellranger-7.1.0:$PATH"
cellranger --version
```

### 8. BiocManager and Visual Studio Code

`BiocManager` is installed automatically by `install_R_packages.R` if it is absent. It is used to manage Bioconductor packages:

- https://bioconductor.org/install/

Visual Studio Code was used as a code editor but is optional and is not part of the computational runtime environment:

- https://code.visualstudio.com/

### 9. Verify the complete environment

```bash
conda activate xenium-analysis

Rscript -e 'cat(R.version.string, "\n"); print(sapply(c("data.table", "readxl", "Matrix", "ggplot2", "pheatmap", "optparse", "scales", "grid", "irlba", "RANN", "Seurat", "RColorBrewer", "magrittr", "dplyr", "igraph", "MASS", "mclust", "openxlsx"), packageVersion))'

python --version
python -m pip check
python -m pip list
```

For Snakemake and Cell Ranger:

```bash
conda run --name xenium-snakemake snakemake --version
cellranger --version
```

### 10. Expected runtime

Package installation generally requires only a few minutes in the tested environment when binaries or cached packages are available. 

Actual runtimes depend on read depth, data size, processor performance, available memory, storage speed, filesystem performance, and whether software packages must be compiled from source. The approximate analysis runtimes were:

| Workflow step | Input and resources | Approximate wall-clock time |
|---|---|---:|
| Most downstream R and Python scripts | Study datasets; standard workstation or compute node | A few minutes per script |
| Initial clustering of 5-µm Xenium bins across multiple clustering settings | Full Xenium bin-level dataset | Less than 5 hours |
| Cell Ranger preprocessing and downstream Seurat analysis | Approximately 20,000 nuclei; 48 CPU cores and 128 GB RAM | Approximately 12 hours in total |

For the reported single-nucleus RNA-seq workflow, Cell Ranger resource limits were set equivalently to:

```bash
--localcores=48 \
--localmem=128
```

## ToyData

https://www.dropbox.com/scl/fo/y740a9schxxmadct86ixp/AN1r1zqCBWoSPBmzRczyKaQ?rlkey=jbxsxdqcpuso5dzr560ou4qz4&dl=0


## Run
RUNNING_AND_EXPECTED_OUTPUTS.md
