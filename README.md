# Pipeline_Seurat2Pyscenic_celltype

Nextflow pipeline for performing pySCENIC analysis separately on each cell type within a Seurat object. Relies on the existence of conda environment that has pySCENIC installed, (see https://pyscenic.readthedocs.io/en/latest/installation.html).

## Parameters

All essential parameters can be set in config file.

### Key files

  *	ct_file: Path to .csv file containing names of cell types to run pySCENIC on, as found in the Seurat object metatdata column described by the ctcol paramter (see below).
  * seuratobject: Path to Seurat object in .rds format.

### Paths

  * rconda: Path conda environment containing R install with dependencies for generating cell type loom files from the Seurat object.
  * psconda: Path to pySCENIC conda environment.

### Accessory files (see https://pyscenic.readthedocs.io/en/latest/installation.html - Auxilliary datasets heading for acquiring these files)

  * tffile: Path to file of transcription factors from relevant genome for grn step of workflow.
  * fthrfiles: Path with wildcard pattern matching to databases ranking the whole genome of your species of interest based on regulatory features (i.e. transcription factors) in feather format.
  * motif2tf: Path to motif to TF annotations database providing the missing link between an enriched motif and the transcription factor that binds this motif.

### Script options

  * ctcol: Name of the Seurat object metatdata column containing the cell type assignments for each cell.

### Executor 

  * Cluster job options, edit accordingly

