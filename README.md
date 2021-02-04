# rBCS

An open-source package that can write a Seurat object a BioTuring Compressed Study (`bcs`) file.
This format can be loaded directly into BBrowser.

It is compatible with R 3.x and 4.x.

## Installation

```r
library(devtools)
install_github("bioturing/rBCS", INSTALL_opts = "--no-multiarch")
# Remove INSTALL_opts you want to use on arch i386
```

## Example

```r
library(rBCS)
ConvertSeuratObject(seurat.object, "/mnt/example/path.bcs")
```