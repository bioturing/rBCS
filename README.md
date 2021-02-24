# rBCS

An open-source package that can write BioTuring Compressed Study (`bcs`) file from a Seurat object.

`bcs` files can be loaded directly into [BBrowser](https://bioturing.com/bbrowser), a software for single cell data.

This package is compatible with R 3.x and 4.x.

## Installation

```r
require(devtools)
install_github("bioturing/rBCS", INSTALL_opts="--no-multiarch")
# Remove INSTALL_opts if you want to use on arch i386
```

## Example

```r
require(rBCS)
ExportSeuratObject(seurat.object, "/mnt/example/path.bcs")
```

Note: Avoid using `~/...` in Windows
