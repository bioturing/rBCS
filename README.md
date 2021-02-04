# rBCS

An open-source package that can write a `SeuratObject` in BioTuring Compressed Study (`bcs`) format.
This format can be loaded directly into BBrowser.

## Example

```r
library(rBCS)
ConvertSeuratObject(x, "/mnt/example/path.bcs") # x is a Seurat object
```