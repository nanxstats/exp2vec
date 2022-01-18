# exp2vec

Tissue-specific gene embeddings trained on GTEx data.

## Installing dependencies

To restore the R package dependencies used by this project:

```r
renv::restore()
```

This pipeline uses an older version of text2vec (0.5.1),
which will require compilation (and thus Rtools on Windows).
You can also install it manually by running:

```r
remotes::install_version("text2vec", version = "0.5.1")
```
