if(!require(GEOquery)){BiocManager::install("GEOquery")}
library(GEOquery)
library(fs)

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
# https://www.nature.com/articles/nmeth.4380#article-info

gse_id <- "GSE100866"

gse <- getGEO(gse_id, GSEMatrix = TRUE)
show(gse)

filePaths <- getGEOSuppFiles(gse_id)

fs::file_delete(fs::path("data", gse_id, fs::path_file(fs::dir_ls(gse_id))))
fs::file_move(fs::path(gse_id), "data")
