# Following tutorial: https://lulushang.org/SpatialPCA_Tutorial/DLPFC.html
# tests in sandbox/08

library(SpatialPCA)
library(ggplot2)
library(Seurat)
library(SingleCellExperiment)
library(optparse)
library(Matrix)

option_list <- list(
  make_option(c("-a", "--adata_dir"), type = "character", help = "path to adata Rdata files"),
  make_option(c("-w", "--workdir"), type = "character", help = "Directory to save outputs"),
  make_option(c("-c", "--copy_rdata_only"), action="store_true", default=FALSE, help = "Directory to save outputs")
)
# Create the parser object
opt_parser <- OptionParser(option_list = option_list)
# Parse the arguments
opt <- parse_args(opt_parser)

print(cat("adata_dir:", opt$adata_dir, "\n"))
print(cat("workdir:", opt$workdir, "\n"))

sample_names=c("151507", "151508", "151509", "151510", "151669", "151670", "151671" ,"151672","151673", "151674" ,"151675" ,"151676")
clusterNum=c(7,7,7,7,5,5,5,5,7,7,7,7) # each sample has different ground truth cluster number

print(getwd())

if (!file.exists(opt$workdir)) {
    dir.create(opt$workdir)
}

for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    cluster_num <- clusterNum[i]
    print(sample_name)
    print(cluster_num)

    sample_output_dir <- paste(opt$workdir , sample_name, sep="/")

    if (!file.exists(sample_output_dir)) {
        dir.create(file.path(opt$workdir, sample_name))
    }
    print(sample_output_dir)

    rdata_path <- paste(opt$adata_dir, paste(sample_name, ".RData",  sep=""), sep="/")
    print(rdata_path)

    load(rdata_path) 
    print(dim(count_sub)) # The count matrix
    print(dim(xy_coords)) # The x and y coordinates. We flipped the y axis for visualization.

    obs_table_fn <- paste(sample_output_dir, file="Spatial_PCA.obs_table.tsv", sep="/")
    print(obs_table_fn)
    write.table(KRM_manual_layers_sub, file=obs_table_fn, sep="\t") # from loading the Rdata

    if (!opt$copy_rdata_only){
        xy_coords = as.matrix(xy_coords)
        rownames(xy_coords) = colnames(count_sub) # the rownames of location should match with the colnames of count matrix
        LIBD = CreateSpatialPCAObject(
            counts=count_sub, location=xy_coords, project = "SpatialPCA",
            gene.type="spatial",
            sparkversion="spark",
            numCores_spark=5, 
            gene.number=3000, customGenelist=NULL,
            min.loctions = 20, min.features=20
        )

        LIBD = SpatialPCA_buildKernel(LIBD, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
        LIBD = SpatialPCA_EstimateLoading(LIBD,fast=FALSE,SpatialPCnum=20) 
        LIBD = SpatialPCA_SpatialPCs(LIBD, fast=FALSE)

        ## Save
        LIBD_fn <- paste(sample_output_dir, file="Spatial_PCA.LIBD_obj.RData",sep = "/")
        print(LIBD_fn)
        saveRDS(LIBD, file=LIBD_fn)

        spatialpcs_fn <- paste(sample_output_dir, file="Spatial_PCA.spatialPCs.tsv",sep = "/")
        print(spatialpcs_fn)
        write.table(LIBD@SpatialPCs, spatialpcs_fn, sep="\t")

        counts_mtx <- paste(sample_output_dir, file="Spatial_PCA.counts.mtx",sep = "/")
        print(counts_mtx)
        Matrix::writeMM(LIBD@counts, file=counts_mtx)

        counts_features_fn <- paste(sample_output_dir, file="Spatial_PCA.counts_features.tsv",sep = "/")
        print(counts_features_fn)
        writeLines(rownames(LIBD@counts), counts_features_fn)

        counts_barcodes_fn <- paste(sample_output_dir, file="Spatial_PCA.counts_barcodes.tsv",sep = "/")
        print(counts_barcodes_fn)
        writeLines(colnames(LIBD@counts), counts_barcodes_fn)

        sparse_normalized_expr <- Matrix(LIBD@normalized_expr, sparse=TRUE)

        normalized_expr_mtx <- paste(sample_output_dir, file="Spatial_PCA.normalized_expr.mtx",sep = "/")
        print(normalized_expr_mtx)
        Matrix::writeMM(sparse_normalized_expr, file=normalized_expr_mtx)

    }

}