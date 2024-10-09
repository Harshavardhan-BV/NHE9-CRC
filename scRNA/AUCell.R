library(Seurat)
library(AUCell)
library(GSEABase)
library(data.table)
library(dplyr)

read_imputed <- function(GSE){
    # Read the imputed counts
    counts = fread(paste0('Data/',GSE,'-magic/X.csv'), sep=',', header = F)
    counts = as.matrix(counts) %>% t(.)
    r_nam = fread(paste0('Data/',GSE,'-magic/var.csv'), sep=',', header = T)
    c_nam = fread(paste0('Data/',GSE,'-magic/obs.csv'), sep=',',header = T)
    rownames(counts) = r_nam[[1]]
    colnames(counts) = c_nam[[1]]
    return(counts)
}
run_AUCell <- function(GSE, gmt){
    print(GSE)
    # Read the counts
    matrix_AUCell = read_imputed(GSE)
    # Get the gene sets
    geneSets = getGmt(paste0('Signatures/',gmt,".gmt"))
    # Subset only the expressed genes
    geneSets <- subsetGeneSets(geneSets, rownames(matrix_AUCell)) 
    # Build cell rankings
    matrix_AUCell <- AUCell_buildRankings(matrix_AUCell, plotStats = F)
    gc(reset = TRUE, verbose = FALSE)
    # Calculate AUC metrics
    matrix_AUCell <- AUCell_calcAUC(geneSets, matrix_AUCell, nCores = parallel::detectCores()-2)
    # AUCmatrix has scores, row = pathway, col = cell
    matrix_AUCell = t(matrix_AUCell@assays@data$AUC) 
    matrix_AUCell = as.data.table(cbind(rownames(matrix_AUCell),matrix_AUCell))
    colnames(matrix_AUCell)[1] = 'Cellname'
    gc(reset = TRUE, verbose = FALSE)
    if (any(dim(matrix_AUCell))){
        # Save the AUCell scores as a csv
        fwrite(matrix_AUCell, paste0("Output/",GSE,"-",gmt,"-AUCell.csv"))
    }
}

GSEs = c('GSE132257')
gmt = 'NHE9'
# Iterate over GSM samples and generate rds
for (GSE in GSEs){
    # Select the file
    run_AUCell(GSE, gmt)
}
