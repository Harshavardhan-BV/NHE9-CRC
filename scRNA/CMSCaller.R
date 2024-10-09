library(CMScaller)
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

subset_epi <- function(GSE,counts){
    # Subset the counts to only include the Epi samples
    fname = list.files('./Data/GSE', pattern = paste0('_annotation.txt.gz'))
    metdat = fread(fname, header = T)
    epi = metdat[metdat$Cell_type == 'Epithelial cells'][[1]]
    # Select only the Epi samples
    counts = counts[,colnames(counts) %in% epi]
    return(counts)
}

run_CMS <- function(GSE){
    print(GSE)
    counts = read_imputed(GSE)
    print('Counts read')
    gc(reset=T)
    counts = subset_epi(GSE,counts)
    print('subsetted')
    gc(reset=T)
    png(paste0('figures/', GSE ,"_CMShmap.png"),  width = 10, height = 10, units = "cm", res = 500)
    cms_subt = CMScaller(counts, rowNames = "symbol", RNAseq = TRUE)
    gc(reset = TRUE, verbose = FALSE)
    dev.off()
    if (any(dim(cms_subt))){
        # Save the AUCell scores as a csv
        fwrite(cms_subt, paste0("Output/",GSE,"-CMS.csv"), row.names = T)
    }
}

GSEs = c('GSE178341')

# Iterate over GSM samples and generate rds
for (GSE in GSEs){
    run_CMS(GSE)
}