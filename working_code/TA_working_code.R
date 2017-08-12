#!/usr/bin/env R

# Copyright © 2017 Tatsiana Aneichyk <taneichyk@mgh.harvard.edu>
# Distributed under terms of the MIT license.


TA_working_code <- function(){
  #  What does this function do
  #     This function does the tests 
  # Args:
  #
  # Returns:
  #

  options(stringsAsFactors = FALSE)
  workDir = "/data/talkowski/Samples/rCNVmap/"
  setwd(workDir)

  #### Create annotation with gene names, coordinates, IDs, type etc.

    # this is where we will save the file
    annotFile <- paste0(workDir, "data/RData/annot.RData") 
    if(!file.exists(annotFile)) {
      warning("annot.RData does not exist. Re-generating ...") 
      suppressPackageStartupMessages(library(biomaRt))
      # get the data from appropriate ENSEMBL archive
      ens37 <- useMart(host = "dec2013.archive.ensembl.org", 
        biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
      annot <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", 
        "start_position", "end_position", "strand", "entrezgene", "external_gene_id", 
        "gene_biotype"),  mart = ens37)
      #rename columns to our liking
      colnames(annot) = c("geneId", "chr", "start", "end", "strand", "entrez",
        "name", "type")
      # set rownames as gene IDs
      annot$geneId = as.character(annot$geneId)

      # there are several ENSEMBL IDs associated with several genes. 
      # We will collapse those into one entry
      collapsedDups =  t(sapply(unique(sapply(annot$geneId[duplicated(annot$geneId)], 
        function(id) {which(annot$geneId == id)}, USE.NAMES = TRUE)), function(ind) {
          c( paste(unique(annot$geneId[ind]), collapse = ","),
            paste(unique(annot$chr[ind]), collapse = ","),
            paste(unique(annot$start[ind]), collapse = ","),
            paste(unique(annot$end[ind]), collapse = ","),
            paste(unique(annot$strand[ind]), collapse = ","),
            paste(unique(annot$entrez[ind]), collapse = ","),
            paste(unique(annot$name[ind]), collapse = ","),
            paste(unique(annot$type[ind]), collapse = ","))
      }))
      # remove duplicates and replace entries with collapsed
      annot = annot[!duplicated(annot$geneId),]
      rownames(annot) = annot$geneId
      annot[collapsedDups[,1],] = collapsedDups

      # exclude all patches etc, leave only actual chromosomes
      annot = subset(annot, chr %in% c(1:22,"X","Y","MT"))

      save(annot, file = annotFile)
      assign("annot",annot, envir = globalenv())
    } else {
      load(annotFile, .GlobalEnv)
    }
  ##### end create annotation

  source("bin/rCNVmap/working_code/GOtest.R") 

  # get the names of files with gene lists to be tested
  listDir = paste0(workDir,"analysis/perGene_burden/signif_genes/merged/")
  listFiles  = list.files(listDir, pattern = ".+E4.+exonic.+Bonferroni.+unique", 
    recursive = TRUE) 
  genelists = sapply(listFiles, function(fileName) {
    tmp = read.table(paste0(listDir,fileName))[,1]
    rownames(subset(annot, name %in% tmp))        
  })
  names(genelists) = gsub(".genes.list","",names(genelists))      

  # get universe of genes to be tested against
  geneUniverse = read.table("/scratch/miket/rlc47temp/tmp.files/all_tested_genes.list")[,1]
  annotUniv = subset(annot, name %in% geneUniverse)

  # run the tests. THe result for each test will be save into analysis/GO_enrichments/perGene_burden/merged/
  tmp = sapply(names(genelists), function(listName){
    GOtest(intGenesAnnot = annot[genelists[[listName]],], allGenesAnnot = annotUniv, 
           fileName = listName, outputDir = "analysis/GO_enrichments/perGene_burden/merged/")
  })

}


################# Summarizing result
  # Below is reading through the results and summarizing it into a single table
  resDir = "analysis/GO_enrichments/perGene_burden/merged" 
  resFiles = list.files(resDir)
 
  # make a list of data frames with resulst
  resList = sapply(resFiles, function(fileName){
    resGO =  read.table(paste(resDir, fileName, sep = "/"), sep = "\t", skip = 10,
                      header = TRUE)
    # appears some tests had such a low p-value, that it was written as character
    # and not a number, so if thats the case, we just output top 10
    if (class(resGO$weight) == "factor"){
      return(cbind(fileName, nrow(resGO), resGO[1:10,1:6]))
    }
    # subset the result only ot the significant ones
    resGO = subset(resGO, weight < 0.01)
    # if nothing is significant, output a row on NAs
    if (nrow(resGO) == 0){
      out = rbind(resGO, rep(NA,7))
      colnames(out) = colnames(resGO)
      cbind(fileName, 0, out[,1:6])
    }
    # if more than 10 are significant, output only top 10
    if (nrow(resGO) > 10){
      cbind(fileName, nrow(resGO), resGO[1:10,1:6])
    } else { # otherwise output all
      cbind(fileName, nrow(resGO), resGO[,1:6])
    }
  }, simplify = FALSE)

  # convert list of data.frames into a single table
  resDF = do.call("rbind",resList)

  write.table(resDF, file = paste0(resDir,"/Top10.txt"), sep = "\t", row.names = FALSE)


############STRING-DB

  # just a test for plotting string-db networks through R
 library(STRINGdb)

 string_db <- STRINGdb$new( version="10", species=9606,
                              score_threshold=0, input_directory="" )


 genes = data.frame("gene" = genelists[[1]])
 exampleMapped <- string_db$map( genes, "gene", removeUnmappedRows = TRUE )

 hits = exampleMapped$STRING_id[1:400]
 pdf("STRINGdb_test.pdf")
 string_db$plot_network( hits )
 dev.off()


 
######## geneScore model
  source("bin/rCNVmap/bin/run_geneScore_model.R")
  infile <- "data/plot_data/perGene_burden/NDD_DEL_E4_exonic.geneScore_data.txt" 
  dat <- readGeneScores(infile)

  path = infile


