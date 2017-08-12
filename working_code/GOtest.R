
GOtest <- function(intGenesAnnot, allGenesAnnot = annot, fileName, outputDir = ".",
                   ontology = "BP", algorithm = "classic")
{
  # intGenesAnnot - a data.frame with rownames as gene ENSEMBL IDs (it is a data.frame 
  #    for a historical reason, but it should just plain be a list of genes)
  # allGenesAnnot - a data.frame with rownames as ensembl IDs of universe of genes
  # to be tested against (again, god knows why I made taht a data frame)
  
    library(org.Hs.eg.db)
    library(topGO)

  # thi sfunction is needed to define GOdat object below
  signif.genes.f <- function(allGenes)
  {
    allGenes ==1
  }
    
  int.genes = rownames(intGenesAnnot)
  all.genes.names = rownames(allGenesAnnot)

  all.genes = as.integer(all.genes.names %in% int.genes)
  names(all.genes) = all.genes.names

  ####### make our own mapping and exclude electronic annotation infered
    x <- org.Hs.egGO
    # Get the entrez gene identifiers that are mapped to a GO ID
    mapped_genes <- mappedkeys(x)
    # Convert to a list
    xx <- as.list(x[mapped_genes])
    gene2GOlist = lapply(xx, function(GOlist) {
      res =   sapply(GOlist, function(GO) { 
        if (!(GO$Evidence %in% c("ND","IEA","NR")) )
            GO$GOID
        else NA
      }, USE.NAMES = FALSE)
      res[!is.na(res)]  
    })
    gene2GOlist = gene2GOlist[names(gene2GOlist) %in% annot$entrez]
    names(gene2GOlist) =  
      sapply(sapply(names(gene2GOlist), function(x){
        annot[annot$entrez == x & !is.na(annot$entrez),"geneId"]
      }), function(x){x[1]})
  ###########

  # create topGO data object
  sampleGOdata <- new("topGOdata", 
    description = paste(fileName, ontology), 
    ontology = ontology,
    allGenes = all.genes,
    geneSel = signif.genes.f,
    nodeSize = 10, ### remove the terms that have less then 10 annotated genes
    annot = annFUN.gene2GO,
    gene2GO = gene2GOlist)

  # run GO enrichment test
  resultFisher <- runTest(sampleGOdata, algorithm = algorithm, statistic = "fisher")

  # get the list of significant genes per GO term
  sel.terms = usedGO(sampleGOdata)
  ann.genes = genesInTerm(sampleGOdata, sel.terms)
  signifPerGOterm =  sapply(ann.genes, function(genesInGO){
    entrezIDs =  int.genes[int.genes %in% genesInGO]
    symbol =  sapply(entrezIDs, function(entrez) {unique(annot[rownames(annot)==entrez,"name"])} )
    paste(symbol, collapse = " ") 
  }, USE.NAMES = TRUE)


  # create table of results
  allRes <- GenTable(sampleGOdata, weight = resultFisher, orderBy = "weight", 
    ranksOf = "weight", topNodes = length(score(resultFisher)), numChar = 1000 )
  # add significant genes per term
  allRes = cbind(allRes, signif.genes = signifPerGOterm[allRes[,"GO.ID"]] )
    
  #write results into text file
  file.name =   paste0(outputDir, fileName,"-", algorithm ,"-", ontology, ".txt")
  # first write the summary
  sink(file.name) 
  print(resultFisher)
  cat("\n")
  sink()
  # then append the results table
  write.table(allRes, file = file.name, append = TRUE, sep = "\t", row.names = FALSE)


  # return a list of GO terms with weigth (p-value) <0.01
  GOterms = allRes[allRes$weight < 0.01,"GO.ID" ]
  
  sapply(ann.genes[GOterms], function(genesInGO){
       int.genes[int.genes %in% genesInGO]
  }, USE.NAMES = TRUE)

# the end
}




