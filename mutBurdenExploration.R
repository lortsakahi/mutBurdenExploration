library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

del <- readRDS("C:/Users/Localadmin_karlohti/Downloads/mutBurdenTabs/mutBurden_del.RDS")
syn <- readRDS("C:/Users/Localadmin_karlohti/Downloads/mutBurdenTabs/mutBurden_synonymous.RDS")
ptv <- readRDS("C:/Users/Localadmin_karlohti/Downloads/mutBurdenTabs/mutBurden_ptv.RDS")

#Looking for the genes that are the least and most mutated in each of the tables: deletions=del, synonymous=syn and loss-of-function=ptv:

MinMax <- function(x) {
  max <- which(apply(x, 1, sum)==max(apply(x, 1, sum)))
  min <- which(rowSums(x)==min(rowSums(x)))

  result <- data.frame(
    Genes_with_most_mutations = rownames(x)[max],
    Genes_with_least_mutations = rownames(x)[min]
  )
  
  return(result)
}
MinMax(del)
MinMax(syn)
MinMax(ptv)

# ggplot2 package to plot the distribution of mutation burden for all genes:

del_mutation_freq <- data.frame(Mutation_Freq = rowSums(del))
ggplot(del_mutation_freq, aes(x=Mutation_Freq)) + geom_histogram(bins=50, colour="violet", fill="pink")+ xlim(0,7500) + ylim(0,2250)

syn_mutation_freq <- data.frame(Mutation_Freq=rowSums(syn))
ggplot(syn_mutation_freq, aes(x=Mutation_Freq)) + geom_histogram(bins=50, colour="black", fill="white")+ xlim(0,7500) + ylim(0,2000)

ptv_mutation_freq <- data.frame(Mutation_Freq=rowSums(ptv))
ggplot(ptv_mutation_freq, aes(x=Mutation_Freq)) + geom_histogram(bins=50, colour="green", fill="lightgreen")+ xlim(0,6000) + ylim(0,300)

# Computing the mutation burden per cell line and also plotting the distribution for all cell lines:

del_burden <- data.frame(apply(del, 2, sum))
colnames(del_burden) <- "Cell_line_burden"
ggplot(del_burden, aes(x=Cell_line_burden)) + geom_histogram(colour="orange", fill="yellow")

syn_burden <- data.frame(apply(syn, 2, sum))
colnames(syn_burden) <- "Cell_line_burden"
ggplot(syn_burden, aes(x=Cell_line_burden)) + geom_histogram(colour="lightpink", fill="white")

ptv_burden <- data.frame(apply(ptv, 2, sum))
colnames(ptv_burden) <- "Cell_line_burden"
ggplot(ptv_burden, aes(x=Cell_line_burden)) + geom_histogram(colour="blue", fill="lightblue")

# DDD genes (downloaded from: https://panelapp.genomicsengland.co.uk/panels/484/):

DDD_genes <- read_tsv("C:/Users/Localadmin_karlohti/Downloads/DDG2P.tsv")

# Next we find out the 5% of genes that seem to be mostly mutated in the DDD gene set:
# First a new data frame with mutation sum and gene code information needs to be created to use intersect:

del_mutation_freq$SYMBOL <- row.names(del)
syn_mutation_freq$SYMBOL <- row.names(syn)
ptv_mutation_freq$SYMBOL <- row.names(ptv)

#Descending order:

ord_del <- del_mutation_freq[order(del_mutation_freq$Mutation_Freq, decreasing = TRUE), ]
ord_syn <- syn_mutation_freq[order(syn_mutation_freq$Mutation_Freq, decreasing = TRUE), ]
ord_ptv <- ptv_mutation_freq[order(ptv_mutation_freq$Mutation_Freq, decreasing = TRUE), ]

# 
Intersection <- function(x) {
  top5_quantile <- quantile(x$Mutation_Freq, probs=0.95)
  top5_genes <- data.frame(SYMBOL=x[x$Mutation_Freq>top5_quantile,]$SYMBOL, drop=FALSE)
  matches <- intersect(DDD_genes$`Entity Name`, top5_genes$SYMBOL)
  return(matches)
}

common_del <- Intersection(ord_del)
common_syn <- Intersection(ord_syn)
common_ptv <- Intersection(ord_ptv)


#Normalizing the mutation burden by the gene length("cds"), so we get total number of mutations per Mb for each gene:

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Counting the gene lengths and matching DDD_genes and their lengths with Entrez_IDs and cds_ids columns:
cds_by_gene <- cdsBy(txdb, by="gene")

computeCdsLengthPerGene <- function(geneVector, cds_by_gene){
  
  df <- select(org.Hs.eg.db, keys=geneVector, columns="ENTREZID", keytype = "SYMBOL")
  df <- subset(df, !is.na(ENTREZID))
  df <- df[df$ENTREZID %in% names(cds_by_gene),]
  list_df <- split(df$ENTREZID, df$SYMBOL)
  
  
  cds_vec <- sapply(list_df, function(x){
    sum(width(cds_by_gene[x]))
  }, simplify=T)
  
  names(cds_vec) <- gsub("\\..+","",names(cds_vec))
  
  #link names of the "cds_vec" to the data frame
  
  df$cdsLength_bp <- cds_vec[match(df$SYMBOL, names(cds_vec))]
  return(df)
  
}

bpCDSTab1 <- computeCdsLengthPerGene(rownames(del), cds_by_gene)
# bpCDSTab2 <- computeCdsLengthPerGene(rownames(syn), cds_by_gene)
# bpCDSTab3 <- computeCdsLengthPerGene(rownames(ptv), cds_by_gene)

stopifnot(all(rownames(del)==rownames(syn)))
stopifnot(all(rownames(ptv)==rownames(syn)))

normaliseAndPlot <- function(ptv, bpCDSTab){
  
  tag=substitute(ptv)
  tmp <- rowSums(ptv)
  tmpDf <- data.frame(symbolId=names(tmp),
                      mutBurden=unname(tmp))
  
  tmpDf <- merge(tmpDf, bpCDSTab, by.x="symbolId", by.y="SYMBOL")
  
  tmpDf$pseudoCounts <- ifelse(tmpDf$mutBurden==0,0.001,tmpDf$mutBurden)
  tmpDf$Normalized_Burden <- tmpDf$pseudoCounts* 1e6/tmpDf$cdsLength_bp

  ### Return object
  return(tmpDf)
}

ptvNorm <- normaliseAndPlot(ptv, bpCDSTab1)
delNorm <- normaliseAndPlot(del, bpCDSTab1)
synNorm <- normaliseAndPlot(syn, bpCDSTab1)

ggplot(delNorm, aes(x=Normalized_Burden))+ geom_histogram(bins=50,colour="darkblue", fill="magenta") + xlim(0,5e6) + ylim(0, 300)
ggplot(synNorm, aes(x=Normalized_Burden))+ geom_histogram(bins =50,colour="turquoise", fill="cyan") + xlim(0,5e6) + ylim(0, 300)
ggplot(ptvNorm, aes(x=Normalized_Burden))+ geom_histogram(bins =50,colour="coral", fill="beige") + xlim(0,3e6) + ylim(0, 150)

