library("ChIPseeker")
library("org.At.tair.db")
library("TxDb.Athaliana.BioMart.plantsmart28")
library("clusterProfiler")
library("pathview")
library("enrichplot")

args <- commandArgs(trailingOnly = T)
peak.file <- args[[1]]
chromosomes <- as.character(args[[2]])
tssup <- as.numeric(args[[3]])
tssdown <- as.numeric(args[[4]])
peak.type <- as.numeric(args[[5]])

## Defining promoter regions.
peaks <- readPeakFile(peakfile = peak.file,header= FALSE)
txdb <- TxDb.Athaliana.BioMart.plantsmart28
promoter <- getPromoters(TxDb=txdb, upstream=tssup, downstream=tssdown)

## Calculating the peak distribution along the genome.
covplot(peaks,weightCol = "V5")

## Annotating peaks according to the types of DNA regions they bind to.
peakAnno <- annotatePeak(peak = peaks, 
                         tssRegion=c(-(tssdown), tssup),
                         TxDb=txdb,annoDb = "org.At.tair.db")
plotAnnoPie(peakAnno,main="Distribution of TF binding regions",outer=T,line=-2)
plotAnnoBar(peakAnno, main="Distribution of TF binding regions")

plotDistToTSS(peakAnno,
              title="Distribution of genomic loci relative to TSS",
              ylab = "Genomic Loci (%) (5' -> 3')")

## Saving peaks that bind to promoter regions.
df.annotation <- as.data.frame(peakAnno)
{
  if (peak.type == 1)
  {
    promoter.df.annotation <- subset(df.annotation,annotation=="Promoter")
  }
  else if (peak.type == 2)
  {
    annotation.to.exclude <- c("Distal Intergenic", "Downstream (1-2kb)","Downstream (2-3kb)","Downstream (<1kb)")
    annotation.to.include <- setdiff(unique(df.annotation$annotation),annotation.to.exclude)
    promoter.df.annotation <- subset(df.annotation,df.annotation$annotation %in% annotation.to.include)
  }
}

## Listing genes affected by the TF (its regulome).
regulome <- unique(promoter.df.annotation$geneId)
write(regulome,file = "regulome.txt")

## Defining universe for GO & kegg terms enrichment.
genes.atha <- as.data.frame(genes(txdb))
{
  if (chromosomes == "all")
  {
    my.universe<-genes.atha$gene_id
  }
  else
  {
    chromosomes.for.universe <- as.vector(strsplit(chromosomes,",")[[1]])
    genes.of.my.universe <- subset(genes.atha,seqnames==chromosomes.for.universe)
    my.universe <- genes.of.my.universe$gene_id
  }
}


## GO terms BP enrichment.
ego <- enrichGO(gene          = regulome,
                universe      = my.universe,
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                keyType = "TAIR", pvalueCutoff = 0.05)
GO.enrichment <- as.data.frame(ego)
write.table(GO.enrichment,file = "go_terms.tsv",sep="\t",row.names = F)

{
  
  if(nrow(GO.enrichment) == 0)
  {
    print("No enrichment of GO terms for biological processes detected.")
  }
  
  else
  {
    pdf(file = "plots_go_bp.pdf",width = 14, height = 14)
    
    if (nrow(GO.enrichment) > 5)
    {
      b <- goplot(ego)
      plot(b)
    }
    
    c <- barplot(ego,showCategory = 30,title = "Barplot of detected GO terms for biological processes")
    plot(c)
    
    d <- dotplot(ego, showCategory=30,title = "Dotplot of detected GO terms for biological processes")
    plot(d)
    
    e <- cnetplot(ego,colorEdge=T,layout = "circle",circular=T,cex_category=2,cex_label_category=2)
    plot(e)
    
    dev.off()
  }
}

## KEGG terms enrichment.
pathway.enrich <- as.data.frame(enrichKEGG(gene=regulome, keyType = "kegg", 
                                           organism="ath", pvalueCutoff = 0.1, pAdjustMethod = "BH"))
write.table(pathway.enrich,file = "kegg_terms.tsv",sep="\t",row.names = F)
{
  if (nrow(pathway.enrich) == 0)
  {
    print("No significant enrichment of KEGG pathways was observed.")
  }
  else
  {
    for (i in 1:nrow(pathway.enrich))
    {
      pathway.current.id <- pathway.enrich[i,1]
      pathview(gene.data = regulome, pathway.id = pathway.current.id,species = "ath",
               gene.idtype = "TAIR", kegg.dir = "kegg_images/")
    }
  }
}


