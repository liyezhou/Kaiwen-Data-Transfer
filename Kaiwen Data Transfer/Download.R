library(GDCRNATools)
gdcCorPlot <- function(gene1, gene2, rna.expr, metadata) {
  
  samples = intersect(colnames(rna.expr), metadata$sample)
  
  lncDa=rna.expr[gene1,samples]
  pcDa=rna.expr[gene2,samples]
  sampleType=as.factor(metadata$sample_type)
  
  x <- cor.test(x=lncDa, y=pcDa, alternative = 'greater')
  
  c <- format(x$estimate, digits=3)
  p <- format(x$p.value, digits=3)
  
  corDa <- data.frame(lncDa=rna.expr[gene1,], pcDa=rna.expr[gene2,], 
                      sampleType=as.factor(metadata$sample_type))
  
  
  xpos <- (min(lncDa)+max(lncDa))/2
  ypos <- as.numeric(summary(pcDa)[6])+0.5
  
  x <- sampleType
  ggplot(corDa, aes(x=lncDa, y=pcDa)) + 
    geom_point(aes(shape=sampleType, color=sampleType)) + 
    xlab(paste(gene1,' (',GDCRNATools:::ensembl2symbolFun(gene1),')',sep='')) +
    ylab(paste(gene2,' (',GDCRNATools:::ensembl2symbolFun(gene2),')',sep='')) + 
    geom_smooth(method="lm",se=FALSE, col='darkgreen', size=0.5) + 
    scale_colour_manual(breaks = levels(sampleType), 
                        values = c('chocolate1', 'blue')) +
    ggplot2::annotate("text", x = xpos, y = ypos, 
                      label = paste('cor=', c, ', p=', p, sep=''), size = 5) +
    theme_bw()+theme(legend.title = element_blank(),
                     legend.text = element_text(size=14),
                     axis.line = element_line(colour = "black"),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_rect(colour='white'),
                     panel.background = element_blank(),
                     axis.text = element_text(size=14),
                     axis.title = element_text(size=14))
}
library(DT)

project <- 'TCGA-LUAD'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

####### Download RNAseq data #######
gdcRNADownload(project.id     = 'TCGA-LUAD', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)

####### Download mature miRNA data #######
gdcRNADownload(project.id     = 'TCGA-LUAD', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)

####### Download clinical data #######
clinicaldir <- paste(project, 'Clinical', sep='/')

gdcClinicalDownload(project.id     = 'TCGA-LUAD', 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)

####### Parse RNAseq metadata #######
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-LUAD',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

####### Filter duplicated samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata #######
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)

####### Parse miRNAs metadata #######
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-LUAD',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)

####### Filter duplicated samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)

####### Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata #######
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)


####### Merge RNAseq data #######
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'RNAseq')

####### Merge miRNAs data #######
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, # if the data are in separate folders
                         data.type = 'miRNAs')

####### Merge clinical data #######
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = FALSE)
clinicalDa[1:6,5:10]

####### Normalization of RNAseq data #######
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

####### Normalization of miRNAs data #######
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)

save.image(file = "TCGA-LUAD-workspace.RData")
