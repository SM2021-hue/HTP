#setting the working directory
setwd("~/Desktop/HTP_2019")
#GDCRNATools package installation
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("GDCRNATools")
library(GDCRNATools)
#Downloading the data
library(DT)
#loading RNA-seq counts data
data(rnaCounts)
#loading miRNAs counts data
data(mirCounts)
#Normalization of RNAseq data 
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
#Normalization of miRNAs data
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)
#Parseing RNAseq data
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
datatable(as.data.frame(metaMatrix.RNA[1:5,]))#, extensions = 'Scroller',
#options = list(scrollX = TRUE, deferRender = TRUE, scroller = TRUE)
#ceRNAs network analysis. lncRNAs act as competing RNAs to regulate other transcripts
#Identification of DEGs
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
datatable(as.data.frame(DEGAll), 
        options = list(scrollX = TRUE, pageLength = 5))
#All DEGs
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
#DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')

#DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
#ceRNAs network analysis of DEGs
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
#Analyses!
datatable(as.data.frame(ceOutput), 
          options = list(scrollX = TRUE, pageLength = 5))
#Export ceRNAs network to Cytoscape
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 
                      & ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
#Export edges
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
datatable(as.data.frame(edges), 
          options = list(scrollX = TRUE, pageLength = 5))
#Export nodes
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
datatable(as.data.frame(nodes), 
          options = list(scrollX = TRUE, pageLength = 5))
#Download RNA and mature miRNA expression data
project <- 'TCGA-CHOL'
rnadir <- paste(project, 'RNAseq', sep='/')
mirdir <- paste(project, 'miRNAs', sep='/')

#Downloading RNAseq data 
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)

#Downloading mature miRNA data #######
gdcRNADownload(project.id     = 'TCGA-CHOL', 
               data.type      = 'miRNAs', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = mirdir)
#Downloading clinical data
clinicaldir <- paste(project, 'Clinical', sep='/')
gdcClinicalDownload(project.id     = 'TCGA-CHOL', 
                    write.manifest = FALSE,
                    method         = 'gdc-client',
                    directory      = clinicaldir)
#Parsing RNAseq metadata
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)

#Filtering duplicated samples in RNAseq metadata
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

#Filter non-Primary Tumor and non-Solid Tissue Normal samples in RNAseq metadata 
metaMatrix.RNA <- gdcFilterSampleType(metaMatrix.RNA)
#Parse miRNAs metadata
metaMatrix.MIR <- gdcParseMetadata(project.id = 'TCGA-CHOL',
                                   data.type  = 'miRNAs', 
                                   write.meta = FALSE)
#Filter duplicated samples in miRNAs metadata 
metaMatrix.MIR <- gdcFilterDuplicate(metaMatrix.MIR)
#Filter non-Primary Tumor and non-Solid Tissue Normal samples in miRNAs metadata
metaMatrix.MIR <- gdcFilterSampleType(metaMatrix.MIR)
#merge raw counts data
#Merge RNAseq data 
rnaCounts <- gdcRNAMerge(metadata  = metaMatrix.RNA, 
                         path      = rnadir, # the folder in which the data stored
                         organized = FALSE, 
                         data.type = 'RNAseq')

#Merge miRNAs data
mirCounts <- gdcRNAMerge(metadata  = metaMatrix.MIR,
                         path      = mirdir, # the folder in which the data stored
                         organized = FALSE, 
                         data.type = 'miRNAs')
#Merge clinical data
#Merge clinical data 
clinicalDa <- gdcClinicalMerge(path = clinicaldir, key.info = TRUE)
#removing NA values
clinicalDa[1:6,5:10]
#TMM normalization and voom transformation
#Normalization of RNAseq data 
rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)

#Normalization of miRNAs data
mirExpr <- gdcVoomNormalization(counts = mirCounts, filter = FALSE)

#Differential gene expression analysis
DEGAll <- gdcDEAnalysis(counts     = rnaCounts, 
                        group      = metaMatrix.RNA$sample_type, 
                        comparison = 'PrimaryTumor-SolidTissueNormal', 
                        method     = 'limma')
data(DEGAll)
#All DEGs#
deALL <- gdcDEReport(deg = DEGAll, gene.type = 'all')
View(deALL)
#DE long-noncoding
deLNC <- gdcDEReport(deg = DEGAll, gene.type = 'long_non_coding')
View(deLNC)
#DE protein coding genes
dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
View(dePC)
#DEG visualization
#volcano plot
gdcVolcanoPlot(DEGAll)
#Barplot
gdcBarPlot(deg = deALL, angle = 45, data.type = 'RNAseq')
#Heatmap,generated based on the heatmap.2() function in gplots package which is already an add on in the GDCRNAtools.
degName = rownames(deALL)
gdcHeatmap(deg.id = degName, metadata = metaMatrix.RNA, rna.expr = rnaExpr)
# Hypergeometric test
#Hypergenometric test is performed to test whether a lncRNA and mRNA share many miRNAs significantly.
#using the gdcCEanalysisfunction for the test
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = 'starBase', 
                          pc.targets  = 'starBase', 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
#gdcCEAnalysis() can also take user-provided miRNA-mRNA and miRNA-lncRNA interaction datasets, so, i also downloaded datasets like miRNA-target interactions predicted by TargetScan, miRanda, and Diana Tools, etc. for the ceRNAs network analysis.
#load miRNA-lncRNA interactions
data(lncTarget)
#load miRNA-mRNA interactions
data(pcTarget)
pcTarget[1:3]
ceOutput <- gdcCEAnalysis(lnc         = rownames(deLNC), 
                          pc          = rownames(dePC), 
                          lnc.targets = lncTarget, 
                          pc.targets  = pcTarget, 
                          rna.expr    = rnaExpr, 
                          mir.expr    = mirExpr)
#Network visulization in Cytoscape
#lncRNA-miRNA-mRNA interactions can be reported by the gdcExportNetwork() and visualized in Cytoscape. edges should be imported as network and nodes should be imported as feature table.
ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                        ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
#This command ran well in the first time but didn't run after that. So, I couldn't save the network image 
write.table(edges, file='edges.txt', sep='\t', quote=F)
write.table(nodes, file='nodes.txt', sep='\t', quote=F)

#Functional enrichment analysis
#Gene ontology (GO) 

enrichOutput <- gdcEnrichAnalysis(gene = rownames(deALL), simplify = TRUE)
View(enrichOutput)
#Barplot
data(enrichOutput)
gdcEnrichPlot(enrichOutput, type = 'bar', category = 'GO', num.terms = 10)

