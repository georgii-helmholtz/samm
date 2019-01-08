# The systems architecture of molecular memory in poplar after abiotic stress (samm)

# Sample sheet and raw RNA-seq data can be downloaded at http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-6121
# The following code assumes that these files are stored in the subfolder "0_data" of the res.dir folder 

# Prerequisites
# Installations of tophat2 and stringtie (see https://ccb.jhu.edu/software/tophat/index.shtml and http://ccb.jhu.edu/software/stringtie/#install)
# Installations of the following R packages:
install.packages("hash")
install.packages("car", dependencies = TRUE)
install.packages("mixOmics", dependencies = TRUE)
install.packages("venn", dependencies = TRUE)
install.packages("gplots", dependencies = TRUE)
install.packages("WGCNA", dependencies = TRUE)
install.packages("dynamicTreeCut", dependencies = TRUE)
install.packages("flashClust", dependencies = TRUE)
install.packages("pheatmap", dependencies = TRUE) # imports also grid, gtable
install.packages("igraph", dependencies = TRUE)
install.packages("gridExtra", dependencies = TRUE)
install.packages("ggplot2")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("impute") # needed for WGCNA
biocLite("preprocessCore")
biocLite("GENIE3")

res.dir <- "results/"
source("functions_samm.R") # contains the functions called below and helper functions

# create result directories
# dir.create(res.dir)
# dir.create(sprintf("%s0_data", res.dir)) # should already exist # store here also the files downloaded from the RNA-seq repository

# 1. Alignments
## build reference index
## reference genome file (ref) can be downloaded at https://phytozome.jgi.doe.gov/pz/portal.html
ref <- "Ptrichocarpa_210_v3.0.hardmasked"
setwd(sprintf("%s0_data", res.dir))
system(sprintf("bowtie2-build %s.fa %s", ref, ref))
setwd("../..")
## do alignments
dir.create(sprintf("%s1_alignments", res.dir))
rd <- read.delim(sprintf("%s0_data/E-MTAB-6121.sdrf.txt", res.dir), stringsAsFactors=F, check.names=F)
files <- rd[, "Comment[SUBMITTED_FILE_NAME]"]
for (i in 1:length(files)){ # recommended to run in parallel on your cluster system of choice
  doAlignment(sprintf("%s0_data/%s", res.dir, files[i]), sprintf("%s1_alignments/%s/", res.dir, extractPart(files[i], part=1, spl=".")), sprintf("%s0_data/%s", res.dir, ref)) 
}

# 2. Quantification
# infofile for samples can be downloaded from publication and saved as tab-delimited text
# genefile can be downloaded at https://phytozome.jgi.doe.gov/pz/portal.html (here annotation version 3.1 for assembly version 3.0)
runStringtie(ids=extractPart(files, part=1, spl="."), infofile=sprintf("%s0_data/SuppTableS9_Samples.txt", res.dir), indir=sprintf("%s1_alignments/", res.dir), outdir=sprintf("%s2_quantification/", res.dir), genefile=sprintf("%s0_data/Ptrichocarpa_444_v3.1.gene.gff3", res.dir))

# 3. Differential analysis
# uses script available at http://ccb.jhu.edu/software/stringtie/dl/prepDE.py (store in the same folder as these R scripts)
runDESeq2(curr.dir=sprintf("%s2_quantification/", res.dir), outdir=sprintf("%s3_diffAnalysis/", res.dir))
# overview figure differentially expressed genes
sfile <- sprintf("%s3_diffAnalysis/DEG.Rdata", res.dir)
runStressComparison(sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir), sprintf("%s3_diffAnalysis/overviewDEG.pdf", res.dir), sfile, myx=c(-700,-700), myy=c(22, 22))
runStressComparisonZoom(sfile, sprintf("%s3_diffAnalysis/recoveryzoomDEG.pdf", res.dir), myx=c(-700,-700), myy=c(22, 22))
# stress regulation relative to recovery up-regulation 
plotTransition(sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir), sprintf("%s3_diffAnalysis/barplotSvsR.pdf", res.dir), mycex=(4/3)*1.7)
# volcano plots
volcanoWrapper(sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir), sprintf("%s3_diffAnalysis/", res.dir), sprintf("%s3_volcano/", res.dir))

# 4. PCA 
# with all tissue samples
load(sprintf("%s2_quantification/TPM_filtered.Rdata", res.dir))
plotPcaTissue(log2(TPM_filtered+1), sprintf("%s4_pca_samples/", res.dir))
# with whole-tree view
load(sprintf("%s2_quantification/treedata.R", res.dir))
w <- which(apply(treedata, 1, noNA))
treedata.alltissues <- treedata[w,] # still contains all-zero genes because before data were filtered only across all tissue samples
w0 <- which(apply(treedata.alltissues,2,sum)==0)
plotPcaTreeWithEllipse(log2(treedata.alltissues[,-w0]+1), sprintf("%s4_pca_trees/", res.dir), myx=61, myy=364, xd=445, yd=250)

# 5. CCA
doCCAWithEllipse(resdir=sprintf("%s5_cca/", res.dir), infile1=sprintf("%s2_quantification/TPM_filtered.Rdata", res.dir), infile2=sprintf("%s0_data/gas_exchange_mean.txt", res.dir))

# 6. Venn diagrams and memory heatmap
plotVenn(resdir=sprintf("%s6_venn/", res.dir), infile=sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir))
load(sprintf("%s2_quantification/TPM_filtered.Rdata", res.dir))
createMemoryHeatmapBothStresses(genes=colnames(TPM_filtered), indir=sprintf("%s3_diffAnalysis/", res.dir), outdir=sprintf("%s6_memorymap/", res.dir), lfc.thr=1, selectionfunction=selectUnion) 

# 7. Coexpression network
out.dir <- sprintf("%s7_coexpression/", res.dir)
infofile <- sprintf("%splotinfo.Rdata", out.dir)
doCoexpressionNetwork(infile1=sprintf("%s2_quantification/TPM_filtered.Rdata", res.dir), outdir=out.dir, deepSplitParam=0, minNrTissues=5, savefile=infofile, omitGrey=T)
doCoexpressionFigure(infofile, out.dir, my.cex=1)

# 8. Between-tissue self-correlation of each gene (across the different trees)
# fetch list of candidate memory genes to overlap it later on with self-correlated genes
load(sprintf("%s6_memorymap/combined.matrix.PSorCS.Rdata", res.dir))
s1 <- apply((combined.matrix.PS>0)+0,1,sum)
w1 <- which(s1>0)
mgenes <- rownames(combined.matrix.PS)[w1]
out.dir <- sprintf("%s8_selfcorrelation/", res.dir)
if (!file.exists(out.dir)){
  dir.create(out.dir)
}
mfile <- sprintf("%smemgenes.txt", out.dir)
write.table(mgenes, file=mfile, sep="\t", col.names=F, row.names=F)
# self-correlation analysis
doCorrelationAnalysis(infile=sprintf("%s2_quantification/treedata.R", res.dir), outdir=out.dir, isfile=mfile)

# 9. Gene regulatory networks
sc <- read.delim(sprintf("%scorrmemgenes.txt", out.dir), header=F, stringsAsFactors=F) # self-correlated genes: similar response functions in different tissues
out.dir.new <- sprintf("%s9_grn/", res.dir)
# obtain candidate TFs for GRNs: P.trichocarpa gene and best matching Arabidopsis gene both have TF annotation 
# Ath_TF_list and Ptr_TF_list can be downloaded at http://planttfdb.cbi.pku.edu.cn
# Ptrichocarpa_444_v3.1.annotation_info.txt can be downloaded at https://phytozome.jgi.doe.gov/pz/portal.html
rd <- read.delim(sprintf("%s0_data/Ptrichocarpa_444_v3.1.annotation_info.txt", res.dir), header=T, stringsAsFactors=F, check.names=F)
geneCol <- 2
araCol <- 11
atf <- read.delim(sprintf("%s0_data/Ath_TF_List.txt", res.dir), check.names=F, stringsAsFactors=F)
ptf <- read.delim(sprintf("%s0_data/Ptr_TF_List.txt", res.dir), check.names=F, stringsAsFactors=F)
matched.locus <- extractPart(rd[, araCol], part=1, spl=".")
w <- which(matched.locus %in% atf[,2])
tfs <- intersect(rd[w, geneCol], ptf[,2])
# GRN learning and visualization
findNeighborsEach(goi=sc[,1], sprintf("%s2_quantification/TPM_filtered.Rdata", res.dir), tfs, out.dir.new, method="GENIE3", maxNrZeros=20, top=5, tfInput=F)
graph.file <- sprintf("%sgraphdata.Rdata", out.dir.new)
colorGRN(graph.file, out.dir.new) # node and edge cooccurrence coloring
# intersecting edges across tissues
load(graph.file)
ap <- apply(edge.occurrence, 1, sum)
w <- which(ap>1)
makeIntersectingEdgeGraph(edge.occurrence[w,], sprintf("%sintersectEdges", out.dir.new), sprintf("%s0_data/Ptr_TF_List.txt", res.dir))
# tissue-specific core graphs
xgenes <- c("Potri.001G083700", "Potri.014G103000", "Potri.009G037300", "Potri.001G092100", "Potri.T111300", "Potri.016G046400")
xcolors <- c("gray80", "gray80", "lightgoldenrod1", "lightgoldenrod1", "cornsilk", "cornsilk")
makeTissueGraphCoreIgraph(tissue.edges.list[[5]], sprintf("%sXYL_core", out.dir.new), "X", sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir), xgenes, xcolors)
makeTissueGraphCoreIgraph(tissue.edges.list[[2]], sprintf("%sLE2_core", out.dir.new), "A", sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir), c("Potri.011G123300", xgenes, "Potri.010G193000"), c("darkseagreen1", xcolors, "darkslategray1"))
# ppi graph
# TairProteinInteraction.20090527.txt can be downloaded at https://www.arabidopsis.org/download_files/Proteins/Protein_interaction_data/TairProteinInteraction.20090527.txt
# Arabidopsisinteractome_SOM_TableS4.txt can be downloaded at http://science.sciencemag.org/content/333/6042/601/tab-figures-data
# pnas.1603229113.sd04.txt can be downloaded at http://www.pnas.org/content/113/29/E4238.long
edgefile <- sprintf("%s9_grn/Ath_allTF_alledges.Rdata", res.dir)
collectPPI(annfile=sprintf("%s0_data/Ptrichocarpa_444_v3.1.annotation_info.txt", res.dir), tfFile=sprintf("%s0_data/Ath_TF_List.txt", res.dir), ppifile1=sprintf("%s0_data/TairProteinInteraction.20090527.txt", res.dir), ppifile2=sprintf("%s0_data/Arabidopsisinteractome_SOM_TableS4.txt", res.dir), ppifile3=sprintf("%s0_data/pnas.1603229113.sd04.txt", res.dir), outfile=edgefile, dotfile=sprintf("%s9_grn/Ath_allTF.dot", res.dir))
load(edgefile) # all.edges
w1 <- which(all.edges[,1]=="AT1G32640")
w2 <- which(all.edges[,2]=="AT1G32640")
w3 <- which(all.edges[,1]=="AT5G65210")
w4 <- which(all.edges[,2]=="AT5G65210")
intersect(c(all.edges[w1,2], all.edges[w2,1]), c(all.edges[w3,2], all.edges[w4,1])) # "AT1G66350"
makePPIgraph(edgefile, sprintf("%s3_diffAnalysis/decisions.Rdata", res.dir), sprintf("%sppigraph", out.dir.new))

