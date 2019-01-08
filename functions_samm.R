extractPart <- function(x, part, spl="_") {
  # splits each element of character vector x at spl and returns one part
  s <- strsplit(x, split=spl, fixed=T)
  extractCond <- function(x){
    return(x[part])
  }
  s <- sapply(s, extractCond)
  return(s)
}

doAlignment <- function(infile, outdir, ref) {
  dir.create(outdir)
  callTophat <- sprintf("tophat2 -p 2 -N 5 --read-edit-dist 5 --max-multihits 20 -o %s %s %s> %slog.txt", outdir, ref, infile, outdir)
  print(callTophat)
  system(callTophat)  
}

mapTissue <- function(x) {
  library(hash)
  tissue.name <- c("J", "A", "P", "X", "W")
  annot <- c("LE1", "LE2", "PHL", "XYL", "ROO")
  tissuehash <- hash(keys=tissue.name, values=annot)
  getTissueAnnot <- function(x){
    return(tissuehash[[sprintf('%s', x)]])
  }
  return(sapply(x, getTissueAnnot))
}

backmapTissue <- function(x) {
  library(hash)
  tissue.name <- c("J", "A", "P", "X", "W")
  annot <- c("LE1", "LE2", "PHL", "XYL", "ROO")
  tissuehash <- hash(keys=annot, values=tissue.name)
  getTissueAnnot <- function(x){
    return(tissuehash[[sprintf('%s', x)]])
  }
  return(sapply(x, getTissueAnnot))
}

adjustIds <- function(ids, addon){
  idsTrunc <- substr(ids,1,5)
  return(paste(idsTrunc, backmapTissue(addon), sep="_"))
}

runStringtie <- function(ids, infofile, indir, outdir, genefile) {
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  rd <- read.delim(infofile, stringsAsFactors=F, check.names=F)
  m <- match(ids, rd[,"SourceID"])
  ids.new <- adjustIds(rd[m, "TreeID"], rd[m, "Tissue"])
  paths <- sprintf("%s%s/accepted_hits.bam", indir, ids)
  for (i in 1:length(paths)){
    curr.command <- sprintf("stringtie %s -o %s%s.gtf -A %s%s.tab -G %s -e", paths[i], outdir, ids.new[i], outdir, ids.new[i], genefile)
    print(curr.command)
    system(curr.command) 
  }
  
  # collect TPM data from all samples into single file
  curr.files <- sprintf("%s%s.tab", outdir, ids.new)
  all.genes <- c()
  for (i in 1:length(curr.files)){
    rd <- read.delim(curr.files[i], stringsAsFactors=F, check.names=F)
    all.genes <- union(all.genes, rd[,"Gene Name"])
  }
  TPM <- matrix(0, length(curr.files), length(all.genes))
  colnames(TPM) <- all.genes
  rownames(TPM) <- ids.new
  for (i in 1:length(curr.files)){
    rd <- read.delim(curr.files[i], stringsAsFactors=F, check.names=F)
    m <- match(rd[,"Gene Name"], colnames(TPM)) # rd does not necessarily contain all genes
    TPM[i,m] <- rd[, "TPM"]
  }
  save(TPM, file=sprintf("%sTPM.Rdata", outdir))
  w0 <- which(apply(TPM,2,sum)==0)
  TPM_filtered <- TPM[,-w0]
  save(TPM_filtered, file=sprintf("%sTPM_filtered.Rdata", outdir))
  
  # whole-tree view of data
  ## Re-organization of expression data: TPM (genes x samples from all tissues) to concatenated cross-tissue matrix (plants x genes.tissueA, genes.tissueB)
  w <- which(rownames(TPM_filtered) %in% c("1_3_R_W", "3_4_S_X")) # tissue sample not from same tree
  TPM_filtered2 <- TPM_filtered[-w,]
  plant.id <- substr(rownames(TPM_filtered2), 1, nchar(rownames(TPM_filtered2)[1])-2)
  plants <- unique(plant.id)
  tissues <- c("J", "A", "W", "P", "X")
  m <- matrix(0, length(plants), ncol(TPM_filtered2)*length(tissues))
  rownames(m) <- plants
  colnames(m) <- c(sprintf("%s.%s", colnames(TPM_filtered2), tissues[1]), sprintf("%s.%s", colnames(TPM_filtered2), tissues[2]), sprintf("%s.%s", colnames(TPM_filtered2), tissues[3]), sprintf("%s.%s", colnames(TPM_filtered2), tissues[4]), sprintf("%s.%s", colnames(TPM_filtered2), tissues[5]))
  for (i in 1:length(plants)){
    for (j in 1:length(tissues)){
      w <- which(rownames(TPM_filtered2)==sprintf("%s_%s", plants[i], tissues[j]))
      if (length(w)==0){
        m[i, ((j-1)*ncol(TPM_filtered2)+1):(j*ncol(TPM_filtered2))] <- NA
      }
      else {
        m[i, ((j-1)*ncol(TPM_filtered2)+1):(j*ncol(TPM_filtered2))] <- TPM_filtered2[w[1], ]
      }
    }
  }
  treedata <- m
  save(treedata, TPM_filtered2, file=sprintf("%streedata.R", outdir))
}

mapCond <- function(x) {
  library(hash)
  cond.name <- c("1", "2", "3", "4")
  annot <- c("CS", "EC", "AC", "PS")
  condhash <- hash(keys=cond.name, values=annot)
  getCondAnnot <- function(x){
    return(condhash[[sprintf('%s', x)]])
  }
  return(sapply(x, getCondAnnot))
}

runDESeq2 <- function(curr.dir, outdir){
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  load(sprintf("%sTPM.Rdata", curr.dir))
  curr.names <- rownames(TPM)
  tempfile <- sprintf("%stemp_samples.txt", curr.dir)
  curr.files <- sprintf("%s%s.gtf", curr.dir, curr.names)
  write.table(data.frame(n=curr.names, f=curr.files), file=tempfile, row.names=F, col.names=F, sep="\t", quote=F)
  gfile <- sprintf("%sgene_counts.Rdata", curr.dir)
  tfile <- sprintf("%stranscript_counts.Rdata", curr.dir)
  tempgfile <- gsub(".Rdata", ".csv", gfile)
  temptfile <- gsub(".Rdata", ".csv", tfile)
  command <- sprintf("./prepDE.py -i %s -g %s -t %s -l 100", tempfile, tempgfile, temptfile)
  system(command) 
  adjustAndSave <- function(infile, outfile) {
    rc <- read.csv(infile)
    Counts <- t(rc[,2:ncol(rc)])
    colnames(Counts) <- rc[,1]
    rownames(Counts) <- substr(colnames(rc)[2:ncol(rc)], start=2, stop=nchar(colnames(rc)[2]))
    save(Counts, file=outfile)
  }
  adjustAndSave(tempgfile, gfile)
  adjustAndSave(temptfile, tfile)
  unlink(c(tempfile, tempgfile, temptfile))
  
  load(gfile) # Counts
  w0 <- which(apply(Counts,2,sum)==0)
  Counts_filtered <- Counts[,-w0]
  unlink(c(gfile, tfile))
  
  library(DESeq2)
  mys <- rownames(Counts_filtered)
  ty <- sapply(extractPart(mys, 2), mapCond)
  ph <- extractPart(mys, 3)
  myp <- data.frame(id=mys, group=substr(mys, start=3, stop=nchar(mys[1])), type=ty, phase=ph, tissue=extractPart(mys, 4), typePhase=paste(ty, ph, sep="."))
  v <- DESeqDataSetFromMatrix(t(Counts_filtered), myp, ~0+group)
  ds <- DESeq(v)
  # PS vs. EC and CS vs. EC
  coef.of.interest <- c()
  short <- c()
  tissues <- unique(myp$tissue)
  phases <- c("S", "R")
  stresses <- c("PS", "CS")
  stressesOld <- c("4", "1")
  for (i in 1:length(tissues)) {
    for (j in 1:length(phases)) {
      for (k in 1:length(stresses)) {
        coef.of.interest <- c(coef.of.interest, sprintf("group%s_%s_%s-group2_%s_%s", stressesOld[k], phases[j], tissues[i], phases[j], tissues[i]))
        short <- c(short, sprintf("comparison_%s_%s_%s_vsEC", tissues[i], phases[j], stresses[k]))
      }
    }
  }
  comp1 <- extractPart(coef.of.interest, part=1, spl="-")
  comp2 <- extractPart(coef.of.interest, part=2, spl="-")
  comp1 <- gsub("group", "", comp1)
  comp2 <- gsub("group", "", comp2)
  decisions <- matrix(0, ncol(Counts_filtered), length(coef.of.interest))
  colnames(decisions) <- coef.of.interest
  rownames(decisions) <- colnames(Counts_filtered)
  for (i in 1:length(coef.of.interest)){
    comparison <- results(ds, contrast=c("group", comp1[i], comp2[i]))
    save(comparison, file=sprintf("%s%s.Rdata", outdir, short[i]))
    m <- match(rownames(decisions), rownames(comparison)) # to be sure
    pos <- which(comparison[m, "padj"]<=0.05 & comparison[m, "log2FoldChange"]>=1)
    neg <- which(comparison[m, "padj"]<=0.05 & comparison[m, "log2FoldChange"]<=-1)
    decisions[pos,i] <- 1
    decisions[neg,i] <- -1
  }
  save(decisions, file=sprintf("%sdecisions.Rdata", outdir))  
  
  # PS vs. CS
  coef.of.interest <- c()
  short <- c()
  tissues <- unique(myp$tissue)
  phases <- c("S", "R")
  stresses <- c("PS")
  stressesOld <- c("4")
  for (i in 1:length(tissues)) {
    for (j in 1:length(phases)) {
      for (k in 1:length(stresses)) {
        coef.of.interest <- c(coef.of.interest, sprintf("group%s_%s_%s-group1_%s_%s", stressesOld[k], phases[j], tissues[i], phases[j], tissues[i]))
        short <- c(short, sprintf("comparison_%s_%s_%s_vsCS", tissues[i], phases[j], stresses[k]))
      }
    }
  }
  comp1 <- extractPart(coef.of.interest, part=1, spl="-")
  comp2 <- extractPart(coef.of.interest, part=2, spl="-")
  comp1 <- gsub("group", "", comp1)
  comp2 <- gsub("group", "", comp2)
  decisions <- matrix(0, ncol(Counts_filtered), length(coef.of.interest))
  colnames(decisions) <- coef.of.interest
  rownames(decisions) <- colnames(Counts_filtered)
  for (i in 1:length(coef.of.interest)){
    comparison <- results(ds, contrast=c("group", comp1[i], comp2[i]))
    save(comparison, file=sprintf("%s%s.Rdata", outdir, short[i]))
    m <- match(rownames(decisions), rownames(comparison)) # to be sure
    pos <- which(comparison[m, "padj"]<=0.05 & comparison[m, "log2FoldChange"]>=1)
    neg <- which(comparison[m, "padj"]<=0.05 & comparison[m, "log2FoldChange"]<=-1)
    decisions[pos,i] <- 1
    decisions[neg,i] <- -1
  }
  save(decisions, file=sprintf("%sdecisionsPSvsCS.Rdata", outdir))  
  
}

mapTissueColor <- function(x) {
  library(hash)
  tissuehash <- hash(keys=c("LE1", "LE2", "PHL", "XYL", "ROO"), values=c("greenyellow", "green4", "mediumblue", "deeppink2", "gold"))
  return(tissuehash[[sprintf('%s', x)]])
}

getTissuecolor <- function(x) {
  return(mapTissueColor(mapTissue(x)))
}

plotPcaTissue <- function(expr.data.2, out.dir, i=1, j=2){
  if (!file.exists(out.dir)){
    dir.create(out.dir)
  }
  my.tissues <- extractPart(rownames(expr.data.2), 4)
  mycol <- sapply(my.tissues, getTissuecolor)
  data <- t(expr.data.2)
  my.col <- list(rep(21, length(mycol)), mycol)
  normalize <- T
  my.cex <- 1.7
  my.offset <- 3
  annot.list <- c("LE1", "LE2", "PHL", "XYL", "ROO")
  tissue.color <- sapply(annot.list, mapTissueColor)
  annot.col <- list(rep(21, length(tissue.color)), tissue.color)
  
  plotLegend <- TRUE
  x <- na.omit(data)
  pc <- prcomp(t(x), scale.=normalize)
  
  s <- pc$sdev
  v <- s^2
  su <- sum(v)
  vs <- 100*v/su
  
  maxX <- max(pc$x[,i])
  minX <- min(pc$x[,i])
  maxY <- max(pc$x[,j])
  minY <- min(pc$x[,j])
  margin1 <- c(minX, maxX, minY, maxY) 
  pdf(sprintf("%sfig_%d_%d.pdf", out.dir, i, j))
  par(mar=c(4.1,5.1,4.1,2.1))
  plot(pc$x[,c(i,j)], col="black", pch=my.col[[1]], bg=mycol, xlim=margin1[c(1,2)], ylim=margin1[c(3,4)], cex=my.cex, cex.axis=my.cex, cex.lab=my.cex, xlab=sprintf("Component %d (%.0f%% of variance)", i, vs[i]), ylab=sprintf("Component %d (%.0f%% of variance)", j, vs[j]))  # double %% to escape % in sprintf
  legend("bottomleft", annot.list, col="black", pch=annot.col[[1]], pt.bg=annot.col[[2]], pt.cex=my.cex, cex=my.cex, title="Tissue")
  dev.off()
}

noNA <- function(x){
  l <- length(which(is.na(x)))
  if (l>0){
    return(F)
  }
  else {
    return(T)
  }
}

library(hash)
phasehash <- hash(keys=c("S", "R"), values=c(22,24))
getPhaseSymbol <- function(x){
  return(phasehash[[sprintf('%s', x)]])
}
treatmenthash3 <- hash(keys=c("AC", "EC", "PS", "CS"), values=c("lightblue1", "grey89", "lightcoral","lightgoldenrod1"))
getTreatmentCol3 <- function(x){
  return(treatmenthash3[[sprintf('%s', x)]])
}

plotPcaTreeWithEllipse <- function(expr.data.2, out.dir, myx, myy, xd, yd, myLegendPos="topright", scaleFactor=1.25, mylevel=0.75, i=1, j=2) {
  if (!file.exists(out.dir)){
    dir.create(out.dir)
  }
  my.phase <- extractPart(rownames(expr.data.2), 3)
  my.treatments <- extractPart(rownames(expr.data.2), 2)
  my.treatments.new <- sapply(my.treatments, mapCond)
  
  cond.pch <- sapply(my.phase, getPhaseSymbol)
  cond.col <- sapply(my.treatments.new, getTreatmentCol3)
  my.cond <- list(cond.pch, cond.col)
  cond.legend <- rep(c("AC", "EC", "PS", "CS"),2)
  annot.cond <- paste(cond.legend, c(rep(" (S)",4), rep(" (R)",4)), sep="")
  my.annot <- list(c(rep(22,4), rep(24,4)), sapply(cond.legend, getTreatmentCol3) )
  
  data <- t(expr.data.2)
  my.col <- my.cond
  normalize <- T
  my.cex <- 1.7
  my.offset <- 3
  annot.list <- annot.cond
  annot.col <- my.annot
  plotLegend <- TRUE
  
  x <- na.omit(data)
  pc <- prcomp(t(x), scale.=normalize)
  s <- pc$sdev
  v <- s^2
  su <- sum(v)
  vs <- 100*v/su

  maxX <- max(pc$x[,i])
  minX <- min(pc$x[,i])
  maxY <- max(pc$x[,j])
  minY <- min(pc$x[,j])
  margin10percent <- scaleFactor*c(minX, maxX, minY, maxY) 
  ell.groups <- rep(0, nrow(pc$x))
  w1 <- which(my.treatments.new %in% c("PS", "CS") & my.phase %in% c("S"))
  ell.groups[w1] <- 1
  library(car)
  pdf(sprintf("%sfig_%d_%d.pdf", out.dir, i, j))
  par(mar=c(5.1,5.1,4.1,2.1))
  dataEllipse(pc$x[,i], pc$x[,j], groups=as.factor(ell.groups), group.labels=rep("",2), levels=c(mylevel), col=rep("gray77",3), fill=T, fill.alpha=0.1, xlab=sprintf("Component %d (%.0f%% of variance)", i, vs[i]), ylab=sprintf("Component %d (%.0f%% of variance)", j, vs[j]), xlim=margin10percent[c(1,2)], ylim=margin10percent[c(3,4)], plot.points=FALSE, cex.axis=my.cex, cex.lab=my.cex, draw=T, add=F, grid=F, center.pch=F) 
  points(pc$x[,c(i,j)], col="black", pch=my.col[[1]], bg=my.col[[2]], cex=1.2*my.cex, cex.axis=my.cex, cex.lab=my.cex) 
  legend(x=myx, y=myy, annot.list[1:4], col="black", pch=annot.col[[1]][1:4], pt.bg=annot.col[[2]][1:4], cex=my.cex, pt.cex=1.2*my.cex, title="Type (phase)", bty="n") 
  legend(myLegendPos, annot.list[5:8], col="black", pch=annot.col[[1]][5:8], pt.bg=annot.col[[2]][5:8], cex=my.cex, pt.cex=1.2*my.cex, title="", bty="n") 
  lines(c(myx,myx), c(myy, myy-yd))
  lines(c(myx,myx+xd), c(myy-yd, myy-yd))
  text(220,10, "Stress\nstate", cex=my.cex, pos=1)
  dev.off()
}

matchData <- function(x,y){
  ov <- intersect(rownames(x), rownames(y))
  mx <- match(ov, rownames(x))
  xnew <- x[mx, , drop=FALSE]
  my <- match(ov, rownames(y))
  ynew <- y[my, , drop=FALSE]
  return(list(xnew, ynew))
}

computeRegularizationFactor <- function(x){
  # assumes that x is standardized
  # from Strimmer et al. 2005 
  x <- as.matrix(x)
  n <- nrow(x)
  cp <- crossprod(x,x)
  s <- (1/(n-1))*cp
  w.mean <- (1/n)*cp
  s.var <- matrix(0,ncol(x),ncol(x))
  for (i in 1:n){
    s.var <- s.var+(outer(x[i,],x[i,])-w.mean)^2
  }
  s.var <- (n/((n-1)^3))*s.var
  s2 <- s^2
  denom <- sum(apply(s2,1,sum))
  denom <- denom-sum(diag(s2))
  denom <- denom+sum((diag(s)-1)^2)
  lambda <- sum(apply(s.var,1,sum))/denom
  lambda <- lambda/(1-lambda) 
  return(lambda)
}

runMixomicsCCA <- function(file.dataset1, file.dataset2, var.subset1=NA, var.subset2=NA, outprefix, nrcomp=2, my.cex=1.2) {
  prepareData <- function(datafile, varset){
    rl <- load(datafile)
    data <- get(rl)
    
    if (!is.na(varset[1])){ 
      var.match <- (colnames(data) %in% varset)
      data <- data[, var.match, drop=FALSE] 
    }
    return(data)
  }
  
  dat1 <- prepareData(file.dataset1, var.subset1)
  dat2 <- prepareData(file.dataset2, var.subset2)
  matched.dat <- matchData(dat1, dat2)
  dat1 <- scale(matched.dat[[1]])
  dat2 <- scale(matched.dat[[2]])
  library(mixOmics)
  rcca.res <- rcc(dat1, dat2, ncomp=nrcomp, lambda1=computeRegularizationFactor(dat1), lambda2=computeRegularizationFactor(dat2)) # alternatively to giving lambda values, mixOmics now also has the option method="shrinkage"
  save(rcca.res,file=sprintf("%sresults.Rdata", outprefix))
}

doCCAWithEllipse <- function(resdir, infile1, infile2, my.cex=1.7) {
  if (!file.exists(resdir)){
    dir.create(resdir)
  }
  # average within group (because samples are not paired)
  load(infile1) # TPM_filtered
  groupid <- substr(rownames(TPM_filtered), start=3, stop=nchar(rownames(TPM_filtered)[1]))
  groups <- unique(groupid)
  TPM_filtered_logavg <- matrix(0, length(groups), ncol(TPM_filtered))
  for (i in 1:length(groups)){
    w <- which(groupid==groups[i])
    TPM_filtered_logavg[i,] <- apply(log2(TPM_filtered[w,]+1), 2, mean)
  }
  colnames(TPM_filtered_logavg) <- colnames(TPM_filtered)
  rownames(TPM_filtered_logavg) <- sprintf("m_%s", groups)
  # extract mature leaf data
  my.tissues <- extractPart(rownames(TPM_filtered_logavg), 4)
  wa <- which(my.tissues=="A")
  expr.dat <- TPM_filtered_logavg[wa,]
  
  rd.mean <- read.delim2(infile2, check.names = FALSE, stringsAsFactors=F)
  
  # name mapping
  datafile1 <- sprintf("%scca_expr_red.Rdata", resdir)
  datafile2 <- sprintf("%scca_ps.Rdata", resdir)
  my.phase <- extractPart(rownames(expr.dat), 3)
  my.treatments <- extractPart(rownames(expr.dat), 2)
  my.treatments.new <- sapply(my.treatments, mapCond)
  rownames(expr.dat) <- paste(my.treatments.new, my.phase, sep=".")
  var.s <- apply(expr.dat,2,sd)
  var.s.sort <- sort(var.s, decreasing=T, index.return=T)
  expr.data.old.leaf.3 <- expr.dat[, var.s.sort$ix[1:100] ]
  save(expr.data.old.leaf.3, file=datafile1)
  
  ps.data <- as.matrix(rd.mean[,-c(1,2)])
  rownames(ps.data) <- paste(rd.mean[,1], rd.mean[,2], sep=".") # row mapping done within runMixomicsCCA
  save(ps.data, file=datafile2)
  
  runMixomicsCCA(datafile1, datafile2, outprefix=resdir)
  
  load(sprintf("%sresults.Rdata", resdir))
  
  # compute projections of replicates to plot ellipse
  TPM_filtered_log <- log2(TPM_filtered+1)
  wa2 <- which(extractPart(rownames(TPM_filtered_log), 4, spl="_")=="A")
  expr.data.old.leaf.3rep <- TPM_filtered_log[wa2, var.s.sort$ix[1:100] ]
  ssc <- scale(expr.data.old.leaf.3)
  expr.data.old.leaf.3rep.proj <- scale(expr.data.old.leaf.3rep, center=attr(ssc, "scaled:center"), scale=attr(ssc, "scaled:scale")) %*% rcca.res$loadings[[1]] 
  rep.cond <- sapply(extractPart(rownames(expr.data.old.leaf.3rep), 2, spl="_"), mapCond)
  rep.phase <- extractPart(rownames(expr.data.old.leaf.3rep), 3, spl="_")
  ell.groups <- rep(0, nrow(expr.data.old.leaf.3rep))
  w1 <- which(rep.cond %in% c("PS", "CS") & rep.phase %in% c("S"))
  ell.groups[w1] <- 1
  w2 <- which(rep.cond %in% c("PS", "CS") & rep.phase %in% c("R"))
  ell.groups[w2] <- 2
  
  # plot
  samples <- rownames(rcca.res$variates[[1]])
  my.treatments <- sapply(samples, extractPart, part=1, spl=".")
  my.phases <- sapply(samples, extractPart, part=2, spl=".")
  my.bg <- sapply(my.treatments, getTreatmentCol3)
  my.pch <- sapply(my.phases, getPhaseSymbol)
  samples.legend <- gsub(".", " (", samples, fixed=T)
  samples.legend <- paste(samples.legend, ")", sep="")
  my.ix <- c(6,4,8,2,5,3,7,1)
  
  library(car)
  pdf(sprintf("%splot_projection.pdf", resdir))
  par(mar=c(5.1,5.1,4.1,2.1))
  dataEllipse(expr.data.old.leaf.3rep.proj[,1], expr.data.old.leaf.3rep.proj[,2], groups=as.factor(ell.groups), group.labels=rep("",3), levels=c(0.75), col=rep("gray77",3), fill=T, fill.alpha=0.1, xlab=sprintf("Component 1 (correlation %.2f)", rcca.res$cor[1]), ylab=sprintf("Component 2 (correlation %.2f)", rcca.res$cor[2]), plot.points=FALSE, cex.axis=my.cex, cex.lab=my.cex, draw=T, add=F, grid=F, center.pch=F, ylim=c(-2.7,2), xaxp=c(-1,1,2))
  points(rcca.res$variates[[1]][,1], rcca.res$variates[[1]][,2], col="black", pch=my.pch, bg=my.bg, cex=1.2*my.cex)
  legend(x=-1.25, y=0.45, samples.legend[my.ix], col="black", pch=my.pch[my.ix], pt.bg=my.bg[my.ix], cex=my.cex, pt.cex=1.2*my.cex, title="Type (phase)")
  text(x=c(0.5, -1.55, 0.4), y=c(-1.1,-0.3,1.1), labels=c("Untreated state", "Stress\nstate", "Post-recovery state"), col="black", cex=my.cex, pos=3) 
  dev.off()
  
}

runVenn <- function(infile, plotfile, intersectfile, mycex=1.2, stressId, phases, direction, mytitle="") { 
  load(infile) # decisions
  tissues <- c("J", "A", "P", "X", "W")
  tissue.color <- sapply(tissues, getTissuecolor)
  intersectall <- rownames(decisions)
  tissueDE <- list(intersectall, intersectall, intersectall, intersectall, intersectall)
  for (i in 1:length(tissues)) {
    for (j in 1:length(phases)) {
      wcomp <- which(colnames(decisions)==sprintf("group%d_%s_%s-group2_%s_%s", stressId, phases[j], tissues[i], phases[j], tissues[i]))
      w <- which(decisions[, wcomp]==direction)
      intersectall <- intersect(intersectall, rownames(decisions)[w])
      tissueDE[[i]] <- intersect(tissueDE[[i]], rownames(decisions)[w])
    }
  }
  write.table(data.frame(is=intersectall), file=intersectfile, sep="\t", col.names=F, row.names=F)
  
  library(venn)
  pdf(plotfile)
  venn(x=tissueDE, snames=mapTissue(tissues), col=tissue.color, cexil=mycex, cexsn=mycex, lwd=2) 
  text(x=20, y=950, labels=mytitle, pos=4, cex=mycex, font=1)
  dev.off()
}

plotVenn <- function (resdir, infile) {
  if (!file.exists(resdir)){
    dir.create(resdir)
  }
  stressIds <- c(4,1)
  stressIdNames <- c("PS", "CS")
  phases <- c("S", "R")
  directions <- c(1,-1)
  directionNames <- c("Up", "Down")
  for (i in 1:length(stressIds)){
    for (j in 1:length(phases)){
      for (k in 1:length(directions)){
        runVenn(infile, sprintf("%s%s_%s_%s.pdf", resdir, stressIdNames[i], phases[j], directionNames[k]), sprintf("%s%s_%s_%s_intersect.txt", resdir, stressIdNames[i], phases[j], directionNames[k]), mycex=1.7, stressId=stressIds[i], phases=c(phases[j]), direction=directions[k], mytitle=sprintf("%s (%s) %s", stressIdNames[i], phases[j], directionNames[k] ))
        
      }
    }
  }
}

selectUnion <- function (X, Y) {
  s1 <- apply((X>0)+0,1,sum)
  s2 <- apply((Y>0)+0,1,sum)
  w1 <- which(s1>1)
  w2 <- which(s2>1)
  w <- union(w1, w2)
  return(cbind(X[w,], Y[w,]))
}

classifyDE <- function (X, fdr.thr=0.05, lfc.thr=1) {
  r <- rep(0, nrow(X))
  w1 <- which(X[, "padj"]<=fdr.thr & X[, "log2FoldChange"]>=lfc.thr)
  w2 <- which(X[, "padj"]<=fdr.thr & X[, "log2FoldChange"]<=(-1)*lfc.thr)
  r[w1] <- 1
  r[w2] <- -1
  return(r)
}

classifyPattern <- function(x) {
  if (x[1]==1 & x[2]==1) {
    return(1)
  }
  else {
    if (x[1]==1 & x[2]==-1) {
      return(2)
    }
    else {
      if (x[1]==0 & x[2]==1) {
        return(3)
      }
      else {
        if (x[1]==0 & x[2]==-1) {
          return(4)
        }
        else {
          if (x[1]==-1 & x[2]==1) {
            return(5)
          }
          else {
            if (x[1]==-1 & x[2]==-1) {
              return(6)
            }
            else {
              return(0)
            }
          }
        }
      }
    }
  }
}

createMemoryHeatmapBothStresses <- function(genes, indir, outdir, fdr.thr=0.05, lfc.thr=1, selectionfunction=selectUnion) {
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  tissues <- c("J", "A", "P", "X", "W")
  phases <- c("S", "R")
  ntissues <- length(tissues)
  combined.matrix.PS <- matrix(0, length(genes), ntissues)
  colnames(combined.matrix.PS) <- tissues
  rownames(combined.matrix.PS) <- genes
  combined.matrix.CS <- combined.matrix.PS
  for (i in 1:length(tissues)) {
    curr.PS <- matrix(0, length(genes), length(phases))
    curr.CS <- curr.PS
    for (j in 1:length(phases)) {
      load(sprintf("%scomparison_%s_%s_PS_vsEC.Rdata", indir, tissues[i], phases[j])) # comparison
      m <- match(genes, gsub(".v3.1", "", rownames(comparison), fixed=T)) 
      curr.PS[,j] <- classifyDE(comparison[m,], fdr.thr=fdr.thr, lfc.thr=lfc.thr)
      load(sprintf("%scomparison_%s_%s_CS_vsEC.Rdata", indir, tissues[i], phases[j])) # comparison
      m <- match(genes, gsub(".v3.1", "", rownames(comparison), fixed=T))
      curr.CS[,j] <- classifyDE(comparison[m,], fdr.thr=fdr.thr, lfc.thr=lfc.thr)
    }
    combined.matrix.PS[,i] <- apply(curr.PS, 1, classifyPattern)
    combined.matrix.CS[,i] <- apply(curr.CS, 1, classifyPattern)
  }
  
  combined.matrix <- selectionfunction(combined.matrix.PS, combined.matrix.CS)
  save(combined.matrix.PS, combined.matrix.CS, file=sprintf("%scombined.matrix.PSorCS.Rdata", outdir))
  
  # sort according to PS and plot
  crit1 <- apply((combined.matrix[,1:ntissues]>0)+0, 1, sum) 
  crit2 <- apply(combined.matrix[,1:ntissues], 1, function(x) return(max(table(x))))
  p <- paste(combined.matrix[,1], combined.matrix[,2], combined.matrix[,3], combined.matrix[,4], combined.matrix[,5], sep=",")
  sortindex <- c()
  row.breaks <- c(0)
  count <- 0
  s.down <- seq(ntissues, 2, -1) # focus on cross-tissue effect of PS # if second parameter=0, shows all CS genes as well
  for (i in s.down) {
    for (j in s.down){
      w <- which(crit1==i & crit2==j)
      curr.p <- p[w]
      unp <- sort(unique(curr.p))
      for (k in 1:length(unp)) {
        w.in <- which(curr.p==unp[k])
        sortindex <- c(sortindex, w[w.in])
        count <- count+length(w.in)
        if (length(w.in)>0) {
          row.breaks <- c(row.breaks, count)
        }
        
      }
      
    }
  }
  combined.matrix.sorted <- combined.matrix[sortindex,]
  
  library(gplots)
  patternCol <- c("white", "darkred", "lightblue", "red", "blue", "coral", "darkblue")
  colnames(combined.matrix.sorted) <- paste(c(rep("PS", ntissues), rep("CS", ntissues)), mapTissue(colnames(combined.matrix.sorted)), sep=" ") #changed
  pdf(sprintf("%s/summaryHeatmap.pdf", outdir), width=5, height=23) 
  heatmap.2(combined.matrix.sorted, Rowv=FALSE, Colv=FALSE, dendrogram="none", scale="none", col=patternCol, breaks=seq(-0.5,6.5,1), rowsep=row.breaks, colsep=c(0, ntissues, 2*ntissues), sepcolor="black", sepwidth=c(0.01,0.01), key=F, symkey=F, density.info="none", trace="none", margins=c(5,10), cexCol=1.1, cexRow=1.1, srtRow=0, srtCol=90, fg="black", bty="o", offsetRow=0.05, offsetCol=2.39, adjCol=c(NA,0.1))
  lix <- c(1,3,5,2,4,6)
  legend("top", legend=c("Up (R), up (S)", "Down (R), up (S)", "Up (R), unchanged (S)", "Down (R), unchanged (S)", "Up (R), down (S)", "Down (R), down (S)")[lix], fill=patternCol[c(2:length(patternCol))[lix] ], title="Gene expression memory", cex=1.1, ncol=1, xpd=T, x.intersp=1) 
  dev.off()
  
} 

doCoexpressionNetwork <- function (infile1, outdir, thr=0.3, my.cex=1.7, deepSplitParam=1, deepSplitParam2=T, minNrTissues=2, savefile="", omitGrey=F) {
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  # group average of log(TPM+1) 
  tissues <- c("J", "A", "P", "X", "W")
  load(infile1) # TPM_filtered
  groupid <- substr(rownames(TPM_filtered), start=3, stop=nchar(rownames(TPM_filtered)[1]))
  groups <- unique(groupid)
  # sort groups such that treatments and phases in the same order for each tissue
  matchStr.raw <- c("2_S", "1_S", "4_S", "2_R", "1_R", "4_R", "3_S", "3_R")
  matchStr <- c()
  for (i in 1:length(tissues)){
    matchStr <- c(matchStr, sprintf("%s_%s", matchStr.raw, tissues[i]))
  }
  groups <- matchStr
  
  TPM_filtered_logavg <- matrix(0, length(groups), ncol(TPM_filtered))
  for (i in 1:length(groups)){
    w <- which(groupid==groups[i])
    TPM_filtered_logavg[i,] <- apply(log2(TPM_filtered[w,]+1), 2, mean)
  }
  colnames(TPM_filtered_logavg) <- colnames(TPM_filtered)
  rownames(TPM_filtered_logavg) <- sprintf("m_%s", groups)
  
  # leave out AC conditions ("3")
  my.treatments <- extractPart(rownames(TPM_filtered_logavg), 2)
  wr <- which(my.treatments!="3")
  # Fig. S5
  data <- list()
  cum.thr <- seq(0, 1.9, 0.1)
  
  cum <- matrix(0, length(tissues), length(cum.thr)) # cumulative
  my.tissues <- extractPart(rownames(TPM_filtered_logavg[wr, ]), 4)
  for (i in 1:length(tissues)){
    mt <- which(my.tissues==tissues[i])
    coef.of.var <- apply(TPM_filtered_logavg[wr[mt], ], 2, function (x){ return(sd(x)/mean(x))})
    wc <- which(coef.of.var>=thr)
    data[[i]] <- TPM_filtered_logavg[wr[mt], wc]
    cum[i,] <- sapply(cum.thr, function(x){return(length(which(coef.of.var>=x)))})
    
  }
  pdf(sprintf("%splotCoefOfVar.pdf", outdir),20,7)
  par(mar=c(5.1,5.1,3.1,2.1))
  barplot(cum, names.arg=c("0", sprintf("%.1f", cum.thr[2:length(cum.thr)])), legend.text=mapTissue(tissues), beside=T, col=sapply(tissues, getTissuecolor), cex.axis=my.cex, cex.names=my.cex, cex.lab=my.cex, args.legend=list(cex=my.cex), xlab="Minimum coefficient of variation", ylab="Total number of genes")
  dev.off()
  
  # Compute module eigengene (ME) network according to Langfelder and Horvath, 2007
  library(WGCNA)
  library(flashClust)
  library(dynamicTreeCut)
  wgcna.tissue.MEs <- NULL
  wgcna.tissue.ModuleAssignment <- list()
  wgcna.tissue.ColorAssignment <- list()
  for (i in 1:length(tissues)){
    t.wgcna.sft <- pickSoftThreshold(data[[i]], networkType="signed hybrid", RsquaredCut=0.7)
    t.wgcna.sft.values <- -sign(t.wgcna.sft$fitIndices[,3]) * t.wgcna.sft$fitIndices[,2]
    t.wgcna.sft.powerEstimate <- which(t.wgcna.sft.values>=0.7)[1]
    t.wgcna.adjacency <- adjacency(data[[i]], type="signed hybrid", power=t.wgcna.sft.powerEstimate)
    t.wgcna.TOM <- TOMdist(t.wgcna.adjacency, TOMType = "signed") 
    t.wgcna.geneTree <- flashClust(as.dist(t.wgcna.TOM), method ="average")
    
    t.wgcna.cutreeDynamic <- cutreeDynamic(dendro=t.wgcna.geneTree, distM=t.wgcna.TOM, deepSplit=deepSplitParam, pamRespectsDendro = FALSE, minClusterSize=50)
    names(t.wgcna.cutreeDynamic) <- colnames(data[[i]])
    wgcna.tissue.ModuleAssignment[[i]] <- t.wgcna.cutreeDynamic
    t.wgcna.init.colors <- labels2colors(t.wgcna.cutreeDynamic)
    names(t.wgcna.init.colors) <- colnames(data[[i]])
    wgcna.tissue.ColorAssignment[[i]] <- t.wgcna.init.colors
    
    # MEs for initial modules
    t.wgcna.init.MEList <- moduleEigengenes(data[[i]], colors=t.wgcna.init.colors, excludeGrey=omitGrey)
    
    t.wgcna.init.MEs    <- t.wgcna.init.MEList$eigengenes
    colnames(t.wgcna.init.MEs) <- sub("^ME", "", colnames(t.wgcna.init.MEs))
    
    my.names <- c(colnames(wgcna.tissue.MEs), sprintf("%s.%s", mapTissue(tissues[i]),colnames(t.wgcna.init.MEs)))
    wgcna.tissue.MEs <- cbind(wgcna.tissue.MEs, as.matrix(t.wgcna.init.MEs))
    colnames(wgcna.tissue.MEs) <- my.names
  }
  my.cond <- extractPart(rownames(data[[1]]), part=2, spl="_")
  my.phase <- extractPart(rownames(data[[1]]), part=3, spl="_")
  rownames(wgcna.tissue.MEs) <- sprintf("%s (%s)", sapply(my.cond, mapCond), my.phase)
  
  wgcna.tissue.MEDiss <- 1 - cor(wgcna.tissue.MEs) # dist measure suggested by Langfelder and Horvath, 2007
  wgcna.tissue.METree <- flashClust(as.dist(wgcna.tissue.MEDiss), method ="average")
  comm <- cutreeDynamic(dendro=wgcna.tissue.METree, distM=wgcna.tissue.MEDiss, deepSplit=deepSplitParam2, pamRespectsDendro = FALSE, minClusterSize=2)
  ## comm is 0 for not assigned members, 1 for biggest cluster etc.
  nrComm <- max(comm)
  
  # plot wgcna.tissue.MEs in blocks according to communities (to check community profiles)
  most.extreme <- max(abs(min(min(wgcna.tissue.MEs))), abs(max(max(wgcna.tissue.MEs))))
  print(most.extreme)
  nrSteps <- 75
  my.breaks <- seq((-1)*most.extreme, most.extreme, most.extreme/nrSteps) # in all plots the same color scheme
  library(gplots)
  selectedComm <- c()
  for (i in 1:nrComm){
    w <- which(comm==i)
    curr.tissues <- extractPart(colnames(wgcna.tissue.MEs)[w], part=1, spl=".")
    if (length(unique(curr.tissues))>=minNrTissues) {
      pdf(sprintf("%sHeatmap_comm_%d.pdf", outdir, i))
      heatmap.2(wgcna.tissue.MEs[,w], col=bluered(length(my.breaks)-1), Rowv=F, Colv=F, density.info="none", trace="none", dendrogram="none", symm=F, symkey=F, symbreaks=T, breaks=my.breaks, scale="none") 
      dev.off()
      selectedComm <- c(selectedComm, i)
      
    }
  }
  # save all information for network plot
  if (savefile != "") {
    wgcna.tissue.MECor <- cor(wgcna.tissue.MEs)
    wgcna.tissue.MECorPval <- corPvalueFisher(wgcna.tissue.MECor, dim(wgcna.tissue.MEs)[1])
    save(selectedComm, comm, wgcna.tissue.ModuleAssignment, wgcna.tissue.ColorAssignment, wgcna.tissue.MEs, wgcna.tissue.MECor, wgcna.tissue.MECorPval, file=savefile)
  }
}

doCoexpressionFigure <- function (infile, outdir, my.cex=1.7, minNodeSize=10, maxNodeSize=20) {
  load(infile)
  
  # heatmaps in two rows
  memory.heatmaps <- c(2,3,8,9,10,12,14) 
  other.heatmaps <- setdiff(selectedComm, memory.heatmaps)
  
  # adjust community and module numbering (module numbering per tissue) 
  # get module size
  commIdNew <- 1:(length(other.heatmaps)+length(memory.heatmaps))
  names(commIdNew) <- c(other.heatmaps, memory.heatmaps)
  w <- which(comm %in% c(other.heatmaps,memory.heatmaps))

  moduleIdOld <- colnames(wgcna.tissue.MEs)[w]
  my.tissues <- extractPart(colnames(wgcna.tissue.MEs), part=1, spl=".")[w] 
  moduleIdNew <- moduleIdOld # initialization, content changed below
  names(moduleIdNew) <- moduleIdOld
  annot <- c("LE1", "LE2", "PHL", "XYL", "ROO")
  
  moduleSizeNew <- rep(0, length(moduleIdNew))
  allocationCommand <- sprintf("list(%s)", paste(rep("rep(0,5000)", length(moduleIdNew)), collapse=","))
  moduleGenesNew <- eval(parse(text=allocationCommand)) 
  moduleColorNew <- rep("black", length(moduleIdNew))
  labelColorNew <- rep("black", length(moduleIdNew))
  for (i in 1:length(moduleIdNew)){
    curr.tissue <- extractPart(names(moduleIdNew)[i], part=1, spl=".")
    curr.tissue.ix <- which(annot==curr.tissue)
    curr.col <-  extractPart(names(moduleIdNew)[i], part=2, spl=".")
    wc <- which(wgcna.tissue.ColorAssignment[[curr.tissue.ix]]==curr.col)
    moduleGenesNew[[i]] <- names(wgcna.tissue.ColorAssignment[[curr.tissue.ix]])[wc]
    moduleSizeNew[i] <- length(moduleGenesNew[[i]])
    moduleColorNew[i] <- mapTissueColor(curr.tissue)
    if (moduleColorNew[i] %in% c("green4", "mediumblue", "deeppink2")) {
      labelColorNew[i] <- "white"
    }
  }
  # set new moduleIds, sort according to size in each tissue
  for (i in 1:length(annot)){
    wt <- which(my.tissues==annot[i])
    sizesort <- sort(moduleSizeNew[wt], decreasing=T, index.return=T)
    moduleIdNew[wt[sizesort$ix] ] <- sprintf("%s.%d", annot[i], 1:length(wt))
  }
  
  # for heatmap map old module names to new module names
  mapModuleName <- function(x){
    wm <- which(names(moduleIdNew)==x)
    return(moduleIdNew[wm[1] ])
  }
  mapCommIndex <- function(x){
    wm <- which(names(commIdNew)==x)
    return(commIdNew[wm[1] ])
  }
  
  library(gplots) # for bluered
  library(pheatmap)
  library(grid)
  library(gtable)
  source("pheatmapEdited.R") # changed fonts: fontface="plain", fontsize = 1.0 * fontsize
  most.extreme <- max(abs(min(min(wgcna.tissue.MEs))), abs(max(max(wgcna.tissue.MEs))))
  print(most.extreme)
  nrSteps <- 75
  my.breaks <- seq((-1)*most.extreme, most.extreme, most.extreme/nrSteps) # in all plots the same color scheme
  
  plot_list <- list()
  count <- 1
  for (i in memory.heatmaps) {
    w <- which(comm==i)
    my.mat <- wgcna.tissue.MEs[,w]
    colnames(my.mat) <- sapply(colnames(my.mat), mapModuleName)
    if (i==14){
      row.TF=T
    }
    else {
      row.TF=F
    }
    ph <- pheatmap(my.mat, color=bluered(length(my.breaks)-1), breaks=my.breaks, border_color="black", cellwidth=10, cellheight=10, cluster_rows = F, cluster_cols = F, legend=F, show_rownames=row.TF, gaps_row=c(3), gaps_col=1:(ncol(my.mat)-1), main=sprintf("C%d", mapCommIndex(i)), cex.main=1, rot_col = 90, hjust_col = 1, vjust_col = 0.5)
    plot_list[[count]] <- ph[[4]]
    count <- count + 1
    a <- gtable(unit(0.45, c("inch")), unit(rep(0.1,5), c("inch")))
    gr <- rectGrob(x=0, y=0, width = unit(0.2, "npc"), height = unit(0.2, "npc"), gp=gpar(col="white"))
    a <- gtable_add_grob(a, gr, 1, 1)
    plot_list[[count]] <- a
    count <- count + 1
    
  }
  for (i in other.heatmaps) {
    w <- which(comm==i)
    my.mat <- wgcna.tissue.MEs[,w]
    colnames(my.mat) <- sapply(colnames(my.mat), mapModuleName)
    if (i==13){
      row.TF=T
    }
    else {
      row.TF=F
    }
    ph <- pheatmap(my.mat, color=bluered(length(my.breaks)-1), breaks=my.breaks, border_color="black", cellwidth=10, cellheight=10, cluster_rows = F, cluster_cols = F, legend=F, show_rownames=row.TF, gaps_row=c(3), gaps_col=1:(ncol(my.mat)-1), main=sprintf("C%d", mapCommIndex(i)), cex.main=1, rot_col = 90, hjust_col = 1, vjust_col = 0.5) 
    plot_list[[count]] <- ph[[4]]
    count <- count + 1
    a <- gtable(unit(0.45, c("inch")), unit(rep(0.1,5), c("inch")))
    gr <- rectGrob(x=0, y=0, width = unit(0.2, "npc"), height = unit(0.2, "npc"), gp=gpar(col="white"))
    a <- gtable_add_grob(a, gr, 1, 1)
    plot_list[[count]] <- a
    count <- count + 1
  }
  library(ggplot2)
  library(gridExtra)
  
  g1 <- do.call(cbind, c(plot_list[1:13], size="first"))
  pdf(sprintf("%sheatmaps_memory.pdf", outdir), width=12, height=7)
  do.call(grid.draw, list(g1))
  dev.off()
 
  g1 <- do.call(cbind, c(plot_list[15:23], size="first"))
  pdf(sprintf("%sheatmaps_other.pdf", outdir), width=8.6, height=7)
  do.call(grid.draw, list(g1))
  dev.off()
  
  # separate figure for legend
  pdf(sprintf("%sheatmap_legend.pdf", outdir))
  pheatmap(my.mat, color=bluered(length(my.breaks)-1), breaks=my.breaks, border_color="black", cellwidth=10, cellheight=10, cluster_rows = F, cluster_cols = F, legend=T, legend_breaks=my.breaks, legend_labels=c(sprintf("%.1f", (-1)*most.extreme), rep("", length(my.breaks)-2),sprintf("%.1f", most.extreme)), main=expression(paste("Gene expression [", 'log'[2], "(TPM+1)]", sep="")), cex.main=my.cex, font.main=1)
  dev.off()
  
  # separate figure for graph
  library(igraph)
  
  # map adjacency matrix to new order, select significant edges
  m <- match(names(moduleIdNew), colnames(wgcna.tissue.MEs))
  wps <- which(wgcna.tissue.MECorPval[m,m]>1 | wgcna.tissue.MECor[m,m]<=0.7, arr.ind=T) # not (positive and significant)
  adj.matrix <- wgcna.tissue.MECor[m,m]
  adj.matrix[wps] <- 0
  diag(adj.matrix) <- 0
  ga <- graph_from_adjacency_matrix(adj.matrix, weighted=T, mode="undirected")
  ## vertex groups from comm 
  my.comm <- comm[m]
  my.comm.list <- list()
  for (i in 1:length(names(commIdNew))){
    wcomm <- which(my.comm==as.numeric(names(commIdNew)[i]))
    my.comm.list[[i]] <- wcomm
  }
  
  moduleSizeNewScaled <- minNodeSize+(moduleSizeNew-min(moduleSizeNew))/(max(moduleSizeNew)-min(moduleSizeNew))*(maxNodeSize-minNodeSize)
  
  # Uncomment the following lines to interactively adjust the network plot and save its coordinates
  #id <- tkplot(ga, canvas.width = 600, canvas.height = 450, mark.groups=my.comm.list, mark.expand=10, mark.shape=0.5, edges.curved=F, vertex.size=moduleSizeNewScaled, vertex.color=moduleColorNew, vertex.label="", vertex.label.cex=my.cex, vertex.label.color=labelColorNew, vertex.label.family="sans", edge.width=3, edge.color="black") # tkplot does not accept vector for edge.width
  #my.coords <- tk_coords(id)
  #save(my.coords, file=sprintf("%scoords.Rdata", outdir))
  load(sprintf("%scoords.Rdata", outdir))
  # rescale coords to [-1,1] 
  ## first rescale to [0,2] and then subtract 1
  my.max.coord <- max(max(my.coords[,1]), max(my.coords[,2]))
  my.coords.rescaled.raw <- my.coords*2/my.max.coord
  my.coords.rescaled <- my.coords.rescaled.raw-1
  transformEdgeWeight <- function(x){
    return(1+((x-0.7)/(1-0.7))*(4-1)) # range [0.7,1] to range [1,4]
  }
  
  pdf(sprintf("%sgraph.pdf", outdir), width=14, height=12)
  plot(ga, mark.groups=my.comm.list, mark.expand=4, mark.col = gray(c(1:length(my.comm.list))/(length(my.comm.list)*2), alpha = 0.3), mark.border = gray(c(1:length(my.comm.list))/(length(my.comm.list)*2), alpha = 1), vertex.size=moduleSizeNewScaled, vertex.color=moduleColorNew, vertex.label=moduleIdNew, vertex.label.cex=my.cex, vertex.label.color=labelColorNew, vertex.label.family="sans", layout=my.coords.rescaled, rescale=F, asp=0.8, margin=0.1, edge.color=rep("black",length(E(ga)$weight)), edge.width=transformEdgeWeight(E(ga)$weight)) 
  ## add community labels
  e <- 0.1
  markernodes <- sprintf("ROO.%d", 1:12)
  markernodes.shiftx <- c(0,-e,0,e,-e,-(1.2*e),0,e,-e,e,0,-e)
  markernodes.shifty <- c(-(1.2*e),0,-e,0,0,0,-e,0,0,0,-e,0)
  mn <- match(markernodes, moduleIdNew)
  my.comm.labels <- sprintf("C%d",sapply(comm[ m[mn] ], mapCommIndex))
  text(my.coords.rescaled[mn,1]+markernodes.shiftx, my.coords.rescaled[mn,2]+markernodes.shifty, my.comm.labels, cex=my.cex)
  dev.off()

  # separate legend for node size
  sf <- 3
  pdf(sprintf("%sgraph_sizeLegend.pdf", outdir))
  plot(graph_from_adjacency_matrix(matrix(c(0,1),2,2)), vertex.size=c(sf*min(moduleSizeNewScaled), sf*max(moduleSizeNewScaled)), edge.color="white", vertex.label=sprintf("%.0f", c(min(moduleSizeNew), max(moduleSizeNew))), vertex.label.dist=0, vertex.color="white", vertex.label.cex=my.cex, main="", cex.main=my.cex, font.main=1, layout=matrix(c(0.1, 0.1,-0.01,0.01), 2,2,byrow=F), vertex.label.family="sans", vertex.label.color="black", vertex.label.degree=0)
  dev.off()
  
  # separate legend for edge width 
  pdf(sprintf("%sgraph_edgeLegend.pdf", outdir))
  plot(c(1,1), c(0.8,1.0), type="l", lwd=4, ylim=c(0,2), xlim=c(0,2), ljoin=1)
  lines(c(1,1), c(0.4,0.6), type="l", lwd=1, ljoin=1)
  text(c(1.1,1.1),c(0.9,0.5), c("1.0", "0.7"))
  dev.off()
  
  # separate legend for node color
  pdf(sprintf("%sgraph_tissueLegend.pdf", outdir))
  plot(graph_from_adjacency_matrix(matrix(0,5,5)), vertex.color=sapply(annot, mapTissueColor), vertex.label=annot, edge.color="white", vertex.label.dist=1, vertex.label.cex=my.cex, main="Tissue", cex.main=my.cex, font.main=1)
  dev.off()
}

doCorrelationAnalysis <- function (infile, outdir, isfile, my.thr=0.8) {
  load(infile) # treedata, TPM_filtered2
  e <- extractPart(colnames(treedata), spl=".", part=2)
  # correlation between the same gene across different tissue pairs (full data)
  ## compute a genes x (tissue pairs) matrix with correlations
  ix <- seq(1,ncol(treedata), ncol(treedata)/5)
  addon <- c(1:(ncol(treedata)/5))-1
  singleGeneTissuePairs <- function(a){ 
    cm <- cor(log2(treedata[, ix+a]+1), use="pairwise.complete.obs") # use log level as in all other analyses
    return(cm[upper.tri(cm)]) # goes first through columns
  }
  res <- sapply(addon, singleGeneTissuePairs)
  val <- t(res)
  rownames(val) <- sprintf("Potri.%s", e[1:nrow(val)])
  colnames(val) <- c("AJ", "WJ", "WA", "PJ", "PA", "PW", "XJ", "XA", "XW", "XP")
  
  a2 <- apply(val, 1, median)
  w2 <- which(!is.na(a2))
  sres2 <- sort(a2[w2], decreasing=T, index.return=T)
  
  glob <- read.delim(isfile, header=F, stringsAsFactors=F) 
  mrank <- match(gsub(".v3.1", "", glob[,1]), names(a2)[w2[sres2$ix] ])
  
  w08 <- which(a2[w2]>my.thr)
  r <- length(w08)/length(w2)
  ws <- which(mrank<=length(w08))
  
  my.bg <- rep("white", length(w08))
  my.bg[mrank[ws] ] <- "black" 

  pdf(sprintf("%smedian_top.pdf", outdir))
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(1:length(w08), sres2$x[1:length(w08)], type="p", col="black", pch=21, bg=my.bg, xlab=sprintf("Top %.1f%% of genes (median self-correlation > %.1f)", r*100, my.thr), ylab="Median self-correlation across tissue pairs", cex=1.7, cex.lab=1.7, cex.axis=1.7, xlim=c(0, length(w08)+5)) 
  txtpos <- mrank[ws] 
  add.y <- rep(0, length(txtpos)) 
  add.x <- add.y
  add.y[c(6,2,10,4,5,3,8,14,17,16,7,1,11)] <- c(0.0099+0.002, 0.0069+0.002, 0.004, 0.0039+0.002,0.0006+0.002,-0.002+0.002,-0.001+0.002,-0.0042+0.002,-0.0079+0.002,-0.011+0.002,0.004, 0.004,-0.002)
  wtxt <- which(sres2$x[txtpos]<0.833)
  add.x[wtxt] <- -31
  add.x[c(8,14,17,16)] <- 8
  add.y[wtxt] <- add.y[wtxt]-0.002
  text(txtpos+add.x, sres2$x[txtpos]+add.y, labels=gsub(".v3.1", "", glob[,1])[ws], cex=1.7, pos=4, offset=1)
  dev.off()

  write.table(gsub(".v3.1", "", glob[,1])[ws], file=sprintf("%scorrmemgenes.txt", outdir), sep="\t", col.names=F, row.names=F)
  
  ## pairwise tissue heatmap  # black color scheme
  getReadout <- function(x){ # number of genes with high self-correlation across tissues
    return(length(which(x>my.thr))) 
  }
  valcondensed <- apply(val,2,getReadout)
  tissues <- c("J", "A", "W", "P", "X")
  tissuematrix1 <- matrix(0, length(tissues), length(tissues))
  tissuematrix1[upper.tri(tissuematrix1)] <- valcondensed # need to make complete before reordering; better leave like that, then extremes in corners
  rownames(tissuematrix1) <- mapTissue(tissues) 
  colnames(tissuematrix1) <- rownames(tissuematrix1)
  library(gplots)
  pdf(sprintf("%sheatmap_crossTissueSelfCorrelatedGenes.pdf", outdir))
  heatmap.2(tissuematrix1[1:4,2:5], Rowv=NA, Colv=NA, dendrogram="none", scale="none", col=colorpanel(125,"white", "black"), trace="none", density.info="none", srtCol=0, adjCol=c(NA,-25)) 
  dev.off()
  
}

makeEdges <- function(x, yvec, val) {
  return(data.frame(pcor=sign(val), node1=rep(x, length(yvec)), node2=yvec)) # GeneNet format
}

# GRN only for genes of interest (goi)
findNeighbors <- function(goi, tpmfile, tfs, outdir, maxNrZeros=0, method="glmnet", qval.thr=0.05, myalpha=1, top=0, tfInput=T, tfOutput=T, myk="sqrt", myntrees=1000) {
  tissues <- c("J", "A", "P", "X", "W")
 
  ## do GRN for each tissue for each goi
  load(tpmfile)
  
  # should output or input of the GRN learning function be restricted to candidates contained in tfs?
  if (tfOutput) {
    m <- match(tfs, colnames(TPM_filtered))
    if (tfInput){
      data <- log2(TPM_filtered[,m]+1)
    }
    else {
      m2 <- setdiff(match(goi, colnames(TPM_filtered)),m)
      data <- log2(TPM_filtered[,c(m2,m)]+1)
    }
  }
  else {
    data <- log2(TPM_filtered+1)
  }
  
  countZero <- function(x) {
    return(length(which(x==0 | is.na(x))))
  }
  a <- apply(data, 2, countZero)
  w <- which(a<=maxNrZeros)
  input.data.raw <- data[,w]
  sample.tissue <- extractPart(rownames(data), part=4, spl="_")
  
  mg <- match(goi, colnames(input.data.raw))
  mg <- mg[which(!is.na(mg))]
  
  for (i in 1:length(tissues)) {
    wt <- which(sample.tissue==tissues[i])
    input.data.tissue <- input.data.raw[wt,]
    
    if (method=="glmnet"){ 
      library(glmnet)
      df <- data.frame(pcor=c(), node1=c(), node2=c(), qval=c())
      for (j in 1:length(mg)){
        input2 <- input.data.tissue[,-mg[j]]
        glmodel <- cv.glmnet(input2, input.data.tissue[,mg[j]], alpha=myalpha, nfolds=nrow(input2), grouped=F)
        glmodel.coef <- coef(glmodel)
        wp <- which(glmodel.coef[,1]!=0)
        mp <- match(colnames(input2)[wp], colnames(input.data.tissue))
        df <- rbind(df, makeEdges(mg[j], mp, glmodel.coef[wp,1]))
      }
    }
    else {
      if (method=="GENIE3"){
        library(GENIE3)
        df <- data.frame(pcor=c(), node1=c(), node2=c(), qval=c())
        for (j in 1:length(mg)){
          set.seed(123)
          myreg <- setdiff(1:ncol(input.data.tissue), c(mg[j]))
          genmodel <- GENIE3(t(input.data.tissue), regulators=colnames(input.data.tissue)[myreg], targets=colnames(input.data.tissue)[mg[j] ], K=myk, nTrees=myntrees)
          # resort genmodel according to myreg
          mm <- match(colnames(input.data.tissue)[myreg], rownames(genmodel))
          genmodel <- genmodel[mm, , drop=F]
          if (top>0) {
            sr <- sort(genmodel[,1], decreasing=T, index.return=T)
            wp <- sr$ix[1:min(top, nrow(genmodel))]
          }
          else {
            wp <- which(genmodel[,1]!=0)
          }
          mp <- myreg[wp]
          df <- rbind(df, makeEdges(mg[j], mp, genmodel[wp,1]))
        }
      }
    }
  
    node1names <- colnames(input.data.tissue)[df[, "node1"] ]
    node2names <- colnames(input.data.tissue)[df[, "node2"] ]
    save(node1names, node2names, file=sprintf("%s%s_%s.Rdata", outdir, method, tissues[i]))
  }
}

sortNodes <- function(X) {
  wx <- which(X$node1>X$node2)
  tmp <- X$node1[wx]
  X$node1[wx] <- X$node2[wx]
  X$node2[wx] <- tmp
  return(X)
}

intersectEdges <- function(X,Y){
  # direction of edge not taken into account
  X <- sortNodes(X)
  Y <- sortNodes(Y)
  px <- paste(X$node1, X$node2, sep="_")
  py <- paste(Y$node1, Y$node2, sep="_")
  is <- intersect(px, py)
  wx <- which(px %in% is)
  wy <- which(py %in% is)
  if (length(wx)>0) {
    w <- which(sign(X$pcor[wx])==sign(Y$pcor[wy]))
    
    return(X[wx[w],])
  }
  else {
    return(X[wx,])
  }
  
  
}

combineEdges <- function(X,Y){
  X <- sortNodes(X)
  Y <- sortNodes(Y)
  px <- paste(X$node1, X$node2, sep="_")
  py <- paste(Y$node1, Y$node2, sep="_")
  is <- intersect(px, py)
  wy <- which(py %in% is)
  if (length(wy)>0) {
    return(rbind(X, Y[-wy,]))
  }
  else {
    return(rbind(X,Y))
  }
  
}

findNeighborsEach <- function(goi, tpmfile, tfs, outdir, maxNrZeros=0, method="glmnet", qval.thr=0.05, myalpha=1, top=0, tfInput=T, tfOutput=T) {
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  for (i in 1:length(goi)){
    if (!dir.exists(sprintf("%sgoi_%d/",outdir,i))){
      dir.create(sprintf("%sgoi_%d/",outdir,i))
    }
    findNeighbors(goi[i], tpmfile, tfs, sprintf("%sgoi_%d/",outdir,i), method=method, maxNrZeros=maxNrZeros, top=top, tfInput=tfInput, tfOutput=tfOutput)
  }
  is.edges <- data.frame() # intersection across tissues
  is.nodes <- c()
  tissue.edges.list <- list() # edges for each tissue
  tissues <- c("J", "A", "W", "P", "X") 
  for (i in 1:length(tissues)){ 
    union.edges <- data.frame() # combining across gois for each tissue
    
    for (j in 1:length(goi)){
      load(sprintf("%sgoi_%d/%s_%s.Rdata", outdir, j, method, tissues[i]))
      df <- data.frame(pcor=rep(1, length(node1names)), node1=node1names, node2=node2names, stringsAsFactors=F)
      if (j==1){
        union.edges <- df
      }
      else {
        is.edges <- intersectEdges(is.edges, df)
        union.edges <- combineEdges(union.edges, df) # overlap possible, mutual best regulators for two goi
      }
      
    }
    if (i==1){
      is.edges <- union.edges
      is.nodes <- union(union.edges$node1, union.edges$node2)
    }
    else {
      is.edges <- intersectEdges(is.edges, union.edges)
      is.nodes <- intersect(is.nodes, union(union.edges$node1, union.edges$node2))
    }
    tissue.edges.list[[i]] <- union.edges
    
  }
  node.list <- list(c(),c(),c(),c(),c())
  all.nodes <- c() # without goi
  all.edges <- data.frame()
  for (i in 1:length(tissues)){
    a <- tissue.edges.list[[i]]
    node.list[[i]] <- setdiff(union(a$node1, a$node2), goi)
    all.nodes <- union(all.nodes, node.list[[i]])
    all.edges <- combineEdges(all.edges, a)
  }
  node.pairwise <- matrix(0, length(tissues), length(tissues))
  for (i in 1:length(tissues)){
    for (j in 1:length(tissues)){
      node.pairwise[i,j] <- length(intersect(node.list[[i]],node.list[[j]]))
    }
  }
  
  node.occurrence <- matrix(0, length(all.nodes), length(tissues))
  rownames(node.occurrence) <- all.nodes
  for (i in 1:length(tissues)){
    m <- match(node.list[[i]], all.nodes)
    node.occurrence[m,i] <- 1
  }
  
  # edge intersections/occurrence
  edge.occurrence <- matrix(0, nrow(all.edges), length(tissues))
  rownames(edge.occurrence) <- paste(all.edges$node1, all.edges$node2, sep="_")
  for (i in 1:length(tissues)){
    a <- tissue.edges.list[[i]]
    m <- match(paste(a$node1, a$node2, sep="_"), paste(all.edges$node1, all.edges$node2, sep="_"))
    edge.occurrence[m,i] <- 1
  }
  
  save(tissues, goi, tissue.edges.list, node.list, node.pairwise, node.occurrence, edge.occurrence, file=sprintf("%sgraphdata.Rdata", outdir))
  #apn <- apply(node.occurrence, 1, sum)
  #wa <- which(apn>2)
  #node.occurrence[wa,] # 4 patterns, see getIntersectColor
}

getIntersectColor <- function(x, tcolor, curr.i, edge=F){
  s <- sum(x)
  if (s==1){
    if (edge){
      return("black")
    }
    else {
      return("white")
    }
  }
  else {
    if (s==2){
      w <- which(x==1)
      wr <- setdiff(w, curr.i)
      return(tcolor[wr])
    }
    else {

      if (length(which(x==c(1,0,1,1,1)))==length(x)){
        return("darkviolet")
      }
      else {
        if (length(which(x==c(1,0,1,1,0)))==length(x)){
          return("deepskyblue1")
        }
        else {
          if (length(which(x==c(1,0,0,1,1)))==length(x)){
            return("orange")
          }
          else {
            if (length(which(x==c(0,1,0,1,1)))==length(x)){
              return("chocolate") 
            }
            else {
              return("bisque") 
            }
          }
        }
      }
    }
  }
}


colorGRN <- function (infile, outdir) {
  library(igraph)
  load(infile)
  tissue.color <- sapply(tissues, getTissuecolor) # tissues from infile
  for (i in 1:length(tissues)){
    curr.nodes <- c(goi, node.list[[i]])
    a <- tissue.edges.list[[i]]
    adj.matrix <- matrix(0, length(curr.nodes), length(curr.nodes))
    m1 <- match(a$node1, curr.nodes)
    m2 <- match(a$node2, curr.nodes)
    adj.matrix[cbind(m1,m2)] <- 1
    adj.matrix[cbind(m2,m1)] <- 1
    ga <- graph_from_adjacency_matrix(adj.matrix, mode="undirected")
    
    # get node colors; goi grey
    # V(ga) # vertex sequence is in the order from curr.nodes
    mo <- match(curr.nodes, rownames(node.occurrence))
    vc <- rep("white", length(curr.nodes))
    vc[1:length(goi)] <- "gray"
    ix <- (length(goi)+1):length(vc)
    vc[ix] <- apply(node.occurrence[mo[ix], ], 1, getIntersectColor, tcolor=tissue.color, curr.i=i)
    
    # get edge colors
    ## E(ga)$weight
    edgeorder <- ends(ga, E(ga), names=F) 
    edgeorderNamesDf <- data.frame(pcor=rep(1, nrow(edgeorder)), node1=curr.nodes[edgeorder[,1] ], node2=curr.nodes[edgeorder[,2] ], stringsAsFactors=F)
    so <- sortNodes(edgeorderNamesDf)
    p <- paste(so$node1, so$node2, sep="_")
    moe <- match(p, rownames(edge.occurrence))
    ec <- apply(edge.occurrence[moe,], 1, getIntersectColor, tcolor=tissue.color, curr.i=i, edge=T)
    
    # draw graph
    pdf(sprintf("%scolorgraph_%s.pdf", outdir, tissues[i]))
    plot(ga, vertex.size=10, vertex.color=vc, vertex.label="", edge.color=ec, edge.width=3, layout=layout_with_graphopt) 
    rect(-1.14, 0.86, -0.86, 1.14, border=tissue.color[i], col=NA, lwd=3)
    text(-1, 1, label=mapTissue(tissues[i]), cex=1.7)
    dev.off()
  }
  # plot legend separately
  pdf(sprintf("%scolorgraph_legend.pdf", outdir))
  plot(c(-1,1),c(-1,1))
  tissues.sorted <- c("LE1", "LE2", "PHL", "XYL", "ROO")
  legendMultiple <- function(x){
    mt <- mapTissue(tissues[which(x==1)])
    
    ma <- match(mt, tissues.sorted)
    ma.sorted <- sort(ma)
    mt.sorted <- tissues.sorted[ma.sorted]
    return(paste(mt.sorted, collapse=", "))
  }
  legend("center", legend=c("Query gene (self-correlated)", sprintf("Additional occurrence in %s", tissues.sorted), sprintf("Occurrence in %s", legendMultiple(c(1,0,1,1,1))), sprintf("Occurrence in %s", legendMultiple(c(1,0,1,1,0))), sprintf("Occurrence in %s", legendMultiple(c(1,0,0,1,1))), sprintf("Occurrence in %s", legendMultiple(c(0,1,0,1,1)))), fill=c("gray", tissue.color[c(1,2,4,5,3)], "darkviolet", "deepskyblue1", "orange", "chocolate"), title="Occurrence across tissues", cex=1.7) # Node/edge
  dev.off()
}

mapTfcolor <- function(x) {
  if (is.na(x)){
    return("white")
  }
  else {
    library(hash)
    famhash <- hash(keys=c("HD-ZIP", "bZIP", "TCP", "AP2", "ERF"), values=c("bisque3", "bisque", "cadetblue3", "darkolivegreen4", "darkkhaki"))
    return(famhash[[sprintf('%s', x)]])
  }
}

myellipse <- function(coords, v=NULL, params) {
  library(plotrix)
  # from https://stackoverflow.com/questions/48457824/how-to-plot-a-ellipse-node-with-igraph?noredirect=1&lq=1
  
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/30 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  draw.ellipse(x=coords[,1], y=coords[,2], a = vertex.size/2.5, b=vertex.size/5, col=vertex.color)
}

myellipse2 <- function(coords, v=NULL, params) {
  library(plotrix)
  # from https://stackoverflow.com/questions/48457824/how-to-plot-a-ellipse-node-with-igraph?noredirect=1&lq=1
  
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/30 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  draw.ellipse(x=coords[,1], y=coords[,2], a = vertex.size, b=vertex.size/2, col=vertex.color)
}

makeTissueGraphTfcolorIgraph <- function(el, prefix, widths=rep(1, nrow(el)), tf.file){
  ptf <- read.delim(tf.file, check.names=F, stringsAsFactors=F)
  curr.nodes <- union(el$node1, el$node2)
  m <- match(curr.nodes, ptf[,2]) # NA for non-TF
  w <- which(!is.na(m))
  library(igraph)
  adj.matrix <- matrix(0, length(curr.nodes), length(curr.nodes))
  m1 <- match(el$node1, curr.nodes)
  m2 <- match(el$node2, curr.nodes)
  adj.matrix[cbind(m1,m2)] <- widths
  adj.matrix[cbind(m2,m1)] <- widths
  ga <- graph_from_adjacency_matrix(adj.matrix, weighted=T, mode="undirected")
  
  add_shape("ellipse", clip=shapes("circle")$clip, plot=myellipse)
  #browser()
  #id <- tkplot(ga, canvas.width = 900, canvas.height = 310, vertex.shape="ellipse", vertex.color=sapply(ptf[m, 3], mapTfcolor), edges.curved=F, vertex.size=25, vertex.label=curr.nodes, layout=layout_nicely, vertex.label.family="sans", vertex.label.cex=1.7, vertex.label.color="black", edge.width=E(ga)$weight, edge.color=rep("black",length(E(ga)$weight)))
  ## can be changed interactively
  #my.coords <- tk_coords(id)
  #save(my.coords, file=sprintf("%scoords.Rdata", prefix))
  load(sprintf("%scoords.Rdata", prefix)) # my.coords
  pdf(sprintf("%s_igraph.pdf", prefix), width=21, height=7)
  plot(ga, vertex.shape="ellipse", vertex.color=sapply(ptf[m, 3], mapTfcolor), edges.curved=F, vertex.size=15, vertex.label=curr.nodes, layout=my.coords, rescale=T, asp=0.3, vertex.label.family="sans", vertex.label.cex=1.7, vertex.label.color="black", edge.width=1.3*E(ga)$weight, edge.color=rep("black",length(E(ga)$weight)))
  dev.off()
  
  # legend
  pdf(sprintf("%s_igraphlegend.pdf", prefix), width=21, height=7)
  plot(c(-1,1),c(-1,1))
  myf <- c("HD-ZIP", "bZIP", "TCP", "AP2", "ERF")
  legend("center", legend=myf, fill=sapply(myf, mapTfcolor), title="Transcription factor family", ncol=5, cex=1.7)
  dev.off()
}

makeIntersectingEdgeGraph <- function (em, prefix, tf.file) {
  ap <- apply(em, 1, sum)
  # extract edges from names
  df <- data.frame(pcor=rep(1, nrow(em)), node1=extractPart(rownames(em), part=1, spl="_"), node2=extractPart(rownames(em), part=2, spl="_"), stringsAsFactors=F)
  makeTissueGraphTfcolorIgraph(df, prefix, widths=ap, tf.file) 
}

makeTissueGraphCoreIgraph <- function (el, prefix, tname, decision.file, my.genes, my.color) {
  curr.nodes <- union(el$node1, el$node2)
  adj.matrix <- matrix(0, length(curr.nodes), length(curr.nodes))
  m1 <- match(el$node1, curr.nodes)
  m2 <- match(el$node2, curr.nodes)
  adj.matrix[cbind(m1,m2)] <- 1
  adj.matrix[cbind(m2,m1)] <- 1
  
  # iteratively reduce graph to core (remove single-edge nodes)
  adj.matrix.deleted <- adj.matrix
  alreadyDeleted <- c()
  curr.ix <- setdiff(c(1:length(curr.nodes)), alreadyDeleted)
  s <- apply(adj.matrix.deleted[curr.ix, curr.ix], 1, sum)
  w <- which(s==1)
  while (length(w)>0){
    adj.matrix.deleted[curr.ix[w],] <- 0
    adj.matrix.deleted[, curr.ix[w]] <- 0
    alreadyDeleted <- c(alreadyDeleted, curr.ix[w])
    curr.ix <- setdiff(c(1:length(curr.nodes)), alreadyDeleted)
    s <- apply(adj.matrix.deleted[curr.ix, curr.ix], 1, sum)
    w <- which(s<=1)
  }
  curr.nodes <- curr.nodes[curr.ix]
  adj.matrix.reduced <- adj.matrix[curr.ix, curr.ix]
  
  # overlap with memory information
  load(decision.file)
  m <- match(sprintf("%s.v3.1", curr.nodes), rownames(decisions))
  wcol1 <- which(colnames(decisions)==sprintf("group4_S_%s-group2_S_%s", tname, tname))
  wcol2 <- which(colnames(decisions)==sprintf("group4_R_%s-group2_R_%s", tname, tname))
  # classifyPattern for specified tissue (tname)
  a <- apply(decisions[m, c(wcol1, wcol2)], 1, classifyPattern)
  
  # plot
  library(igraph)
  ga <- graph_from_adjacency_matrix(adj.matrix.reduced, mode="undirected")
  ## font color according to memory
  
  add_shape("ellipse", clip=shapes("circle")$clip, plot=myellipse)
  patternCol <- c("black", "darkred", "lightblue", "red", "blue", "coral", "darkblue")
  #browser()
  #id <- tkplot(ga, canvas.width = 900, canvas.height = 310, vertex.shape="ellipse", vertex.color="white", edges.curved=F, vertex.size=15, vertex.label=curr.nodes, layout=layout_nicely, rescale=T, asp=0.3, vertex.label.family="sans", vertex.label.cex=1.7, vertex.label.color=patternCol[a+1], edge.width=1.3, edge.color=rep("black",nrow(ends(ga, E(ga), names=F))))
  ## can be changed interactively
  #my.coords <- tk_coords(id)
  #save(my.coords, file=sprintf("%scoords.Rdata", prefix))
  load(sprintf("%scoords.Rdata", prefix)) # my.coords
  my.vc <- rep("white", length(curr.nodes))
  mv <- match(my.genes, curr.nodes)
  my.vc[mv] <- my.color
  pdf(sprintf("%s_igraphMarked.pdf", prefix), width=21, height=7)
  plot(ga, vertex.shape="ellipse", vertex.color=my.vc, edges.curved=F, vertex.size=15, vertex.label=curr.nodes, layout=my.coords, rescale=T, asp=0.3, vertex.label.family="sans", vertex.label.cex=1.7, vertex.label.color=patternCol[a+1], edge.width=1.3, edge.color=rep("black",nrow(ends(ga, E(ga), names=F))))
  dev.off()
  
  # legend
  pdf(sprintf("%s_igraphlegend.pdf", prefix), width=21, height=7)
  plot(c(-1,1),c(-1,1))
  legend("center", legend=c("PP2C", "LEA4-5 co-ortholog (of 2)", "HB7 co-ortholog (of 4)", "NAC019 co-ortholog (of 2)", "MYB co-ortholog (of 2)"), fill=c("lightgoldenrod1", "cornsilk", "gray80", "darkseagreen1", "darkslategray1"), title="Gene annotation (Node color)", ncol=3, cex=1.7)
  dev.off()
  
  pdf(sprintf("%s_igraphlabellegend.pdf", prefix), width=21, height=7)
  plot(c(-1,1),c(-1,1))
  legend("center", legend=c("Up (R), up (S)", "Up (R), unchanged (S)", "Down (R), down (S)"), text.col=c("darkred", "red", "darkblue"), cex=1.7, ncol=2, title="Gene expression memory (Label color)", title.col="black")
  dev.off()
  
}

reverseTissue <- function(x){
  y <- x
  for (i in 1:nrow(x)){
    y[i,] <- x[nrow(x)+1-i,]
  }
  z <- y
  for (i in 1:ncol(y)){
    z[,i] <- y[, ncol(y)+1-i]
  }
  return(z)
}

runStressComparison <- function(infile, plotfile, savefile, mycex=1.7, myx=c(2,2), myy=c(42,21)) {
  load(infile) # decisions
  tissues <- c("J", "A", "P", "X", "W")
  phases <- c("S", "R")
  ntissues <- length(tissues)
  ngroups <- 3
  j <- 1
  barinfo.S.down <- matrix(0, ntissues, ngroups)
  barinfo.S.up <- matrix(0, ntissues, ngroups)
  for (i in 1:ntissues) {
    wcomp1 <- which(colnames(decisions)==sprintf("group4_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    wcomp2 <- which(colnames(decisions)==sprintf("group1_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    S.PS <- rownames(decisions)[which(decisions[, wcomp1]==1)]
    S.CS <- rownames(decisions)[which(decisions[, wcomp2]==1)]
    barinfo.S.up[i,] <- c(length(setdiff(S.PS, S.CS)), length(intersect(S.PS, S.CS)), length(setdiff(S.CS, S.PS)))
    S.PS <- rownames(decisions)[which(decisions[, wcomp1]==-1)]
    S.CS <- rownames(decisions)[which(decisions[, wcomp2]==-1)]  
    barinfo.S.down[i,] <- c(length(setdiff(S.PS, S.CS)), length(intersect(S.PS, S.CS)), length(setdiff(S.CS, S.PS)))
  }
  j <- 2
  barinfo.R.down <- matrix(0, ntissues, ngroups)
  barinfo.R.up <- matrix(0, ntissues, ngroups)
  for (i in 1:ntissues) { 
    wcomp1 <- which(colnames(decisions)==sprintf("group4_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    wcomp2 <- which(colnames(decisions)==sprintf("group1_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    R.PS <- rownames(decisions)[which(decisions[, wcomp1]==1)]
    R.CS <- rownames(decisions)[which(decisions[, wcomp2]==1)]
    barinfo.R.up[i,] <- c(length(setdiff(R.PS, R.CS)), length(intersect(R.PS, R.CS)), length(setdiff(R.CS, R.PS)))
    R.PS <- rownames(decisions)[which(decisions[, wcomp1]==-1)]
    R.CS <- rownames(decisions)[which(decisions[, wcomp2]==-1)]
    barinfo.R.down[i,] <- c(length(setdiff(R.PS, R.CS)), length(intersect(R.PS, R.CS)), length(setdiff(R.CS, R.PS)))
  }
  #print(barinfo.S.down)
  #print(barinfo.S.up)
  #print(barinfo.R.down)
  #print(barinfo.R.up)
  tissues.resorted <- mapTissue(tissues) 
  save(barinfo.S.down, barinfo.S.up, barinfo.R.down, barinfo.R.up, tissues.resorted, ntissues, file=savefile)
  
  ## bar chart with same axis scales
  global.max <- max(max(max(barinfo.S.down)), max(max(barinfo.S.up)), max(max(barinfo.R.down)), max(max(barinfo.R.up)))
  pdf(plotfile, width=9, height=7)
  par(mfrow=c(2,2), mar=c(5, 4, 2, 1.6))  
  tick.positions <- seq(0,global.max,500)
  ## since horizontal barplot goes on y-axis from bottom to top, need to reverse tissue order and condition order
  barplot((-1)*t(reverseTissue(barinfo.S.down)), beside=T, col=rep(rev(c("lightcoral", "black", "lightgoldenrod1")), ntissues), horiz=T, axes=F, axisnames=F, xlim=c((-1)*tick.positions[length(tick.positions)], tick.positions[1]), xaxt="n")
  axis(side=1,at=(-1)*rev(tick.positions),labels=rev(tick.positions), cex.axis=mycex)
  
  bp <- barplot(t(reverseTissue(barinfo.S.up)), beside=T, col=rep(rev(c("lightcoral", "black", "lightgoldenrod1")), ntissues), names.arg=rev(tissues.resorted), horiz=T, las=1, axes=F, xlim=c(tick.positions[1], tick.positions[length(tick.positions)]), xaxt="n", cex.names=mycex)
  axis(side=1,at=tick.positions,labels=tick.positions, cex.axis=mycex)
  text(myx[1], myy[1], labels=c("Stress phase"), cex=mycex, xpd=T, font=1, pos=4)
  
  barplot((-1)*t(reverseTissue(barinfo.R.down)), beside=T, col=rep(rev(c("lightcoral", "black", "lightgoldenrod1")), ntissues), horiz=T, axes=F, axisnames=F, xlim=c((-1)*tick.positions[length(tick.positions)], tick.positions[1]), xaxt="n", xlab="Number of down-regulated genes", cex.lab=mycex)
  axis(side=1,at=(-1)*rev(tick.positions),labels=rev(tick.positions), cex.axis=mycex) 
  
  barplot(t(reverseTissue(barinfo.R.up)), beside=T, col=rep(rev(c("lightcoral", "black", "lightgoldenrod1")), ntissues), legend.text=rev(c("PS-specific", "Overlap", "CS-specific")), args.legend=list(x="topright", cex=mycex), names.arg=rev(tissues.resorted), horiz=T, las=1, axes=F, xlim=c(tick.positions[1], tick.positions[length(tick.positions)]), xaxt="n", cex.names=mycex, xlab="Number of up-regulated genes", cex.lab=mycex)
  axis(side=1,at=tick.positions,labels=tick.positions, cex.axis=mycex) 
  text(myx[2], myy[2], labels=c("Recovery phase"), cex=mycex, xpd=T, font=1, pos=4) 
  dev.off()
}

runStressComparisonZoom <- function(infile, plotfile, mycex=1.7, myx=c(2,2), myy=c(42,21)) {
  load(infile) # barinfo.S.down, barinfo.S.up, barinfo.R.down, barinfo.R.up, tissues.resorted, ntissues
    
  ## bar chart with same axis scales
  global.max <- max(max(max(barinfo.S.down)), max(max(barinfo.S.up)), max(max(barinfo.R.down)), max(max(barinfo.R.up)))
  pdf(plotfile, width=5.8, height=7) 
  par(mfrow=c(2,2), mar=c(5, 4, 2, 1.6))  
  tick.positions <- seq(0,100,20) 
 
  barplot((-1)*t(reverseTissue(barinfo.R.down)), beside=T, col=rep(rev(c("lightcoral", "black", "lightgoldenrod1")), ntissues), horiz=T, axes=F, axisnames=F, xlim=c((-1)*tick.positions[length(tick.positions)], tick.positions[1]), xaxt="n", xlab="# Down-regulated genes", cex.lab=mycex)
  axis(side=1,at=(-1)*rev(tick.positions[1:2]),labels=rev(tick.positions[1:2]), cex.axis=mycex) 
  
  barplot(t(reverseTissue(barinfo.R.up)), beside=T, col=rep(rev(c("lightcoral", "black", "lightgoldenrod1")), ntissues), names.arg=rev(tissues.resorted), horiz=T, las=1, axes=F, xlim=c(tick.positions[1], tick.positions[length(tick.positions)]), xaxt="n", cex.names=mycex, xlab="# Up-regulated genes", cex.lab=mycex) 
  axis(side=1,at=tick.positions,labels=tick.positions, cex.axis=mycex) 
  text(myx[2], myy[2], labels=c("Recovery"), cex=mycex, xpd=T, font=2, pos=4) 
  dev.off()
}

overlapPSandCS <- function(infile, j=2) { # by default only for recovery
  load(infile) # decisions
  tissues <- c("J", "A", "P", "X", "W")
  phases <- c("S", "R")
  ntissues <- length(tissues)
  ngroups <- 3
  barinfo.R.down <- matrix(0, ntissues, ngroups)
  barinfo.R.up <- matrix(0, ntissues, ngroups)
  res <- list()
  for (i in 1:ntissues) { 
    wcomp1 <- which(colnames(decisions)==sprintf("group4_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    wcomp2 <- which(colnames(decisions)==sprintf("group1_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    R.PS <- rownames(decisions)[which(decisions[, wcomp1]==1)]
    R.CS <- rownames(decisions)[which(decisions[, wcomp2]==1)]
    up <- list(setdiff(R.PS, R.CS), intersect(R.PS, R.CS), setdiff(R.CS, R.PS))
    R.PS <- rownames(decisions)[which(decisions[, wcomp1]==-1)]
    R.CS <- rownames(decisions)[which(decisions[, wcomp2]==-1)]
    down <- list(setdiff(R.PS, R.CS), intersect(R.PS, R.CS), setdiff(R.CS, R.PS))
    res[[i]] <- list(up, down)
  }
  return(res)
}

matchAndSort <- function(s, genes) {
  rescol.raw <- rep("white", length(genes))
  m1 <- match(s[[1]], genes)
  m2 <- match(s[[2]], genes)
  m3 <- match(s[[3]], genes)
  rescol.raw[m1] <- "lightcoral"
  rescol.raw[m2] <- "black"
  rescol.raw[m3] <- "lightgoldenrod1"
  resix <- c(m1, m3, m2)
  rescol <- rescol.raw[resix]
  return(list(resix, rescol))
}

plotVolcano <- function(sets, compfiles, short, outdir, mycex=1.7) {
  # plot up and down separately
  # circles with black border; colors according to Fig. 2C (comparison vs. EC)
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  reg <- c("Up", "Down")
  sc <- c(0.49, 0.51, 0.49, 0.5, 0.5)
  for (i in 1:length(compfiles)){
    load(compfiles[i]) # comparison
    for (j in 1:length(reg)){
      sortres <- matchAndSort(sets[[i]][[j]], rownames(comparison))
      myix <- sortres[[1]]
      pdf(sprintf("%s%s_%s.pdf", outdir, short[i], reg[j]))
      par(mar=c(6.7,9.6,4.1,2.1)) 
      plot(comparison[myix,"log2FoldChange"], -log10(comparison[myix,"padj"]), pch=21, bg=sortres[[2]], cex=mycex, cex.axis=mycex, cex.lab=mycex, xlab=expression(paste('log'[2], "(fold change)")), ylab=expression(paste('-log'[10], "(adjusted p-value)")))
      abline(h=-log10(0.05))
      abline(v=1)
      abline(v=-1)
      if (j==2){
        curr.tissue <- extractPart(short[i], part=2, spl="_")
        myx <- min(comparison[myix,"log2FoldChange"], na.rm=T)
        text(myx+sc[i]*myx, max(-log10(comparison[myix,"padj"]), na.rm=T)/2, labels=c(mapTissue(curr.tissue)), cex=mycex, xpd=T, font=1, pos=2)
        
      }
      if (i==1){
        myy <- max(-log10(comparison[myix,"padj"]), na.rm=T)
        text(0, myy+0.2*myy, labels=c(reg[j]), cex=mycex, xpd=T, font=1, pos=1)
      }
      dev.off()
    }
  }
}

mergeSets <- function(s1,s2) {
  return(list(union(s1[[1]], s2[[1]]), union(s1[[2]], s2[[2]]), union(s1[[3]], s2[[3]])))
}

plotVolcanoCombined <- function(sets, compfiles, short, outdir, mycex=1.7, mylabel, printTissue) {
  # plot up and down in one figure
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  sc <- c(0.49, 0.51, 0.49, 0.5, 0.5)
  for (i in 1:length(compfiles)){
    load(compfiles[i]) # comparison
    newsets <- mergeSets(sets[[i]][[1]], sets[[i]][[2]])
    sortres <- matchAndSort(newsets, rownames(comparison))
    myix <- sortres[[1]]
    pdf(sprintf("%s%s.pdf", outdir, short[i]))
    par(mar=c(6.7,9.6,4.1,2.1))
    plot(comparison[myix,"log2FoldChange"], -log10(comparison[myix,"padj"]), pch=21, bg=sortres[[2]], cex=mycex, cex.axis=mycex, cex.lab=mycex, xlab=expression(paste('log'[2], "(fold change)")), ylab=expression(paste('-log'[10], "(adjusted p-value)")))
    abline(h=-log10(0.05))
    abline(v=1)
    abline(v=-1)
    if (printTissue){
      curr.tissue <- extractPart(short[i], part=2, spl="_")
      myx <- min(comparison[myix,"log2FoldChange"], na.rm=T)
      text(myx+sc[i]*myx, max(-log10(comparison[myix,"padj"]), na.rm=T)/2, labels=c(mapTissue(curr.tissue)), cex=mycex, xpd=T, font=1, pos=2)
      
    }
    if (i==1){
      myy <- max(-log10(comparison[myix,"padj"]), na.rm=T)
      text(0, myy+0.2*myy, labels=c(mylabel), cex=mycex, xpd=T, font=1, pos=1)
    }
    dev.off()
  }
}

volcanoWrapper <- function (infile, res.dir, outdir) {
  if (!file.exists(outdir)){
    dir.create(outdir)
  }
  # get sets from Fig 2C, three sets per tissue and direction of regulation (PS-specific, overlap, CS-specific)
  my.sets <- overlapPSandCS(infile)
  
  phases <- c("R")
  tissues <- c("J","A", "P", "X", "W") # must be same order as in overlapPSandCS!
  
  # show sets in direct PS-CS comparison
  treatment.types.out <- c("PS")
  my.files <- c()
  my.short <- c()
  for (i in 1:length(tissues)){
    for (j in 1:length(phases)){
      for (k in 1:length(treatment.types.out)){
        my.files <- c(my.files, sprintf("%scomparison_%s_%s_%s_vsCS.Rdata", res.dir, tissues[i], phases[j], treatment.types.out[k]))
        my.short <- c(my.short, sprintf("volc_%s", tissues[i]))
      }
    }
  }
  plotVolcano(my.sets, my.files, my.short, sprintf("%sPSvsCS/", outdir))
  
  # show sets in PS vs. EC comparison
  treatment.types.out <- c("PS")
  my.files <- c()
  my.short <- c()
  for (i in 1:length(tissues)){
    for (j in 1:length(phases)){
      for (k in 1:length(treatment.types.out)){
        my.files <- c(my.files, sprintf("%scomparison_%s_%s_%s_vsEC.Rdata", res.dir, tissues[i], phases[j], treatment.types.out[k]))
        my.short <- c(my.short, sprintf("volc_%s", tissues[i]))
      }
    }
  }
  plotVolcanoCombined(my.sets, my.files, my.short, sprintf("%sPSvsEC/", outdir), mylabel="PS", printTissue=T)
  
  # show sets in CS vs. EC comparison
  treatment.types.out <- c("CS")
  my.files <- c()
  my.short <- c()
  for (i in 1:length(tissues)){
    for (j in 1:length(phases)){
      for (k in 1:length(treatment.types.out)){
        my.files <- c(my.files, sprintf("%scomparison_%s_%s_%s_vsEC.Rdata", res.dir, tissues[i], phases[j], treatment.types.out[k]))
        my.short <- c(my.short, sprintf("volc_%s", tissues[i]))
      }
    }
  }
  plotVolcanoCombined(my.sets, my.files, my.short, sprintf("%sCSvsEC/", outdir), mylabel="CS", printTissue=F)
  
  # plot legend
  mycex <- 1.7
  pdf(sprintf("%slegend.pdf", outdir), width=14, height=7)
  plot(0,0)
  legend(-0.5, 0.5, legend=c("PS-specific", "Overlapping", "CS-specific"), pt.bg=c("lightcoral", "black", "lightgoldenrod1"), pch=21, title="Differentially expressed recovery genes relative to EC", ncol=3, cex=mycex, pt.cex=mycex) 
  dev.off()
  
}

reverseTissue2 <- function(x){
  y <- x
  for (i in 1:nrow(x)){
    y[i,] <- x[nrow(x)+1-i,]
  }
  z <- y
  return(z)
}

plotTransition <- function (infile, outfile, mycex=1.7) { 
  load(infile) # decisions
  tissues <- c("J", "A", "P", "X", "W")
  phases <- c("S", "R")
  ntissues <- length(tissues)
  ngroups <- 3
  
  barinfo <- matrix(0, 2*ntissues, ngroups)
  
  for (i in 1:ntissues) {
    j <- 2
    wcomp1 <- which(colnames(decisions)==sprintf("group4_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    wcomp2 <- which(colnames(decisions)==sprintf("group1_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    PS.R <- which(decisions[, wcomp1]==1)
    CS.R <- which(decisions[, wcomp2]==1)
    j <- 1
    wcomp1.S <- which(colnames(decisions)==sprintf("group4_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    wcomp2.S <- which(colnames(decisions)==sprintf("group1_%s_%s-group2_%s_%s", phases[j], tissues[i], phases[j], tissues[i]))
    barinfo[(i-1)*2+1,1] <- 100*length(which(decisions[PS.R, wcomp1.S]==-1))/length(PS.R)
    barinfo[(i-1)*2+1,3] <- 100*length(which(decisions[PS.R, wcomp1.S]==1))/length(PS.R)
    barinfo[(i-1)*2+1,2] <- 100-(barinfo[(i-1)*2+1,1]+barinfo[(i-1)*2+1,3])
    
    barinfo[(i-1)*2+2,1] <- 100*length(which(decisions[CS.R, wcomp2.S]==-1))/length(CS.R)
    barinfo[(i-1)*2+2,3] <- 100*length(which(decisions[CS.R, wcomp2.S]==1))/length(CS.R)
    barinfo[(i-1)*2+2,2] <- 100-(barinfo[(i-1)*2+2,1]+barinfo[(i-1)*2+2,3])
  }
  #print(barinfo)
  
  # Barplot overlap of recovery and stress phase regulation
  tissues.resorted <- c("LE1", "LE2", "PHL", "XYL", "ROO")
  pdf(outfile, width=7.8, height=7) 
  par(mar=c(8,1.1,0.01,2.1)) 
  barplot(t(reverseTissue2(barinfo)), beside=F, col=rep(c("blue", "white", "red"), nrow(barinfo)), names.arg=NA, horiz=T, las=1, axes=F, axisnames=F, xaxt="n", cex.lab=mycex, space=rep(c(0.4,0.2), ntissues), xlim=c(-12,111.2), width=0.05, ylim=c(-0.1,0.7) , legend.text=c("Down (%)", "Unchanged (%)", "Up (%)"), args.legend=list(x="bottom", ncol=2, title="Stress phase regulation", cex=mycex, inset=c(-0.2,-0.27), text.width=41.5, x.intersp=0.2, y.intersp=1)) 
  text(rep(-8, ntissues), 0.065*(seq(1.5, 10, 2)+0.15), rev(tissues.resorted), pos=1, cex=mycex, xpd=T) # coordinates must be within xlim to be visible
  text(rep(100, 2*ntissues), 0.065*(c(1:(2*ntissues))-0.42), rep(c("CS", "PS"), ntissues), pos=4, cex=mycex)
  title(xlab="Up-regulated genes in recovery phase", cex.lab=mycex, line=-3.5)
  dev.off()
}

writeGraphviz <- function(edges, graphfile){
  write("graph G {\n node [shape=ellipse fixedsize=false height=0.1 fontsize=14 fontname=\"arial\" style=filled fillcolor=white];", file=graphfile)
  for (i in 1:nrow(edges)){
    write(sprintf("%s -- %s \n", edges[i,1], edges[i,2]), file=graphfile, append=T)
  }
  write("}\n", file=graphfile, append=T)
}

collectPPI <- function(annfile, tfFile, ppifile1, ppifile2, ppifile3, outfile, dotfile){
  atf.raw <- read.delim(tfFile, check.names=F, stringsAsFactors=F)
  araCol <- 11
  rd <- read.delim(annfile, header=T, stringsAsFactors=F, check.names=F)
  wm <- which(atf.raw[,1] %in% rd[, araCol])
  atf <- unique(extractPart(atf.raw[wm,1], part=1, spl="."))
  rp1 <- read.delim(ppifile1, stringsAsFactors=F, check.names=F)
  w1 <- which(rp1[,1] %in% atf & rp1[,3] %in% atf)
  rp2 <- read.delim(ppifile2, stringsAsFactors=F, check.names=F)
  w2 <- which(rp2[,1] %in% atf & rp2[,2] %in% atf)
  rp3 <- read.delim(ppifile3, stringsAsFactors=F, check.names=F)
  w3 <- which(rp3[,1] %in% atf & rp3[,3] %in% atf)
  all.edges <- rbind(as.matrix(rp1[w1,c(1,3)]), as.matrix(rp2[w2, c(1,2)]), as.matrix(rp3[w3,c(1,3)]))
  writeGraphviz(all.edges, dotfile) # can be opened with Graphviz software
  save(all.edges, file=outfile)
}

makePPIgraph <- function(efile, decisionfile, prefix) {
  load(efile) # all.edges
  ppinodes <- c("AT2G46680", "AT4G27410", "AT5G13180", "AT5G65210", "AT4G24660", "AT5G15210", "AT1G69600", "AT3G28920", "AT1G14440", "AT2G02540", "AT1G75240", "AT2G18350", "AT1G32640", "AT1G66350", "AT5G05790") 
  ## name according to ortholog group to which also poplar gene belongs
  ppinodes.names <- c("HB7", "NAC019", "NAC083", "TGA1", "ZHD1", "ZHD10/8", "ZHD10/11", "ZHD10/9", "ZHD3/4", "ZHD3", "ZHD5", "ZHD6", "MYC2", "RGL1", "MYB") 
  w <- which(all.edges[,1] %in% ppinodes & all.edges[,2] %in% ppinodes)
  ppiedges.raw <- all.edges[w,]
  ppiedges.raw.df <- data.frame(pcor=rep(1, nrow(ppiedges.raw)), node1=ppiedges.raw[,1], node2=ppiedges.raw[,2], stringsAsFactors=F)
  is <- intersectEdges(ppiedges.raw.df, ppiedges.raw.df) # to make unique; includes sortNodes
  ppiedges <- cbind(is$node1, is$node2)
  
  adj.matrix <- matrix(0, length(ppinodes), length(ppinodes))
  m1 <- match(ppiedges[,1], ppinodes)
  m2 <- match(ppiedges[,2], ppinodes)
  adj.matrix[cbind(m1,m2)] <- 1
  adj.matrix[cbind(m2,m1)] <- 1
  
  library(igraph)
  ga <- graph_from_adjacency_matrix(adj.matrix, mode="undirected")
  add_shape("ellipse", clip=shapes("circle")$clip, plot=myellipse2)
  
  ppinodes.poplar <- list(c("Potri.014G103000", "Potri.001G083700"), c("Potri.011G123300", "Potri.001G404100"), c("Potri.001G061200", "Potri.003G166500"), c("Potri.007G079900", "Potri.007G085700", "Potri.005G082000"), c("Potri.002G102900", "Potri.005G158800"), c("Potri.004G126500"), c("Potri.008G086000", "Potri.010G169400"), c("Potri.017G082900"), c("Potri.003G000400", "Potri.004G229600", "Potri.003G146700"), c("Potri.004G135100"), c("Potri.005G227900", "Potri.002G035200"), c("Potri.005G122500", "Potri.007G024100"), c("Potri.001G142200", "Potri.003G092200"), c("Potri.001G326000", "Potri.005G095100", "Potri.010G143200", "Potri.010G143400", "Potri.010G143700", "Potri.012G093900", "Potri.015G091200"), c("Potri.010G193000"))
  # overlap with memory information
  ppinodes.color <- rep("white", length(ppinodes)) # white: constitutive (see TPM_filtered.Rdata); gray: multi-tissue memory; "darkseagreen1": tissue-specific memory
  load(decisionfile)
  tname <- "A"
  wcol1 <- which(colnames(decisions)==sprintf("group4_R_%s-group2_R_%s", tname, tname))
  tname <- "X"
  wcol2 <- which(colnames(decisions)==sprintf("group4_R_%s-group2_R_%s", tname, tname))
  for (i in 1:length(ppinodes.poplar)) {
    m <- match(sprintf("%s.v3.1", ppinodes.poplar[[i]]), rownames(decisions))
    wa <- which(decisions[m,wcol1]!=0)
    wx <- which(decisions[m,wcol2]!=0)
    wis <- intersect(wa, wx)
    if (length(wis)>0 ){
      ppinodes.color[i] <- "gray80"
    }
    else {
      if (length(wa)>0 || length(wx>0)){
        ppinodes.color[i] <- "darkseagreen1"
      }
    }
    
  }
  
  #browser() # interactive graph plotting
  #id <- tkplot(ga, canvas.width = 900, canvas.height = 310, vertex.shape="ellipse", vertex.color=ppinodes.color, edges.curved=F, vertex.size=15, vertex.label=ppinodes.names, layout=layout_nicely, rescale=T, asp=0.3, vertex.label.family="sans", vertex.label.cex=1.7, vertex.label.color="black", edge.width=1.3, edge.color=rep("black",nrow(ends(ga, E(ga), names=F))))
  #my.coords <- tk_coords(id)
  #save(my.coords, file=sprintf("%scoords.Rdata", prefix))
  load(sprintf("%scoords.Rdata", prefix))
  pdf(sprintf("%s_igraph.pdf", prefix), width=14)
  plot(ga, vertex.shape="ellipse", vertex.color=ppinodes.color, edges.curved=F, vertex.size=7, vertex.label=ppinodes.names, layout=my.coords, rescale=T, asp=0.55, vertex.label.family="sans", vertex.label.cex=1.7, vertex.label.color="black", edge.width=1.3, edge.color=rep("black",nrow(ends(ga, E(ga), names=F))))
  dev.off()
  
  pdf(sprintf("%s_igraphlegend.pdf", prefix), width=21, height=7)
  plot(c(-1,1),c(-1,1))
  legend("center", legend=c("Multi-tissue memory", "Tissue-specific memory", "Constitutive"), fill=c("gray80", "darkseagreen1", "white"), title="Expression", ncol=3, text.width=0.32, cex=1.7)
  dev.off()
}
