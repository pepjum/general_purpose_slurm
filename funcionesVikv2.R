######################################################################################
# FUNCIONES
######################################################################################

library(geneplotter)
library(genefilter)
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
# library(gregmisc) #trim
library(maxstat)
library(survival)
library(limma)
library(stringr)

#######################################
##
## Quality control
##
#######################################

"multi" = function(type, x, xlim, title1, title2, title3, ...) {

    if (type == "density")
      multidensity(x,
                   xlim = xlim,
                   main = "",
                   xlab = "",
                   ylab = "", ...)
    if (type == "ecdf")
      multiecdf(x,
                xlim = xlim,
                main = "",
                xlab = "",
                ylab = "", ...)
    mtext(title1, side = 2, adj = 0.5, padj = -4 , cex = 0.7)
    mtext(title2, side = 1, adj = 0.5, padj = 4 , cex = 0.7)
    mtext(title3, side = 3, adj = 0.5, padj = -1 , cex = 0.7)
     
}


#################################################
##
## Estudio de VARIABILIDAD (Outliers)
##
#################################################

"detOutliers" <- function(data, typeSamples, filenamePdf, W, H) { 

	colourRange = rgb(seq(0,1,l=256),seq(0,1,l=256),seq(1,0,l=256))
	
	sArray = apply(data, 1, median)
	dat <- data
	for (i in 1:ncol(data))
		dat[,i] <- dat[,i]-sArray
		
	outM = as.dist(dist2(na.omit(dat)))
		
	d.row = as.dendrogram(hclust(outM))
	od.row = order.dendrogram(d.row)
	m = as.matrix(outM)
	namesCol <- colnames(data)
	colnames(m) = namesCol
	rownames(m) = namesCol
   
	covar = typeSamples
	lev = levels(as.factor(covar))
	corres = matrix(0,nrow=length(lev),ncol=2)
	colourCov = brewer.pal(12,"Set3")

	pdf(file = filenamePdf, width = W, heigh = H, colormode = "rgb")
		print(levelplot(m[od.row,od.row],
			scales=list(x=list(rot=90)),
			legend=list(
			top=list(fun=dendrogramGrob,args=list(x=d.row,side="top")),
			right=list(fun=dendrogramGrob,args=list(x=d.row,side="right", size.add= 1, add = list(rect = list(col = "transparent", fill = colourCov[as.factor(covar)])), type = "rectangle"))),
			colorkey = list(space ="left"),
			xlab="",ylab="",
			col.regions=colourRange))
		
		x=0.06
		y=0.98
		
		for(i in 1:length(lev))
		{
		corres[i,] = c(unique(covar[covar == lev[i]]),colourCov[i])
		grid.text(lev[i], x=x, y=y, just="left")
		grid.rect(gp=gpar(fill=corres[i,2],col="transparent"), x=x-0.02, y=y, width=0.02, height=0.02)
		y=y-0.03
		}
	dev.off()

	return(outM)
}


#############################################
##
## ANOTAR sondas del microarray
##
#############################################

#Para arrays de affy 3'UTR
"annotVik" <- function(data, array) {
	data.annot = data.frame(paste(data$ID), paste(lookUp(paste(data$ID), array, "SYMBOL")), paste(lookUp(paste(data$ID), array, "GENENAME")), paste(lookUp(paste(data$ID), array, "MAP")), data);
	colnames(data.annot)[1] = "Probeset";
	colnames(data.annot)[2] = "Name";
	colnames(data.annot)[3] = "Description";
	colnames(data.annot)[4] = "Locus";
	
	data.annot
}

#Para arrays de affy gene st
"parseAnnotGene" <- function(annot) {
	annG <- sapply(annot[,1],FUN=function(x) paste(unique(str_trim(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
	annT <- sapply(annot[,2],FUN=function(x) paste(unique(str_trim(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
	annG_m <- matrix(unlist(sapply(annG, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} )),nrow(annot), 3, byrow = T)
	annT_m <- matrix(paste(unlist(sapply(annT, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} ))),nrow(annot), 3, byrow = T)
	ann <- cbind(annG_m,annT_m)
	colnames(ann) <- c("GeneID", "GeneName", "GeneDescription", "TranscriptID", "DBName", "TranscriptDescription")
	ann
}

#Para arrays de Mouse affy gene 2.0 st
"parseAnnotGene_MoGene20" <- function(annot) {
        annG <- sapply(annot[,2],FUN=function(x) paste(unique(gdata::trim(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
        annT <- sapply(annot[,3],FUN=function(x) paste(unique(gdata::trim(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
	annG_m <- matrix(unlist(sapply(annG, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} )),nrow(annot), 3, byrow = T)
        annT_m <- matrix(paste(unlist(sapply(annT, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} ))),nrow(annot),3, byrow = T)
        ann <- data.frame(annot[,1], annG_m, annT_m, paste(annot[,4]))
	ann <- ann[ann[,8] == "main", -8]
        colnames(ann) <- c("Probeset_id", "GeneID", "GeneName", "GeneDescription", "TranscriptID", "DBName", "TranscriptDescription")
	rownames(ann) <- ann[,1]
	ann
} 

#Para arrays Human affy gene 2.0 st
"parseAnnotGene_HuGene20" <- function(annot) {
        annG <- sapply(annot[,2],FUN=function(x) paste(unique(trimWhiteSpace(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
        annT <- sapply(annot[,3],FUN=function(x) paste(unique(trimWhiteSpace(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
	annG_m <- matrix(unlist(sapply(annG, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} )),nrow(annot), 3, byrow = T)
        annT_m <- matrix(paste(unlist(sapply(annT, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} ))),nrow(annot),3, byrow = T)
        ann <- data.frame(annot[,1], annG_m, annT_m, paste(annot[,4]))
	ann <- ann[(ann[,8] == "main" | ann[,8] == "lncrna"), -8]
        colnames(ann) <- c("Probeset_id", "GeneID", "GeneName", "GeneDescription", "TranscriptID", "DBName", "TranscriptDescription")
	rownames(ann) <- ann[,1]
	ann
} 

#Para arrays HTA 2.0 st
"parseAnnotGene_HTA20" <- function(annot) {
        annG <- sapply(annot[,2],FUN=function(x) paste(unique(trimWhiteSpace(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
        annT <- sapply(annot[,3],FUN=function(x) paste(unique(trimWhiteSpace(unlist(strsplit(paste(x), "//"))))[1:3], collapse = "//"))
	annG_m <- matrix(unlist(sapply(annG, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} )),nrow(annot), 3, byrow = T)
        annT_m <- matrix(paste(unlist(sapply(annT, FUN = function(x) { strsplit(paste(x), "//") [[1]][1:3]} ))),nrow(annot),3, byrow = T)
        ann <- data.frame(annot[,1], annG_m, annT_m, paste(annot[,4]), paste(annot[,5]))
	ann <- ann[(ann[,8] == "main" | ann[,8] == "lncrna"), -8]
        colnames(ann) <- c("Probeset_id", "GeneID", "GeneName", "GeneDescription", "TranscriptID", "DBName", "TranscriptDescription","LocusType")
	rownames(ann) <- ann[,1]
	ann
} 


# LEER GFF
"getAttributeField" <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

"gffRead" <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

#####################################################
##
## FUNCIONES PARA ANÁLISIS DE EXPRESIÓN DIFERENCIAL
##
#####################################################

# microarrays de affy con sondas en 3UTR
"graphContrast_3utr" <- function(data, name, ngenes) {
	hist(data$P.Val, main = name, xlab = "p-value", ylab = "Genes", col=NULL);
	plot(data$adj.P.Val,data$P.Val, type="l", xlab="fdr", ylab="p-value", col=NULL, main="FDR")
	x <- data$logFC
	y <- -log(data$P.Value)
	plot(x, y, xlab="Log2(Fold Change)", ylab = "-log(p-value)", main = name)
	abline(v=-1);abline(v=1);
	text (x[1:ngenes], y[1:ngenes], data[1:ngenes, 2],cex=0.7)
}
# nuevos microarrays
"graphContrast" <- function(data, name, Bth, FCth, namecol) {
	hist(data$P.Value, main = name, xlab = "p-value", ylab = "Genes", col=NULL);
	hist(data$logFC, main = name, xlab="Log2(FoldChange)", ylab="Genes",  col=NULL, 50);
#	hist(treatm_vs_ctrl$B, main = name, xlab="B", ylab="Genes", 50);
	volcanoCol(data, Bth, FCth, name, namecol)
}

"volcanoCol" <- function(res, Bth, FCth, title, namecol) {
	colVP <- rep("black", length(res$B))
	colVP[which(res$B>Bth & res$logFC>FCth)] <- "red"
	colVP[which(res$B>Bth & res$logFC<(FCth*(-1)))] <- "green"
	plot(res$logFC, (-1)*log10(res$P.Value), pch = ".", col = colVP, main = title, xlab = "foldchange", ylab = "-log(pvalue)")
	abline(v = FCth)
	abline(v = (FCth*(-1)))
	abline(h = (-1)*log10(max(res$P.Value[res$B>Bth])))
	selFC <- res$logFC[which(res$B>Bth & abs(res$logFC)>FCth)]
	colVP2 <- rep("red", length(selFC))
	colVP2[which(selFC<((-1)*FCth))] <- "green"
	if (length(res[which(res$B>Bth & abs(res$logFC)>FCth),namecol])>0)
		text(res$logFC[which(res$B>Bth & abs(res$logFC)>FCth)], (-1)*log10(res$P.Value[which(res$B>Bth & abs(res$logFC)>FCth)]), res[which(res$B>Bth & abs(res$logFC)>FCth),namecol], pos = 3, cex=0.7, offset = 0.5, col = colVP2)
} 

"volcanoColP" <- function(res, Bth, FCth, title, namecol) {
        colVP <- rep("black", length(res$P.Value))
        colVP[which(res$P.Value<Bth & res$logFC>FCth)] <- "red"
        colVP[which(res$P.Value<Bth & res$logFC<(FCth*(-1)))] <- "green"
        plot(res$logFC, (-1)*log10(res$P.Value), pch = ".", col = colVP, main = title, xlab = "foldchange", ylab = "-log(pvalue)")
        abline(v = FCth)
        abline(v = (FCth*(-1)))
        abline(h = (-1)*log10(max(res$P.Value[res$P.Value<Bth], na.rm=TRUE)))
        selFC <- res$logFC[which(res$P.Value<Bth & abs(res$logFC)>FCth)]
        colVP2 <- rep("red", length(selFC))
        colVP2[which(selFC<((-1)*FCth))] <- "green"
        if (length(res[which(res$P.Value<Bth & abs(res$logFC)>FCth),namecol])>0)
                text(res$logFC[which(res$P.Value<Bth & abs(res$logFC)>FCth)], (-1)*log10(res$P.Value[which(res$P.Value<Bth & abs(res$logFC)>FCth)]), res[which(res$P.Value<Bth & abs(res$logFC)>FCth),namecol], pos = 3, cex=0.7, offset = 0.5, col = colVP2)
}

"graphContrastP" <- function(data, name, Bth, FCth, namecol) {
        
	hist(data$P.Value, main = name, xlab = "p-value", ylab = "Genes", col=NULL);
    hist(data$logFC, main = name, xlab="Log2(FoldChange)", col=NULL, ylab="Genes", 50);
    volcanoColP(data, Bth, FCth, name, namecol)
}

"compare2List" <- function(listA, listB, nameA, nameB, title) {
	
	listu <- unique(c(listA, listB))
	vennC <- cbind((listu %in% listA)*1, (listu %in% listB)*1)
	colnames(vennC) <- c(nameA, nameB)
	vennDiagram(vennC, main = title)

}

"compare3List" <- function(listA, listB, listC, nameA, nameB, nameC, title) {
	
	listu <- unique(c(listA, listB, listC))
	vennC <- cbind((listu %in% listA)*1, (listu %in% listB)*1, (listu %in% listC)*1)
	colnames(vennC) <- c(nameA, nameB, nameC)
	vennDiagram(vennC, main = title)
}

"compare4List" <- function(listA, listB, listC, listD, nameA, nameB, nameC, nameD, title) {
	
	listu <- unique(c(listA, listB, listC, listD))
	vennC <- cbind((listu %in% listA)*1, (listu %in% listB)*1, (listu %in% listC)*1, (listu %in% listD)*1)
	colnames(vennC) <- c(nameA, nameB, nameC, nameD)
	vennDiagram(vennC, main = title)

}

# usa vennDiagram de limma no nuestro custom vennDiagram
"compare5List" <- function(listA, listB, listC, listD, listE, nameA, nameB, nameC, nameD, nameE, title) {
	
	listu <- unique(c(listA, listB, listC, listD, listE))
	vennC <- cbind((listu %in% listA)*1, (listu %in% listB)*1, (listu %in% listC)*1, (listu %in% listD)*1, (listu %in% listE)*1)
	colnames(vennC) <- c(nameA, nameB, nameC, nameD, nameE)
	limma::vennDiagram(vennC, circle.col=c("blue","green","orange","purple","red"),main = title)

}

#############################################
##
## FUNCIONES PARA ANALISIS DE GO
##
#############################################

"goVik" <- function(selProbes, refProbes, envSyn, hgCutoff, annot, catSize) {

	selectedEntrez<- unlist(mget(paste(selProbes), env = envSyn, ifnotfound = NA))
	geneUniverse <-unlist(mget(refProbes, env = envSyn, ifnotfound = NA))
	selectedEntrez <-selectedEntrez[!is.na(selectedEntrez)]
	geneUniverse <-geneUniverse[!is.na(geneUniverse)]

	params <- new("GOHyperGParams",
		geneIds=selectedEntrez,
		universeGeneIds=geneUniverse,
		annotation= annot,
		ontology="BP",
		pvalueCutoff=hgCutoff,
		conditional=TRUE,
		testDirection="over")		
	hgOverBP <- hyperGTest(params)
	resultGOBP <- resultGO <- summary(hgOverBP, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);
	
	ontology(params) <- "MF";
	hgOverMF <- hyperGTest(params)
	resultGOMF <- resultGO <- summary(hgOverMF, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);
	
	ontology(params) <- "CC";
	hgOverCC <- hyperGTest(params)
	resultGOCC <- resultGO <- summary(hgOverCC, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);
	
	return (res = list(resultGOBP, resultGOMF, resultGOCC))
	
}
# 
# 
"geneGOMat" <- function(accSel, accRef, envSyn, hgCutoff, annot, catSize,  catgo = "ALL", formatMat = FALSE) {

	acc_GO <- goVik(accSel, accRef, envSyn, hgCutoff, annot, catSize)
	
	if (catgo == "ALL") {
		acc_GOCodes <-  c(acc_GO[[1]]$GOBPID, acc_GO[[2]]$GOMFID, acc_GO[[3]]$GOCCID)
		acc_GOTerms <-  c(acc_GO[[1]]$Term, acc_GO[[2]]$Term, acc_GO[[3]]$Term)
		acc_GOPvalue <-  c(acc_GO[[1]]$Pvalue, acc_GO[[2]]$Pvalue, acc_GO[[3]]$Pvalue)
		acc_GOCount <-  c(acc_GO[[1]]$Count, acc_GO[[2]]$Count, acc_GO[[3]]$Count)
	} else if (catgo == "BP") {
		acc_GOCodes <- acc_GO[[1]]$GOBPID
		acc_GOTerms <- acc_GO[[1]]$Term
		acc_GOPvalue <- acc_GO[[1]]$Pvalue
		acc_GOCount <- acc_GO[[1]]$Count
	} else if (catgo == "MF") {
		acc_GOCodes <- acc_GO[[2]]$GOMFID
		acc_GOTerms <- acc_GO[[2]]$Term
		acc_GOPvalue <- acc_GO[[2]]$Pvalue
		acc_GOCount <- acc_GO[[2]]$Count
	} else {
		acc_GOCodes <- acc_GO[[3]]$GOCCID
		acc_GOTerms <- acc_GO[[3]]$Term
		acc_GOPvalue <- acc_GO[[3]]$Pvalue
		acc_GOCount <- acc_GO[[3]]$Count
	}
	
	acc_Genes <- lookUp(unique(acc_GOCodes), annot, "GO2ALLPROBES")
	acc_GenesList <- data.frame()
	for(i in 1:length(acc_Genes)) {
		if(i==1) {
			int <- unique(intersect(acc_Genes[[i]], accSel))
			name_int <- lookUp(int, annot, "SYMBOL");
			acc_GenesList <- cbind(acc_GOCodes[i], acc_GOTerms[i], acc_GOPvalue[i], int, name_int)
		} else {
			int <- unique(intersect(acc_Genes[[i]], accSel))
			name_int <- lookUp(int, annot, "SYMBOL");
			start <- dim(acc_Genes)[1]
			acc_GenesList <- rbind(acc_GenesList, cbind(acc_GOCodes[i], acc_GOTerms[i], acc_GOPvalue[i], int, name_int))
		}
	}
	colnames(acc_GenesList) <- c("GOID", "GOTerm", "p-value", "Probeset", "GeneName")
	rownames(acc_GenesList) <- acc_GenesList[,4]
	
	if (formatMat) {
		
		acc_Mat<-matrix(data=0,nrow=length(unique(acc_GenesList[,4])),ncol=length(unique(acc_GenesList[,1])))
		rownames(acc_Mat)<-unique(paste(acc_GenesList[,4],acc_GenesList[,5],sep="-"))
		colnames(acc_Mat)<-unique(paste(acc_GenesList[,1],acc_GenesList[,2],sep="-"))
		for(i in 1:length(acc_GenesList[,4])){
			acc_Mat[paste(acc_GenesList[i,4],acc_GenesList[i,5],sep="-"),paste(acc_GenesList[grep(acc_GenesList[i,4],acc_GenesList[,4]),1],acc_GenesList[grep(acc_GenesList[i,4],acc_GenesList[,4]),2],sep="-")] <- -log10(as.numeric(acc_GenesList[grep(acc_GenesList[i,4],acc_GenesList[,4]),3]))
		}
		acc_Mat
		
	} else {	
		acc_GenesList
	}
}
# Ha dejado de funcionar summary() así que he modificado GOVik para que siga funcionando. Hay un error al sacar los DEG anotados a las categorías enriquecidas
# gsCollection = gsGOHumanAll or gsGOMouseAll
# GO_eli <- function(sel, ref, gsSet=gsGOHumanAll, annot=annotEnsembl) {
# 	paramsBP <- GSEAGOHyperGParams(name = "GOBPTest",
# 	    geneSetCollection = gsSet,
# 	    geneIds = intersect(sel, ref),
# 	    universeGeneIds = ref,
# 	    ontology="BP",
# 	    pvalueCutoff = 0.01,
# 	    conditional = TRUE,
# 	    testDirection = "over")
# 	hgOverBP <- hyperGTest(paramsBP)
# 	#resultGOBP <- summary(hgOverBP, pvalue=0.01, categorySize=2, htmlLinks=FALSE);
# 	x <- sigCategories(hgOverBP)
# 	GO_BP <- data.frame(GOID=rep(paste(x)[1],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), P.Value=rep(pvalues(hgOverBP)[x[1]],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), OR=rep(oddsRatios(hgOverBP)[x[1]],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), ExpCount=rep(expectedCounts(hgOverBP)[x[1]],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), Count=paste(unlist(geneIdsByCategory(hgOverBP)[x[1]])), Size= rep(length(unlist(geneIdUniverse(hgOverBP)[x[1]])),length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))))
# 	for(i in 2:length(x))
# 	{
# 		GO_BP <- rbind(GO_BP, data.frame(GOID=rep(paste(x)[i],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), P.Value=rep(pvalues(hgOverBP)[x[i]],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), OR=rep(oddsRatios(hgOverBP)[x[i]],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), ExpCount=rep(expectedCounts(hgOverBP)[x[i]],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), Count=paste(unlist(geneIdsByCategory(hgOverBP)[x[i]])), Size= rep(length(unlist(geneIdUniverse(hgOverBP)[x[i]])),length(unlist(geneIdsByCategory(hgOverBP)[x[i]])))))
# 	}
# 	paramsMF <- GSEAGOHyperGParams(name = "GOMFTest",
#     geneSetCollection = gsSet,
#     geneIds = intersect(sel, ref),
#     universeGeneIds = ref,
#     ontology="MF",
#     pvalueCutoff = 0.01,
#     conditional = TRUE,
#     testDirection = "over")
# 	hgOverMF <- hyperGTest(paramsMF)
# 	x <- sigCategories(hgOverMF)
# 	GO_MF <- data.frame(GOID=rep(paste(x)[1],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), P.Value=rep(pvalues(hgOverMF)[x[1]],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), OR=rep(oddsRatios(hgOverMF)[x[1]],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), ExpCount=rep(expectedCounts(hgOverMF)[x[1]],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), Count=paste(unlist(geneIdsByCategory(hgOverMF)[x[1]])), Size= rep(length(unlist(geneIdUniverse(hgOverMF)[x[1]])),length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))))
# 	for(i in 2:length(x))
# 	{
# 		GO_MF <- rbind(GO_MF, data.frame(GOID=rep(paste(x)[i],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), P.Value=rep(pvalues(hgOverMF)[x[i]],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), OR=rep(oddsRatios(hgOverMF)[x[i]],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), ExpCount=rep(expectedCounts(hgOverMF)[x[i]],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), Count=paste(unlist(geneIdsByCategory(hgOverMF)[x[i]])), Size= rep(length(unlist(geneIdUniverse(hgOverMF)[x[i]])),length(unlist(geneIdsByCategory(hgOverMF)[x[i]])))))
# 	}
# 	paramsCC <- GSEAGOHyperGParams(name = "GOCCTest",
# 	    geneSetCollection = gsSet,
# 	    geneIds = intersect(sel, ref),
# 	    universeGeneIds = ref,
# 	    ontology="CC",
# 	    pvalueCutoff = 0.01,
# 	    conditional = TRUE,
# 	    testDirection = "over")
# 	hgOverCC <- hyperGTest(paramsCC)
# 	x <- sigCategories(hgOverCC)
# 	GO_CC <- data.frame(GOID=rep(paste(x)[1],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), P.Value=rep(pvalues(hgOverCC)[x[1]],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), OR=rep(oddsRatios(hgOverCC)[x[1]],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), ExpCount=rep(expectedCounts(hgOverCC)[x[1]],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), Count=paste(unlist(geneIdsByCategory(hgOverCC)[x[1]])), Size= rep(length(unlist(geneIdUniverse(hgOverCC)[x[1]])),length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))))
# 	for(i in 2:length(x))
# 	{
# 		GO_CC <- rbind(GO_CC, data.frame(GOID=rep(paste(x)[i],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), P.Value=rep(pvalues(hgOverCC)[x[i]],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), OR=rep(oddsRatios(hgOverCC)[x[i]],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), ExpCount=rep(expectedCounts(hgOverCC)[x[i]],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), Count=paste(unlist(geneIdsByCategory(hgOverCC)[x[i]])), Size= rep(length(unlist(geneIdUniverse(hgOverCC)[x[i]])),length(unlist(geneIdsByCategory(hgOverCC)[x[i]])))))
# 	}
# 	GO_res <- rbind(cbind(Ontology="BP",GO_BP[GO_BP$Size>=2,]), cbind(Ontology="MF",GO_MF[GO_MF$Size>=2,]), cbind(Ontology="CC",GO_CC[GO_CC$Size>=2,]))
# 	godb_doterm <- Term(GOTERM)
# 	GO_res_annot <- merge(cbind(GO_res[,1:2],GOTerm=godb_doterm[paste(GO_res[,2])], GO_res[3], GeneID=GO_res[,6]), annot, by.x=5, by.y=1)
# 	return(GO_res_annot[,c(2:5,1,6:ncol(GO_res_annot))])
# }

# 
# sampleGOMat <- function(goList, goNames) {
# 	goCats <- c()
# 	for (i in 1:length(goList)) {
# 		goCats <- c(goCats, colnames(goList[[i]]))	
# 	}
# 	
# 	goCats <- unique(goCats)
# 	
# 	sampleGO <-  matrix(data = 0, ncol = length(goList), nrow = length(goCats))		
# 	colnames(sampleGO) <- goNames
# 	rownames(sampleGO) <- goCats
# 		
# 	for (i in 1:length(goList)) {				
# 		sampleGO[colnames(goList[[i]]), i] <- apply((goList[[i]]), 2, max)
# 	}
# 	
# 	sampleGO
# 	
# }
#

# Análisis funcional con geneSet de Ensembl (cargar rda de anotaciones que está en la carpeta de funciones)
library(GSEABase)
library(GOstats)
"goGSEnrich" <- function(selSet, refSet, gsSet, geneAnnot, hgCutoff, catSize) {

    paramsBP <- GSEAGOHyperGParams(name = "GOBPTest",
    geneSetCollection = gsSet,
    geneIds = selSet,
    universeGeneIds = refSet,
    ontology="BP",
    pvalueCutoff = hgCutoff,
    conditional = TRUE,
    testDirection = "over")

    print("Analizing gene list. GO Biological Process enrichment ...")

    hgOverBP <- hyperGTest(paramsBP)
    resultGOBP <- summary(hgOverBP, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);

    paramsMF <- GSEAGOHyperGParams(name = "GOMFTest",
    geneSetCollection = gsSet,    
    geneIds = selSet,
    universeGeneIds = refSet,
    ontology="MF",
    pvalueCutoff = hgCutoff,
    conditional = TRUE,
    testDirection = "over")

    print("Analizing gene list. GO Molecular Function enrichment ...")

    hgOverMF <- hyperGTest(paramsMF)
    resultGOMF <- summary(hgOverMF, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);

    paramsCC <- GSEAGOHyperGParams(name = "GOCCTest",
    geneSetCollection = gsSet,
    geneIds = selSet,
    universeGeneIds = refSet,
    ontology="CC",
    pvalueCutoff = hgCutoff,
    conditional = TRUE,
    testDirection = "over")

    print("Analizing gene list. GO Cellular Component enrichment ...")

    hgOverCC <- hyperGTest(paramsCC)
    resultGOCC <- summary(hgOverCC, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);

  if (nrow(resultGOBP) > 0) {
    colnames(resultGOBP)[1] <- "GOID"
    resultGOBP <- cbind(resultGOBP, "BP")
    colnames(resultGOBP)[8] <- "Ontology"
  }
  if (nrow(resultGOMF) > 0) {
    colnames(resultGOMF)[1] <- "GOID"
    resultGOMF <- cbind(resultGOMF, "MF")
    colnames(resultGOMF)[8] <- "Ontology"
  }
  if (nrow(resultGOCC) > 0) {
    colnames(resultGOCC)[1] <- "GOID"
    resultGOCC <- cbind(resultGOCC, "CC")
    colnames(resultGOCC)[8] <- "Ontology"
  }

  resGO <- rbind(resultGOBP, resultGOMF, resultGOCC)

  if (nrow(resGO) > 0) {
    resGOGenes <- data.frame()
    for(i in 1:nrow(resGO)) {
        genesGOtmp <- unique(intersect(unlist(geneIds(gsSet[resGO$GOID[i]])), selSet))

    if(i==1)
      resGOGenes <- cbind(paste(resGO$Ontology[i]), resGO$GOID[i], resGO$Term[i], resGO$Pvalue[i], genesGOtmp)
    else
      resGOGenes <- rbind(resGOGenes, cbind(paste(resGO$Ontology[i]), resGO$GOID[i], resGO$Term[i], resGO$Pvalue[i], genesGOtmp))
  }
    resGOGenes <- merge(resGOGenes, geneAnnot, by.x = 5, by.y = 1)
    resGOGenes <- resGOGenes[,c(2:5, 1, 6, 7)]
    colnames(resGOGenes) <- c("Ontology", "GOID", "GOTerm", "pvalue", "Gene", "Name", "Description")
  } else {
    resGOGenes = NULL
  }
  resGOGenes
}
 
# #Análisis funcional con geneSet de Ensembl (cargar rda de anotaciones que está en la carpeta de funciones)
# library(GOstats)                                      
# 
# goGSEnrich <- function(selSet, refSet, ontology, gsSet, geneAnnot, hgCutoff, catSize) {
#   if (ontology<2) {
#     paramsBP <- GSEAGOHyperGParams(name = "GOBPTest",
#     		geneSetCollection = gsSet,
#     		geneIds = selSet,
#     		universeGeneIds = refSet,
#     		ontology="BP",
#     		pvalueCutoff = hgCutoff,
#     		conditional = TRUE,
#     		testDirection = "over")
#   
#     print("Analizing gene list. GO Biological Process enrichment ...")
#       
#     hgOverBP <- hyperGTest(paramsBP)
#     resultGOBP <- summary(hgOverBP, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE)
#     colnames(resultGOBP)[1] <- "GOID"
#     resultGOBP <- cbind(resultGOBP, "BP") 
#     colnames(resultGOBP)[8] <- "Ontology"
#     resGO <- resultGOBP
#   
#   } 
#   
#   if (ontology==0 | ontology==2) {
#     paramsMF <- GSEAGOHyperGParams(name = "GOMFTest",
#     		geneSetCollection = gsSet,
#     		geneIds = selSet,
#     		universeGeneIds = refSet,
#     		ontology="MF",
#     		pvalueCutoff = hgCutoff,
#     		conditional = TRUE,
#     		testDirection = "over")
#     
#     print("Analizing gene list. GO Molecular Function enrichment ...")
#     
#     hgOverMF <- hyperGTest(paramsMF)
#     resultGOMF <- summary(hgOverMF, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE);
#     colnames(resultGOMF)[1] <- "GOID"
#     resultGOMF <- cbind(resultGOMF, "MF") 
#     colnames(resultGOMF)[8] <- "Ontology"
#     resGO <- resultGOMF         
# 
#   }
#   
#   if (ontology==0 | ontology==3) {  
#     paramsCC <- GSEAGOHyperGParams(name = "GOCCTest",
#     		geneSetCollection = gsSet,
#     		geneIds = selSet,
#     		universeGeneIds = refSet,
#     		ontology="CC",
#     		pvalueCutoff = hgCutoff,
#     		conditional = TRUE,
#     		testDirection = "over")
#     
#     print("Analizing gene list. GO Cellular Component enrichment ...")
#     
#     hgOverCC <- hyperGTest(paramsCC)
#     resultGOCC <- summary(hgOverCC, pvalue=hgCutoff, categorySize=catSize, htmlLinks=FALSE)
#     colnames(resultGOCC)[1] <- "GOID"
#     resultGOCC <- cbind(resultGOCC, "CC") 
#     colnames(resultGOCC)[8] <- "Ontology"
#     resGO <- resultGOCC
# 
#   }
# 
#   if(ontology==0){
#   	resGO <- rbind(resultGOBP, resultGOMF, resultGOCC)
#   }
#   
#   resGOGenes <- data.frame()
#   for(i in 1:nrow(resGO)) {                        
#       genesGOtmp <- unique(intersect(unlist(geneIds(gsGOHumanAll[resGO$GOID[i]])), selSet))
#     
#   		if(i==1)
#     			resGOGenes <- cbind(paste(resGO$Ontology[i]), resGO$GOID[i], resGO$Term[i], resGO$Pvalue[i], genesGOtmp)
#   		else
#     			resGOGenes <- rbind(resGOGenes, cbind(paste(resGO$Ontology[i]), resGO$GOID[i], resGO$Term[i], resGO$Pvalue[i], genesGOtmp))
# 	} 	
#   resGOGenes <- merge(resGOGenes, geneAnnot, by.x = 5, by.y = 1)
#   resGOGenes <- resGOGenes[,c(2:5, 1, 6, 7)]	
#   colnames(resGOGenes) <- c("Ontology", "GOID", "GOTerm", "pvalue", "Gene", "Name", "Description")
#   resGOGenes
# }

# Gráfico de barras para GO
"multiplot" <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


##########################################
#
# CLUSTERING ROBUSTO
#
##########################################

robClust <- function(data, niter = 1000, sampleIter = 0.8, nclusters = 3, verbose = FALSE) {

	clusters = matrix(0, ncol(data), niter)
	numSamples <- seq(1, ncol(data), 1)
	cdf <- matrix(0, 1, ncol(data)*(ncol(data)-1)/2)
	for (i in 1:niter) {
		indSample <- sample(numSamples, round(sampleIter*ncol(data)))
		sampleMat <- data[,indSample]
		hc <- hclust(dist(t(sampleMat)), method = "complete", members=NULL)
		clusters[indSample,i] <- cutree(hc, k=nclusters)
		if (verbose) print(i);
	}

	consensusClust <- matrix(0, ncol(data), ncol(data))
	colnames(consensusClust) <- colnames(data)
	rownames(consensusClust) <- colnames(data)
	k = 1;
	for (i in 1:(ncol(data)-1)) {
		for (j in (i+1):ncol(data)) {
			consensusClust[j,i] <- length(which((clusters[i,]-clusters[j,])==0))
			consensusClust[j,i] <- consensusClust[j,i] - length(which((clusters[j,which(clusters[j,] == 0)] - clusters[i,which(clusters[j,] == 0)])==0))
			consensusClust[j,i] <- consensusClust[j,i]/(length((clusters[j,which(clusters[j,] != 0)]*clusters[i,which(clusters[j,] != 0)]) != 0))
			cdf[1,k] <- consensusClust[j,i]
			k = k+1
			consensusClust[i,j] <- consensusClust[j,i]
		}
	}

	for (i in 1:ncol(data)) {
		consensusClust[i,i] <- 1
	}
	
	empCDF <- ecdf(cdf)
	aucdf <- sum(empCDF(sort(cdf))[-1]*diff(sort(cdf)))
		 
	rclust <- list(consensusClust = consensusClust, cdf = aucdf)
	rclust
}

##############################################################################
#### R script to
#### 	Produce Venn Diagrams with 1 to 5 groups
####		an extension on the code from the limma package
####		
#### Written By: Matt Settles
####				Postdoctoral Research Associate
####				Washington State University
####
##############################################################################
####
#### Change Log: 
####	Feb 8, 2008: 
####		formalized code
####    Dec 23, 2008:
####	    added mixed type to vennCounts
##############################################################################
####
####	Usage:
####	source("http://bioinfo-mite.crb.wsu.edu/Rcode/Venn.R")
####	can change colors now on 4 way plot
##############################################################################

####################################
## Function ellipse
## Add an elipse to the current plot
"ellipse" <- 
function (center, radius, rotate, 
    segments = 360, add = FALSE, xlab = "", ylab = "", las = par("las"), 
    col = palette()[2], lwd = 2, lty = 1, ...) 
{
	# x' = x cosø + y sinø
	# y' = y cosø - x sinø
    if (!(is.vector(center) && 2 == length(center))) 
        stop("center must be a vector of length 2")
    if (!(is.vector(radius) && 2 == length(radius))) 
        stop("radius must be a vector of length 2")
	
	angles <- (0:segments) * 2 * pi/segments  
	rotate <- rotate*pi/180
	ellipse <- cbind(radius[1] * cos(angles), 
					 radius[2] * sin(angles))
	if(rotate != 0)
		ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate),
						  ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
	ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
    if (add) 
        lines(ellipse, col = col, lwd = lwd, lty = lty, ...)
    else    plot(ellipse, type = "l", xlim = c(-4, 4), ylim = c(-4, 4),
			xlab = "", ylab = "", axes = FALSE, col = col, lwd = lwd,
			lty = lty, ...)
}

###################################
## Function vennCounts
## Produce venn table object
"vennCounts" <- 
function (x, include = "both") 
{
    x <- as.matrix(x)
    include <- match.arg(include, c("both", "up", "down","mixed"))
    x <- sign(switch(include, both = abs(x), up = x > 0, down = x < 0, mixed = x))
    nprobes <- nrow(x)
    ncontrasts <- ncol(x)
    names <- colnames(x)
    if (is.null(names)) 
        names <- paste("Group", 1:ncontrasts)
    noutcomes <- 2^ncontrasts
    if ( include == "mixed" ) noutcomes <- 3^ncontrasts
    outcomes <- matrix(0, noutcomes, ncontrasts)
    colnames(outcomes) <- names
    for (j in 1:ncontrasts) {
    	if( include == "mixed"){
    		outcomes[, j] <- rep(-1:1, times = 3^(j - 1), each = 3^(ncontrasts - j))
    	} else {
    	    outcomes[, j] <- rep(0:1, times = 2^(j - 1), each = 2^(ncontrasts - j))
    	}	
    }
    xlist <- list()
    for (i in 1:ncontrasts) {
    	if( include == "mixed"){
	    	xlist[[i]] <- factor(x[, ncontrasts - i + 1], levels = c(-1,0, 1))
		} else {
	       	xlist[[i]] <- factor(x[, ncontrasts - i + 1], levels = c(0, 1))
	    }    
    }
    counts <- as.vector(table(xlist))
    structure(cbind(outcomes, Counts = counts), class = "vennCounts")
}

"vennDiagram" <-
function (object, include = "both", names, mar = rep(1, 4), cex = 0.7, 
    lwd = 2, circle.col, counts.col, show.include, ...) 
{
    if (!is(object, "vennCounts")) {
        if (length(include) > 2) 
            stop("Cannot plot Venn diagram for more than 2 sets of counts")
        if (length(include) == 2) 
            object.2 <- vennCounts(object, include = include[2])
        object <- vennCounts(object, include = include[1])
    }
    else if (length(include == 2)) 
        include <- include[1]
    nsets <- ncol(object) - 1
    if (nsets > 4) 
        stop("Can't plot Venn diagram for more than 4 sets")
    if (missing(names)) 
        names <- colnames(object)[1:nsets]
    counts <- object[, "Counts"]
    if (length(include) == 2) 
    
    if (length(include) == 2 && nsets < 4) 
        counts.2 <- object.2[, "Counts"]
    #Setup colors
	 if (missing(circle.col) & nsets == 4)
        circle.col <- c("red","blue","orange","green")
    else 
        circle.col <- par("col")
    if (length(circle.col) < nsets) 
        circle.col <- rep(circle.col, length.out = nsets)
    if (missing(counts.col) & nsets == 4) 
        counts.col <- c("red","blue","orange","green")
	 else
	 	  counts.col <- par("col") 
    if (length(counts.col) < length(include) & nsets < 4) 
        counts.col <- rep(counts.col, length.out = length(include))
    else if (length(counts.col) < nsets) 
        counts.col <- rep(counts.col, length.out = nsets)
    
    if (missing(show.include)) 
        show.include <- as.logical(length(include) - 1)
    xcentres <- list(0, c(-1, 1), c(-1, 1, 0), c(-0.2,0.2,-1.05,1.05))[[nsets]]
    ycentres <- list(0, c(0, 0), c(1/sqrt(3), 1/sqrt(3), -2/sqrt(3)),c(.20,.20,-0.35,-0.35))[[nsets]]
    centers <- cbind(xcentres,ycentres)
	r1 <- c(1.5, 1.5, 1.5, 1.5)[nsets]
	r2 <- c(1.5, 1.5, 1.5, 2.7)[nsets]
	radius <- c(r1,r2)
	rotate <- list(0, c(0,0), c(0,0,0), c(-45,45,-45,45))[[nsets]]
	
    xtext <- list(-1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0),c(-3.2,3.2,-3.2,3.2))[[nsets]]
    ytext <- list(1.8, c(1.8, 1.8), c(2.4, 2.4, -3),c(3.2,3.2,-3.2,-3.2))[[nsets]]
    old.par <- par(mar = mar)
    on.exit(par(old.par))
    plot(x = 0, y = 0, type = "n", xlim = c(-4.0, 4.0), ylim = c(-4.0, 
        4.0), xlab = "", ylab = "", axes = FALSE, ...)
    for (circle in 1:nsets) {
		ellipse(centers[circle,],radius,rotate[circle],add=TRUE,
				, lwd = lwd, col = circle.col[circle])
        text(xtext[circle], ytext[circle], names[circle], 
			#pos=ifelse(circle != as.integer(circle/2)*2,4,2),
			offset = 0.5, col = circle.col[circle], cex = cex)
    }
#	rect(-3.5, -3.5, 3.5, 3.5)
 
	switch(nsets, {
        rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
			text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(0, 0, counts[2], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        rect(-3, -2.5, 3, 2.5)
        printing <- function(counts, cex, adj, col, leg) {
            col <- col[1]
			text(2.3, -2.1, counts[1], cex = cex, col = col, 
                adj = adj)
            text(1.5, 0.1, counts[2], cex = cex, col = col, adj = adj)
            text(-1.5, 0.1, counts[3], cex = cex, col = col, 
                adj = adj)
            text(0, 0.1, counts[4], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.3, -2.1, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        rect(-3, -3.5, 3, 3.3)
        printing <- function(counts, cex, adj, col, leg) {
			col <- col[1]
            text(2.5, -3, counts[1], cex = cex, col = col, adj = adj)
            text(0, -1.7, counts[2], cex = cex, col = col, adj = adj)
            text(1.5, 1, counts[3], cex = cex, col = col, adj = adj)
            text(0.75, -0.35, counts[4], cex = cex, col = col, 
                adj = adj)
            text(-1.5, 1, counts[5], cex = cex, col = col, adj = adj)
            text(-0.75, -0.35, counts[6], cex = cex, col = col, 
                adj = adj)
            text(0, 0.9, counts[7], cex = cex, col = col, adj = adj)
            text(0, 0, counts[8], cex = cex, col = col, adj = adj)
            if (show.include) 
                text(-2.5, -3, leg, cex = cex, col = col, adj = adj)
        }
    }, {
        rect(-3.5, -3.5, 3.5, 3.5)
        printing <- function(counts, cex, adj, col, leg) {
            text(0, -3, 		counts[1], cex = cex, col = "black", adj = adj)
            text(2.5, 0, 		counts[2], cex = cex, col = col[4], adj = adj)
			 lines(c(2.25,2.75),c(-0.2,-0.2),col=col[4])
			 
            text(-2.5, 0, 		counts[3], cex = cex, col = col[3], adj = adj)
             lines(c(-2.75,-2.25),c(-0.2,-0.2),col=col[3])
			 
			text(0, -2.0, 		counts[4], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(-2.2,-2.2),col=col[3])
			 lines(c(-0.25,0.25),c(-2.25,-2.25),col=col[4])
			
			text(1.3, 2.1, 		counts[5], cex = cex, col = col[2], adj = adj)
             lines(c(1.05,1.55),c(1.9,1.9),col=col[2])
            
			text(1.7, 1.2, 		counts[6], cex = cex, col = "black", adj = adj)
             lines(c(1.45,1.95),c(1.0,1.0),col=col[2])
             lines(c(1.45,1.95),c(0.95,0.95),col=col[4])
            
			text(-1.6, -1.1, 	counts[7], cex = cex, col = "black", adj = adj)
             lines(c(-1.85,-1.35),c(-1.3,-1.3),col=col[2])
			 lines(c(-1.85,-1.35),c(-1.35,-1.35),col=col[3])
			 
			text(-0.8, -1.55, 	counts[8], cex = cex, col = "black", adj = adj)
			 lines(c(-0.55,-1.05),c(-1.75,-1.75),col=col[2])
			 lines(c(-0.55,-1.05),c(-1.8,-1.8),col=col[3])
			 lines(c(-0.55,-1.05),c(-1.85,-1.85),col=col[4])

			text(-1.3, 2.1, 	counts[9], cex = cex, col = col[1], adj = adj)
			 lines(c(-1.55,-1.05),c(1.9,1.9),col=col[1])
            
			text(1.6, -1.1, 	counts[10], cex = cex, col = "black", adj = adj)
             lines(c(1.85,1.35),c(-1.3,-1.3),col=col[1])
			 lines(c(1.85,1.35),c(-1.35,-1.35),col=col[4])

			text(-1.7, 1.2, 	counts[11], cex = cex, col = "black", adj = adj)
			 lines(c(-1.45,-1.95),c(1.0,1.0),col=col[1])
             lines(c(-1.45,-1.95),c(0.95,0.95),col=col[3])
            
			text(0.8, -1.55, 	counts[12], cex = cex, col = "black", adj = adj)
			 lines(c(0.55,1.05),c(-1.75,-1.75),col=col[1])
			 lines(c(0.55,1.05),c(-1.8,-1.8),col=col[3])
			 lines(c(0.55,1.05),c(-1.85,-1.85),col=col[4])
			 
			text(0, 1.6, 		counts[13], cex = cex, col = "black", adj = adj)
             lines(c(-0.25,0.25),c(1.4,1.4),col=col[1])
			 lines(c(-0.25,0.25),c(1.35,1.35),col=col[2])
			 
			text(0.9, 0.5, 		counts[14], cex = cex, col = "black", adj = adj)
             lines(c(1.15,0.65),c(0.3,0.3),col=col[1])
			 lines(c(1.15,0.65),c(0.25,0.25),col=col[2])
			 lines(c(1.15,0.65),c(0.2,0.2),col=col[4])
			 
			text(-0.9, 0.5, 	counts[15], cex = cex, col = "black", adj = adj)
             lines(c(-1.15,-0.65),c(0.3,0.3),col=col[1])
			 lines(c(-1.15,-0.65),c(0.25,0.25),col=col[2])
			 lines(c(-1.15,-0.65),c(0.2,0.2),col=col[3])
			 
			text(0, -0.5, 		counts[16], cex = cex, col = "black", adj = adj)
			 lines(c(-0.25,0.25),c(-0.7,-0.7),col=col[1])
			 lines(c(-0.25,0.25),c(-0.75,-0.75),col=col[2])
			 lines(c(-0.25,0.25),c(-0.8,-0.8),col=col[3])
			 lines(c(-0.25,0.25),c(-0.85,-0.85),col=col[4])			  
      } 
	} )
    adj <- c(0.5, 0.5)
    if (length(include) == 2 & nsets < 4) 
        adj <- c(0.5, 0)
    printing(counts, cex, adj, counts.col, include[1])
    if (length(include) == 2 & nsets < 4) 
        printing(counts.2, cex, c(0.5, 1), counts.col[2], include[2])
    invisible()
}

#######################################################
#
#      Función para hacer gráficos de t-test...
#
#######################################################

tplot <- function(x, ...) UseMethod("tplot")
tplot.default <- function(x, ..., type="d", dist=NULL, jit=0.05, names, xlim=NULL, ylim=NULL, main=NULL, sub=NULL, xlab=NULL, ylab=NULL, col=par("col"), pch=par("pch"), group.col=FALSE, group.pch=FALSE, median.line=FALSE, mean.line=FALSE, median.pars=list(col=par("col")), mean.pars=median.pars, boxplot.pars=NULL, show.n=FALSE, my.gray=gray(.75), ann=par("ann"), axes=TRUE, frame.plot=axes, add=FALSE, at=NULL, horizontal=FALSE, panel.first=NULL, panel.last=NULL) {
    localAxis <- function(..., bg, cex, lty, lwd) axis(...)
    localBox <- function(..., bg, cex, lty, lwd) box(...)
    localWindow <- function(..., bg, cex, lty, lwd) plot.window(...)
    localTitle <- function(..., bg, cex, lty, lwd) title(...)

    args <- list(x, ...)
    namedargs <- if (!is.null(attributes(args)$names))
        attributes(args)$names != ""
    else logical(length(args))
    groups <- if (is.list(x))
        x
    else args[!namedargs]
    pars <- args[namedargs]
    if ((n <- length(groups)) == 0)
        stop("invalid first argument")
    if (length(class(groups)))
        groups <- unclass(groups)
    if (!missing(names))
        attr(groups, "names") <- names
    else {
        if (is.null(attr(groups, "names")))
            attr(groups, "names") <- 1:n
        names <- attr(groups, "names")
    }

    ng <- length(groups) # number of groups
    l <- sapply(groups, length) # size of each group
    g <- rep(1:ng, l) # groups as.numeric
    nv <- sum(l) # total count

    # set y scale
    ylim <- if (!is.null(ylim))
        ylim
    else {
        r <- range(groups, na.rm=TRUE)
        pm <- diff(r) / 20
        r + pm * c(-1,1)
    }
    # set x scale
    if (is.null(xlim)) xlim <- c(0.5, ng+0.5)

    at <- if (is.null(at)) 1:ng else at
    if (length(at) != ng)
        stop("'at' must have same length as the number of groups")

    xlab <- if (is.null(xlab)) "" else xlab
    ylab <- if (is.null(ylab)) "" else ylab
    main <- if (is.null(main)) "" else main
    sub <- if (is.null(sub)) "" else sub

    type <- match.arg(type, choices=c("d","db","bd","b"), several.ok=TRUE)
    # type of plot for each group
    if ((length(type) > 1) && (length(type) != ng))
        warning("length of 'type' does not match the number of groups")
    type <- rep(type, length.out=ng)

    # Use colors by group
    if (group.col) {
        if (length(col) != ng)
            warning("length of 'col' does not match the number of groups")
        g.col <- rep(col, length.out=ng)
        col <- rep(g.col, l)
    # Use colors by individual or global
    } else {
        if((length(col) > 1) && (length(col) != nv))
            warning("length of 'col' does not match the number of data points")
        col <- rep(col, length.out=nv)
        g.col <- rep(1, length.out=ng)
    }

    # Use plot characters by group
    if (group.pch) {
        if (length(pch) != ng)
            warning("length of 'pch' does not match the number of groups")
        pch <- rep(rep(pch, length.out=ng), l)
    # Use plot characters by individual or global
    } else {
        if((length(pch) > 1) && (length(pch) != nv))
            warning("length of 'pch' does not match the number of data points")
        pch <- rep(pch, length.out=nv)
    }

    # split colors and plot characters into groups
    col <- split(col, g)
    pch <- split(pch, g)
    # remove any NAs from the data and options
    nonas <- lapply(groups, function(x) !is.na(x))
    groups <- mapply("[", groups, nonas, SIMPLIFY=FALSE)
    l <- sapply(groups, length)
    col <- mapply("[", col, nonas, SIMPLIFY=FALSE)
    pch <- mapply("[", pch, nonas, SIMPLIFY=FALSE)

    # whether or not to display a mean and median line for each group
    mean.line <- rep(mean.line, length.out=ng)
    median.line <- rep(median.line, length.out=ng)

    # set defaults for dist and jit
    if (is.null(dist) || is.na(dist)) dist <- diff(range(ylim)) / 100
    if (is.null(jit) || is.na(jit)) jit <- 0.025 * ng

    # 1 2 3 1 3 2 1 1 4 2
    # -------------------
    # 1 1 1 2 2 2 3 4 1 3
    how.many.so.far <- function(g) {
        out <- NULL
        u <- unique(g)
        for (i in 1:length(u)) out[which(g==u[i])] <- 1:sum(g==u[i])
        out
    }
    # turns the values in each group into their plotting points
    grouping <- function(v, dif) {
        vs <- sort(v)
        together <- c(FALSE, diff(vs) <= dif)
        g.id <- cumsum(!together)
        g.si <- rep(x<-as.vector(table(g.id)), x)
        vg <- cbind(vs, g.id, g.si)[rank(v),]
        if (length(v)==1) vg <- as.data.frame(t(vg))
        hmsf <- how.many.so.far(vg[,2])
        data.frame(vg, hmsf)
    }
    groups <- lapply(groups, grouping, dif=dist)

    # set up new plot unless adding to existing one
    if (!add) {
        plot.new()
        if (horizontal)
            do.call("localWindow", c(list(ylim, xlim), pars))
        else
            do.call("localWindow", c(list(xlim, ylim), pars))
    }
    panel.first

    # function to compute the jittering
    jit.f2 <- function(g.si, hm.sf) { hm.sf - (g.si + 1) / 2 }

    out <- list()

    Lme <- 0.2 * c(-1, 1)
    for (i in 1:ng) {
        to.plot <- groups[[i]]
        gs <- to.plot$g.si
        hms <- to.plot$hm
        x <- rep(at[i], nrow(to.plot)) + jit.f2(gs, hms) * jit
        y <- to.plot$vs

        if (type[i] == "bd") { # dots behind
            if (horizontal)
                do.call("points", c(list(x=y, y=x, pch=pch[[i]], col=my.gray), pars))
            else
                do.call("points", c(list(x=x, y=y, pch=pch[[i]], col=my.gray), pars))
        }
        if (type[i] %in% c("bd", "b")) { # boxplot in front
            outliers <- do.call("boxplot", c(list(x=y, at=at[i], add=TRUE, axes=FALSE, border=g.col[i], outline=FALSE, horizontal=horizontal), boxplot.pars))$out
            if (type[i] == "b") {
                toplot <- rowSums(outer(y, outliers, "==")) == 1
                if (horizontal)
                    do.call("points", c(list(x=y[toplot], y=x[toplot], pch=pch[[i]][toplot], col=col[[i]][toplot]), pars))
                else
                    do.call("points", c(list(x=x[toplot], y=y[toplot], pch=pch[[i]][toplot], col=col[[i]][toplot]), pars))
            }
        }
        if (type[i] == "db") # boxplot behind
            do.call("boxplot", c(list(x=y, at=at[i], add=TRUE, axes=FALSE, border=my.gray, outline=FALSE, horizontal=horizontal), boxplot.pars))
        if (type[i] %in% c("db", "d")) { # dots in front
            if (horizontal)
                do.call("points", c(list(x=y, y=x, pch=pch[[i]], col=col[[i]]), pars))
            else
                do.call("points", c(list(x=x, y=y, pch=pch[[i]], col=col[[i]]), pars))
        }
        if (mean.line[i]) { # mean line
            if (horizontal)
                do.call("lines", c(list(rep(mean(y), at[i]+Lme, 2)), mean.pars))
            else
                do.call("lines", c(list(at[i]+Lme, rep(mean(y), 2)), mean.pars))
        }
        if (median.line[i]) { # median line
            if (horizontal)
                do.call("lines", c(list(rep(median(y), at[i]+Lme, 2)), median.pars))
            else
                do.call("lines", c(list(at[i]+Lme, rep(median(y), 2)), median.pars))
        }

        out[[i]] <- data.frame(to.plot, col=col[[i]], pch=pch[[i]])
    }
    panel.last

    # add axes
    if (axes) {
        do.call("localAxis", c(list(side=1+horizontal, at=1:ng, labels=names, tcl=0), pars))
        do.call("localAxis", c(list(side=2-horizontal), pars))
    }
    # optional sample sizes
    if (show.n)
        do.call("localAxis", c(list(side=3+horizontal, at=1:ng, tick=F, labels=paste("n=", l, sep=""), tcl=0), pars, list(mgp=c(3,.5,1))))
    # add bounding box
    if (frame.plot)
        do.call("localBox", pars)
    # add titles
    if (ann) {
        if (horizontal)
            do.call("localTitle", c(list(main=main, sub=sub, xlab=ylab, ylab=xlab), pars))
        else
            do.call("localTitle", c(list(main=main, sub=sub, xlab=xlab, ylab=ylab), pars))
    }

    invisible(out)
}

tplot.formula <- function(formula, data=parent.frame(), ..., subset) {
    if (missing(formula) || (length(formula) != 3))
        stop("'formula' missing or incorrect")

    enquote <- function(x) { as.call(list(as.name("quote"), x)) }

    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)

    args <- lapply(m$..., eval, data, parent.frame())
    nmargs <- names(args)
    if ("main" %in% nmargs) args[["main"]] <- enquote(args[["main"]])
    if ("sub" %in% nmargs) args[["sub"]] <- enquote(args[["sub"]])
    if ("xlab" %in% nmargs) args[["xlab"]] <- enquote(args[["xlab"]])
    if ("ylab" %in% nmargs) args[["ylab"]] <- enquote(args[["ylab"]])

    m$... <- NULL
    m$na.action <- na.pass
    subset.expr <- m$subset
    m$subset <- NULL
    require(stats, quietly=TRUE)
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())
    response <- attr(attr(mf, "terms"), "response")

    ## special handling of col and pch
    n <- nrow(mf)
    # rep as necessary
    col <- if ("col" %in% names(args)) args$col else par("col")
    pch <- if ("pch" %in% names(args)) args$pch else par("pch")
    # pick out these options
    group.col <- if ("group.col" %in% names(args)) args$group.col else FALSE
    group.pch <- if ("group.pch" %in% names(args)) args$group.pch else FALSE
    # reorder if necessary
    if (!group.col) args$col <- unlist(split(rep(col, length.out=n), mf[-response]))
    if (!group.pch) args$pch <- unlist(split(rep(pch, length.out=n), mf[-response]))

    if (!missing(subset)) {
        s <- eval(subset.expr, data, parent.frame())
        dosub <- function(x) { if (length(x) == n) x[s] else x }
        args <- lapply(args, dosub)
        mf <- mf[s,]
    }
    do.call("tplot", c(list(split(mf[[response]], mf[-response])), args))
}
# EXAMPLE
#set.seed(100)
#age <- rnorm(80,rep(c(26,36),c(70,10)),4)
#sex <- sample(c('Female','Male'),80,T)
#group <- paste('Group ', sample(1:4,40,prob=c(2,5,4,1),replace=T), sep='')
#d <- data.frame(age,group)
#tplot(age~group,data=d,show.n=F,type=c('d','d','d','d'), dist=.2,jit=.05,  las=1,pch=19)

# plotting along chromosome (extra bar/s, gff optional, legend at bottom optional)
# thresh argument missing in plotAlongChrom (used for setting global ylim for heatmaps)
# also need ability to change the color scheme for panel and main heatmap
plotAlongChrom = function(segObj, y, probeAnno, gff,
                          isDirectHybe=FALSE, 
                          what = c("dots"), ## "heatmap"
                          chr, coord, highlight,  
                          colors, 
                          doLegend=FALSE,
                          featureExclude=c("chromosome", "nucleotide_match", "insertion"),
                          featureColorScheme=1, extras, 
                          rowNamesHeatmap, rowNamesExtras, ylab, ylabExtras, main,
                          colHeatmap=colorRamp(brewer.pal(9, "YlGnBu")),
                          colExtras=colorRamp(brewer.pal(9, "Reds")), 
                          sepPlots=FALSE, reOrder=TRUE,...) {

  ## set up the viewports of the plot layout.
  VP = c("title"=0.4, "expr+"=5, "gff+"=1, "coord"=1, "gff-"=1, "expr-"=5, "legend"=0.4)

  if(sepPlots) { # matrix
     if(!missing(y))
        n <- ncol(y)
     if(!missing(segObj)) { # S4
        if(is.environment(segObj)) {
          segmentationObjectName = paste(chr, "+", sep=".")
          if(segmentationObjectName %in% ls(segObj)) {
            s <- get(segmentationObjectName, segObj)
            n <- ncol(s@y)
          }
          else { # old style list
            dat = get(paste(chr, strand, "dat", sep="."), segObj)
            n <- ncol(dat$y)
          }
        }
     }
     if(reOrder)
       ordering = seq(n,1, by=-1)
     else
       ordering = seq(1:n)
     if(n<4) {
         exprw <- exprc <- NULL
         for(i in 1:n) {
           exprw <- c(exprw, paste("expr", i, "+", sep=""))
           exprc <- c(exprc, paste("expr", i, "-", sep=""))
           }
         VPnames <- c("title", exprw, "gff+", "coord", "gff-", exprc, "legend")
         VP = c(0.8, rep(5, n), 1, 1, 1, rep(5, n), 0.4)
         names(VP) <- VPnames
      } else{
          cat("More than 4 arrays, plotting averages.\n")
          sepPlots=FALSE
      }
  }

  if(!missing(extras)) {
     indgff <-  grep("gff\\+", names(VP))
     indlegend <-  grep("legend", names(VP)) 
     VP <- c(VP[1:(indgff-1)], "extras+"=1, VP[indgff:(indlegend-1)], "extras-"=1, VP[indlegend:length(VP)])
  }
  if(!doLegend)
     VP = VP[-which(names(VP)=="legend")]
  if(missing(gff))
     VP = VP[-which(names(VP)=="gff+" | names(VP)=="gff-")]
  defaultColors = c("+" = "#00441b", "-" = "#081d58", "duplicated" = "grey",
    "cp" = "#555555", "ci" = "#777777", "highlight" = "red", "threshold" = "grey")
  if(!missing(colors)) {
    mt = match(names(colors), names(defaultColors))
    if(any(is.na(mt)))
    stop(paste("Cannot use color specification for", names(colors)[is.na(mt)]))
    defaultColors[mt] = colors 
  }
  colors = defaultColors

  ## check that either y or segObj is present
  if(!missing(y)) {
    if(missing(probeAnno))
      stop("If 'y' is specified, 'probeAnno' must also be specified.")
    if(!missing(segObj))
      stop("If 'y' is specified, 'segObj' must not be specified.")
  } else {
    if(missing(segObj))
      stop("Please specify either 'y' or 'segObj'")
  }

  pushViewport(viewport(width=0.85, height=0.95)) ## plot margin
  pushViewport(viewport(layout=grid.layout(length(VP), 1, height=VP)))
  for(i in 1:2) {
    strand = c("+", "-")[i]

    ## extract and treat  the data
    threshold = as.numeric(NA)

    ## Three mutually exclusive cases:
    ## 1.) y and probeAnno
    ## 2.) segObj is an environment and contains objects of S4 class "segmentation"
    ##    whose names are obtained by paste(chr, strand, sep=".")
    ## 3.) segObj is an environment and contains lists
    ##    whose names are obtained by paste(chr, strand, "dat", sep=".")
    ##
    if(!missing(y)) {
      ## case 1.
      stopifnot(is.matrix(y))
      index = get(paste(chr, strand, "index", sep="."), envir=probeAnno)
      sta   = get(paste(chr, strand, "start", sep="."), envir=probeAnno)
      end   = get(paste(chr, strand, "end",   sep="."), envir=probeAnno)
      if(!missing(extras))
        dat = list(x = (sta+end)/2,
                   y   = y[index,, drop=FALSE],
                   flag = get(paste(chr, strand, "unique", sep="."), envir=probeAnno),  
                   extras = extras[index,, drop=FALSE]) # extras not currently supported 
							# for case 2 or 3.
      else   
        dat = list(x = (sta+end)/2,
                   y = y[index,, drop=FALSE],
                   flag = get(paste(chr, strand, "unique", sep="."), envir=probeAnno))
      stopifnot(is.numeric(dat$flag))
      lengthChr = end[length(end)]
      
    } else {
      if(!is.environment(segObj))
        stop("'segObj' must be an environment.")
      
      segmentationObjectName = paste(chr, strand, sep=".")
      if(segmentationObjectName %in% ls(segObj)) {
        ## case 2: S4 class
        s = get(segmentationObjectName, segObj)
        if(!inherits(s, "segmentation"))
          stop(sprintf("'%s' must be of class'segmentation'.", segmentationObjectName))
        if(is.na(s@nrSegments))
          stop(sprintf("Slot 'nrSegments' of '%s' must not be NA.", segmentationObjectName))
        bp = s@breakpoints[[s@nrSegments]]
        dat = list(x=s@x, y=s@y, flag=s@flag, estimate = bp[, "estimate"])
        if("upper" %in% colnames(bp)) dat$upper = bp[, "upper"]
        if("lower" %in% colnames(bp)) dat$lower = bp[, "lower"]
        lengthChr <- max(s@x, na.rm=TRUE)

      } else {
        ## case 3: list 'dat' and other stuff
        dat = get(paste(chr, strand, "dat", sep="."), segObj)
        stopifnot(all(c("start", "end", "unique", "ss") %in% names(dat)))
        dat$x = (dat$start + dat$end)/2
        dat$flag = dat$unique
        lengthChr = dat$end[length(dat$end)]
        
        if("segScore" %in% ls(segObj)) {
          sgs = get("segScore", segObj)
          sgs = sgs[ sgs$chr==chr & sgs$strand==strand, c("start", "end") ]
        } else {
	  stop("This option is deprecated")
	  ##nrSegments = ...
          ##seg = get(paste(chr, strand, "seg", sep="."), segObj)
          ##th  = c(1, seg$th[nrSegments+1, 1:(nrSegments+1)])
          ##sgs = list(start  = dat$start[dat$ss][th[-length(th)]],
          ##           end    = dat$end[dat$ss][th[-1]]-1)
        }
        dat$estimate = (sgs$start[-1] + sgs$end[-length(sgs$end)]) / 2

        if("theThreshold" %in% ls(segObj))
          threshold = get("theThreshold", segObj)
      } 
    }
    ## At this point, no matter what the input to the function was,
    ##    we have the list 'dat' with elements
    ## x: x coordinate
    ## y: y coordinate
    ## flag
    ## extras (optional, for plotting an extra panel, such as p-values for each segment)
    ## estimate (optional)
    ## lower (optional)
    ## upper (optional)
    ## and possibly also a non-NA value for 'threshold'.

    ## if no region is specified, plot the whole chromosome
    if(missing(coord))
      coord = c(1, lengthChr)
    
    ## plot the data
    vpr=which(names(VP)==sprintf("expr%s", strand))
    switch(what,
      "dots" = {
      if(sepPlots) {
       ylimdata = quantile(as.vector(dat[["y"]][dat[["x"]]>=coord[1] & dat[["x"]]<=coord[2],]), c(0, 1), na.rm=TRUE)
       ylim=ylimdata 
       if(missing(ylab))
           ylab=colnames(dat$y)
       if(length(ylab)==1)
           ylab=rep(ylab, n)
       for(j in seq(1:n)) {
         datj <- dat
         datj$y <- dat$y[,ordering[j]]
         if(missing(ylab))
           ylab=colnames(dat$y)
         ## plot the data
         vpr=which(names(VP)==sprintf(paste("expr",j,"%s",sep=""), strand))
  
         plotSegmentationDots(datj, xlim=coord, ylim=ylim, ylab=ylab[ordering[j]], 
                         chr=chr, strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                         vpr=vpr, colors=colors, sepPlots=sepPlots,...)
         } 
      } else
         plotSegmentationDots(dat, xlim=coord, ylab=ylab, 
                         chr=chr, strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                         vpr=vpr, colors=colors, sepPlots=sepPlots,...)
 
       if(!missing(extras) & !missing(y)) {
             vpr2=which(names(VP)==sprintf("extras%s", strand))
             dat$y = dat$extras[,,drop=FALSE]
             plotSegmentationDots(dat, xlim=coord, chr=chr, 
                     strand=ifelse(isDirectHybe, otherStrand(strand),strand),
                     vpr=vpr2, colors=colors, colHeatmap=colExtras, 
                     ylab=ylabExtras, rowNames=rowNamesExtras,...)
        }
      },

      ## FIXME: Matt's spaghetti code needs cleanup
           
      "heatmap" = {
        plotSegmentationHeatmap(dat, xlim=coord, 
                                rowNames=rowNamesHeatmap,
                                chr=chr, 
                                strand=ifelse(isDirectHybe, otherStrand(strand), strand),
                                vpr=vpr, colors=colors, ylab=ylab,
                                colHeatmap=colHeatmap,...)
        if(!missing(extras) & !missing(y)) {
             vpr2=which(names(VP)==sprintf("extras%s", strand))
             dat$y = dat$extras[,,drop=FALSE]
             plotSegmentationHeatmap(dat, xlim=coord, chr=chr, 
                     strand=ifelse(isDirectHybe, otherStrand(strand),strand),
                     vpr=vpr2, colors=colors, colHeatmap=colExtras, 
                     ylab=ylabExtras, rowNames=rowNamesExtras,...)
        }
      },
           stop(sprintf("Invalid value '%s' for argument 'what'", what))
    ) ## switch
  
    ## plot the features
    if(!missing(gff))
      plotFeatures(gff=gff, chr=chr, xlim=coord, strand=strand, 
                   featureExclude=featureExclude, featureColorScheme=featureColorScheme,
                   vpr=which(names(VP)==sprintf("gff%s", strand)),...) 
  }

  ## chromosomal coordinates
  pushViewport(dataViewport(xData=coord, yscale=c(-0.4,0.8), extension=0, 
                            layout.pos.col=1, layout.pos.row=which(names(VP)=="coord")))
  grid.lines(coord, c(0,0), default.units = "native")
  tck = alongChromTicks(coord)
  grid.text(label=formatC(tck, format="d"), x = tck, y = 0.2, 
            just = c("centre", "bottom"), gp = gpar(cex=.6), default.units = "native")
  grid.segments(x0 = tck, x1 = tck, y0 = -0.17, y1 = 0.17,  default.units = "native")

  
  if(!missing(highlight)){
    ## this part was modified to draw arrows for transcripts rather than bars
    mt = (match(highlight$strand, c("-", "+"))-1.5)*2
    co = highlight$coord
    if(is.na(mt) || !is.numeric(co))
      stop("Invalid parameter 'highlight'.")
    strand.num <- ifelse(highlight$strand=="-",-1,1)
    grid.segments(x0=co, x1=co+(500*strand.num), y0=c(0.4,0.4)*mt, y1=c(0.4,0.4)*mt, default.units = "native", arrow=arrow(), gp=gpar(col="violetred4", lwd=4))
  }
  popViewport()

  ## title
  pushViewport(viewport(layout.pos.col=1, layout.pos.row=which(names(VP)=="title")))
  grid.text(label=paste("Chr ", chr, sep=""), x=0.5, y=1, just="centre", gp=gpar(cex=1))
  if(!missing(main))
    grid.text(label=main, x=0.05, y=1, just="centre", gp=gpar(cex=1))
  popViewport()

  ## legend
  if(doLegend)
    plotAlongChromLegend(which(names(VP)=="legend"),
         featureColorScheme=featureColorScheme, featureExclude=featureExclude)
  
  popViewport(2)
}

## ------------------------------------------------------------
## plot Features
## ------------------------------------------------------------

plotFeatures = function(gff, chr, xlim, strand, vpr, featureColorScheme=1, featureExclude=c("chromosome", "nucleotide_match", "insertion"), featureNoLabel=c("uORF", "CDS"),...) {

  pushViewport(dataViewport(xData=xlim, yscale=c(-1.2,1.2),  extension=0, clip="on",
    layout.pos.col=1, layout.pos.row=vpr))

  stopifnot(all(gff[,"start"] <= gff[, "end"]))
  sel = which(gff[, "chr"] == chr &
              gff[, "strand"]  == strand &
              gff[, "start"] <= xlim[2] &
              gff[, "end"]   >= xlim[1])

  stopifnot(length(strand)==1, strand %in% c("+", "-"))
  
  ## for label, use "gene" if available, otherwise "Name"
  geneName = gff[sel, "gene"]
  featName = gff[sel, "Name"]
  featName[!is.na(geneName)] = geneName[!is.na(geneName)]

  ## split by feature type (e.g. CDS, ncRNA)
  feature  = as.character(gff[sel, "feature"])
  featsp = split(seq(along=sel), feature)

  ## There are now five different cases, and we need to deal with them:
  ## - ignorable features, given by featureExclude
  ## - genes: a horizontal line + name
  ## - introns: a caret
  ## - CDS: a box + no name
  ## - all others: a colored box + name

  ## in this vector we save those features for which we want to have names
  whnames = integer(0)

  ## 1. drop the ignorable ones
  featsp = featsp[ ! (names(featsp) %in% featureExclude) ]
  
  ## 2. gene: just a horizontal line + name
  wh = ("gene" == names(featsp))
  if(any(wh)) {
    i = featsp[["gene"]]
    s = sel[i]
    grid.segments(x0 = gff$start[s], x1 = gff$end[s], y0 = 0, y1 = 0,
                  default.units = "native", gp = gpar(col="#a0a0a0"))
    whnames = i
    featsp = featsp[!wh]
  }

  ## 3.introns
  wh = ("intron" == names(featsp))
  if(any(wh)) {
    i = featsp[["intron"]]
    s = sel[i]
    mid = (gff$start[s]+gff$end[s])/2
    wid = (gff$end[s]-gff$start[s])/2 
    for(z in c(-1,1))
      grid.segments(x0 = mid,
                    x1 = mid+z*wid,
                    y0 = 1.20*c("+"=1, "-"=-1)[strand],  ## istrand is 1 or 2
                    y1 = 0.95*c("+"=1, "-"=-1)[strand],
                    default.units = "native",
                    gp = gpar(col="black"))
     featsp = featsp[!wh]
  } ## if
  
  ## 4. colors for boxes
  ## check that we know how deal with all features
  featCols = featureColors(featureColorScheme)

  whm = names(featsp) %in% rownames(featCols)
  if(!all(whm))
    warning("Don't know how to handle feature of type(s) '", paste(names(featsp)[!whm], collapse=", "), "' in gff.", sep="")

  sfeatsp  = featsp[rownames(featCols)]
  ll       = listLen(sfeatsp)
  
  if(any(ll>0)) {
    i  = unlist(sfeatsp)
    gp = gpar(col = rep(featCols$col,  ll),
                 fill = rep(featCols$fill, ll))
    s  = sel[i]
    grid.rect(x     = gff$start[s],
              y     = 0,
              width = gff$end[s]-gff$start[s],
              height= 2,
              default.units = "native",
              just  = c("left", "center"),
              gp    = gp)
    whnames = c(whnames, unlist(sfeatsp[!(names(sfeatsp) %in% featureNoLabel)]))
    ## additional potentially useful values for featureNoLabel: "binding_site", "TF_binding_site"
  }

  ## labels
  if( !all(tolower(featureNoLabel)=="all") && (length(whnames)>0)) {

    ## this is a bit of a hack to abbreviate the labels of "binding site" features:
    bindingRegexpr = "binding.?site.*$"
    isBindingSite = (regexpr(bindingRegexpr, featName[whnames]) > 0)
    if(any(isBindingSite)) {
      ## replace long labels
      featName[whnames] = gsub(bindingRegexpr, "bs", featName[whnames])
    }

    ## remove duplicated names that are not binding sites
    whnames = whnames[isBindingSite | !duplicated(featName[whnames])]

    txtcex = 0.6
    txtdy  = 0.7
    s      = sel[whnames]
    txtx   = (gff$start[s]+gff$end[s])/2
    txty   = numeric(length(s))
    ord    = order(txtx)
    whnames = whnames[ord]
    s      = s[ord]
    txtx   = txtx[ord]
    
    strw   = convertWidth(stringWidth(featName[whnames]), "native", valueOnly=TRUE)*txtcex
    rightB = txtx[1] + 0.5*strw[1]
    doText = rep(TRUE, length(whnames))

    # not used so far:
    # textstarts <-  txtx - 0.5*strw
    # textends   <-  txtx + 0.5*strw

    # adjust text labels to be still readable in feature-dense areas:
    if(length(whnames) >1) {
      for(k in 2:length(whnames)) {
        leftB = txtx[k] - 0.5*strw[k]
        if(leftB > rightB) { # all texts not overlapping next to each other?
          rightB = txtx[k] + 0.5*strw[k]
        } else { # any overlaps?
          if(!any(txty[k-(1:2)]==txtdy)) {#  2 previous labels not moved up?
            txty[k]= txtdy                #   then this one 
          } else {                        #  else try move down:
            if(!any(txty[k-(1:2)]== -txtdy)) { 
              txty[k]= -txtdy             #  if 2 previous ones weren't
            } else {
              doText[k] = FALSE           #  otherwise don't put the label
            }
          }
        } ##  else
      } ## for
    }
    
    grid.text(label = featName[whnames][doText],
              x = txtx[doText], y = txty[doText], gp=gpar(cex=txtcex), 
              default.units = "native")
  } ## if
  
  popViewport()

} ## plotFeatures

##------------------------------------------------------------
##
##------------------------------------------------------------
alongChromTicks = function(x){
  rx = range(x)
  lz = log((rx[2]-rx[1])/3, 10)
  fl = floor(lz)
  if( lz-fl > log(5, 10))
    fl = fl +  log(5, 10)
  tw = round(10^fl)
  i0 = ceiling(rx[1]/tw)
  i1 = floor(rx[2]/tw)
  seq(i0, i1)*tw
}

##------------------------------------------------------------
## featureColors
## note that features are drawn in the order in which they appear
## here, this can be used to let important features overdraw less
## important ones (e.g. tRNA is more specific than ncRNA)
## to test, say tilingArray:::plotAlongChromLegend()
##------------------------------------------------------------
featureColors = function(scheme=1){
  
  defaultColors = c(
    "chromosome"  = NA,
    "nucleotide_match" = "#e0e0e0",   ## light gray
    "pseudogene"  = "#e0e0e0",        ## light gray
    "uORF"        =   "#FED976" ,     ## orange
    "nc_primary_transcript" = "#a0a0a0",    ## grey
    "region" = "#cc66cc",           ## light red-violet	
    "repeat_family" = "#CC6666",    ## light red
    "repeat_region" = "#e31a1c",    ## bright red                    
    "transposable_element"  = "#f1b6da",    ## pink
    "transposable_element_gene"= "#f1b6da",
    "ARS"         = "#CC9966",    ## light brown
    "centromere"  = "#FFEDA0",    ## orange
    "telomere"    = "#FFEDA0",    ## orange
    "insertion"   = "#FFEDA0",    ## orange
    "CDS"         = "#addfff",    ## light blue
    "CDS_dubious" = "#e0f1f2",    ## lighter blue
    "ncRNA"       = "#a6d96a",    ## green 
    "tRNA"        = "#a6d96a",    ## green
    "snRNA"       = "#8C6BB1",    ## purple
    "rRNA"        = "#fdae61",    ## meat
    "snoRNA"      = "#7F5A58",    ## red brown
    "binding_site"    = "#C9C299", ## lemon chiffon
    "TF_binding_site" = "#C9C299" ## lemon chiffon
  )
  darkenborder = as.logical(c(rep(1,3),0,rep(1, 17),0,0))
  stopifnot(length(darkenborder)==length(defaultColors))
  
  fill = switch(scheme,
    default  = defaultColors,
    unicolor = ifelse(is.na(defaultColors), NA,  "#addfff"),  ## light blue
    stop("Sapperlot"))
  
  ## calculate hex string for a color that is a little bit darker than the
  ## hex string in the argument
  darken = function(x, factor=0.5) {
    wh = which(!is.na(x))

    hex = sapply(x[wh], substring, first=c(2,4,6), last=c(3,5,7))
    hex = apply(hex, 2, function(h) as.integer(factor*as.integer(paste("0x", h, sep=""))))

    res = rep(as.character(NA), length(x))
    res[wh] = apply(hex, 2, function(h) sprintf("#%02x%02x%02x", h[1], h[2], h[3]))
    return(res)
  }
  
  border = ifelse(darkenborder, darken(fill), fill)
  
  res = data.frame(fill=I(fill),
    col =I(border))
  rownames(res)=names(defaultColors) 
  return(res)
} 

#############################################################
##
## Funciones para analizar lincRNAs
##
############################################################

"drawLinc" <- function(lincStr, geneEnsemblData, resDEG, gname="GeneName", gexp="logFC", gstat="B", ustat=max, rangeLow = 5e5, rangeHigh = 5e5, numcol = 10) {
  
  col <- greenred(numcol+1)
  col <- col[-(round(numcol/2)+1)] 
  resDEG2 <- resDEG[,gexp]
  colMax <- ceiling(((resDEG2+max(abs(resDEG2)))/(2*max(abs(resDEG2))))*length(col))
  
  colMax[which(colMax == 0)] <- 1
  
  colTmp <- col[colMax]
  names(colTmp) <- paste(resDEG[,gname])
  
  for(k in 1:length(lincStr))
  {
	if (length(strsplit(lincStr[k], ":")[[1]]) == 4) {
	
		lincLocChr <- substr(strsplit(lincStr[k], ":")[[1]][3], 4, nchar(strsplit(lincStr[k], ":")[[1]][2]))
		lincLocStart <- as.numeric(strsplit(strsplit(lincStr[k], ":")[[1]][4], "-")[[1]][1])
		lincLocEnd <- as.numeric(strsplit(strsplit(strsplit(lincStr[k], ":")[[1]][4], "-")[[1]][2], "_")[[1]][1])
		lincLocStrandChr <- strsplit(strsplit(strsplit(lincStr[k], ":")[[1]][4], "-")[[1]][2], "_")[[1]][2]
	
	} else {
	
		lincLocChr <- substr(strsplit(lincStr[k], ":")[[1]][2], 4, nchar(strsplit(lincStr[k], ":")[[1]][2]))
		lincLocStart <- as.numeric(strsplit(strsplit(lincStr[k], ":")[[1]][3], "-")[[1]][1])
		lincLocEnd <- as.numeric(strsplit(strsplit(strsplit(lincStr[k], ":")[[1]][3], "-")[[1]][2], "_")[[1]][1])
		lincLocStrandChr <- strsplit(strsplit(strsplit(lincStr[k], ":")[[1]][3], "-")[[1]][2], "_")[[1]][2]
	
	}

	geneRegion <- geneEnsemblData[which((geneEnsemblData[, "ChromosomeName"] == lincLocChr) & (((geneEnsemblData[, "ExonChrStart"] > (lincLocStart-rangeLow)) & (geneEnsemblData[, "ExonChrStart"] < (lincLocEnd+rangeHigh))) | ((geneEnsemblData[, "ExonChrEnd"] > (lincLocStart-rangeLow)) & (geneEnsemblData[, "ExonChrEnd"] < (lincLocEnd+rangeHigh))))),]
	geneRegion <- unique(geneRegion[,c("EnsemblGeneID", "GeneName", "ChromosomeName", "GeneStart", "GeneEnd", "Strand", "GeneBiotype")])
	
	if (length(unique(c(which(geneRegion[, "GeneEnd"] <  (lincLocStart-rangeLow)), which(geneRegion[, "GeneStart"] >  (lincLocEnd+rangeHigh)))))>0)
		geneRegion <- geneRegion[-unique(c(which(geneRegion[, "GeneEnd"] <  (lincLocStart-rangeLow)), which(geneRegion[, "GeneStart"] >  (lincLocEnd+rangeHigh)))),]
	
	# par(mfrow=c(3,1))
	
	gffChrF <- geneRegion[geneRegion[,"Strand"] == 1,]
	gffChrR <- geneRegion[geneRegion[,"Strand"] == -1,]
	xlim <- c(min(geneRegion[,"GeneStart"]), max(geneRegion[,"GeneStart"]))
	
	logFCTmp <- resDEG[paste(resDEG[,gname]) %in% c(paste(geneRegion[, "GeneName"]), lincStr[k]), c(gname, gexp)]
	colnames(logFCTmp) <- c("GeneName","Exp")
	logFCTmp <- summaryBy(. ~ GeneName, data = logFCTmp, FUN = mean, keep.names = TRUE)
	fcData <- logFCTmp[,2]
	names(fcData) <- logFCTmp[,1]
	BTmp <- resDEG[paste(resDEG[,gname]) %in% c(paste(geneRegion[, "GeneName"]), lincStr[k]), c(gname, gstat)]
	colnames(BTmp) <- c("GeneName","Stat")
	BTmp <- summaryBy(. ~ GeneName, data = BTmp, FUN = ustat, keep.names = TRUE)
	bData <- BTmp[,2]
	names(bData) <- BTmp[,1]
			
	plot(range(sum(nrow(gffChrF), nrow(gffChrR))), c(1, 1), xlim = xlim, ylim = c(-2.5, 2.5), main = paste(lincStr[k],gstat,"=", round(bData[lincStr[k]], 3), gexp,"=", round(fcData[lincStr[k]], 3), sep = " "), xlab = "", ylab = "", axes = FALSE, type = "n")
	
	axis(1, pos = 0)
	for (i in 1:nrow(gffChrF)) {
	
		if (nrow(gffChrF)>0) {
			if (is.na(colTmp[paste(gffChrF[i, "GeneName"])])) {
				colThis <- "white"
				bThis = NA
				fcThis = NA
			} else {
				colThis <- colTmp[paste(gffChrF[i, "GeneName"])]
				bThis <- bData[paste(gffChrF[i, "GeneName"])]
				fcThis <- fcData[paste(gffChrF[i, "GeneName"])]
			}
			
			if (gffChrF[i, "GeneBiotype"] == "protein_coding") {
				rect(gffChrF[i, "GeneStart"], 0.15, gffChrF[i, "GeneEnd"], 0.3, col = colThis)
			} else {
				rect(gffChrF[i, "GeneStart"], 0.35, gffChrF[i, "GeneEnd"], 0.5, col = colThis)
			}
			strT <- paste(gffChrF[i, "GeneName"], gffChrF[i, "GeneBiotype"], round(bThis, 3), round(fcThis, 3), sep = "_")
			text(round((gffChrF[i, "GeneStart"]+gffChrF[i, "GeneEnd"])/2), 0.75, strT, srt = 90, pos = 4)
		}
			
	}
	
	for (i in 1:nrow(gffChrR)) {
	
		if (nrow(gffChrR)>0) {
			if (is.na(colTmp[paste(gffChrR[i, "GeneName"])])) {
				colThis <- "white"
				bThis = NA
				fcThis = NA
			} else {
				colThis <- colTmp[paste(gffChrR[i, "GeneName"])]
				bThis <- bData[paste(gffChrR[i, "GeneName"])]
				fcThis <- fcData[paste(gffChrR[i, "GeneName"])]
			}
			
			if (gffChrR[i, "GeneBiotype"] == "protein_coding") {
				rect(gffChrR[i, "GeneStart"], -0.3, gffChrR[i, "GeneEnd"], -0.15, col = colThis)
			} else {
				rect(gffChrR[i, "GeneStart"], -0.5, gffChrR[i, "GeneEnd"], -0.35, col = colThis)
			}
			
			strT <- paste(gffChrR[i, "GeneName"], gffChrR[i, "GeneBiotype"], round(bThis, 3), round(fcThis, 3), sep = "_")
			text(round((gffChrR[i, "GeneStart"]+gffChrR[i, "GeneEnd"])/2), -0.75, strT, srt = 270, pos = 4)
		}
			
	}
	
	if (is.na(colTmp[lincStr[k]])) {
		colThis <- "white"
	} else {
		colThis <- colTmp[paste(lincStr[k])]
	}

	rect(lincLocStart, 0.1, lincLocEnd, -0.1, col = colThis)
  }
}


#############################################################
#
# Funciones para analizar TLDAS - miRNAs
#
#############################################################

"qnNormalize" <- 
function (eData)
{
    if ( is.null(dim(eData)) | ncol(eData) == 1) return(eData)

    foo <- function( y ) seq(0,1, length.out=length(y))

    vector.sort <- function( x ){
    	x.sort <- sort(x, method = "quick") 
	if( length(x.sort) == length(x) ) return(x.sort)    ## If there are no null values
	approx( foo( x.sort ), x.sort, foo( x ), ties = "ordered" )$y
    }

    m <- rowMeans( apply( eData, 2, vector.sort ) )
    i <- foo( m )

    vector.approx <- function( x ){
    	R <- rank( x , na = NA)
	R <- ( R - min(R) ) / (max(R) - min(R))
        x[!is.na(x)] <- approx(i, m, R, ties = "ordered")$y
	x
    }

    apply( eData, 2, vector.approx )
}

"repmat" <- function(a,n,m) {kronecker(matrix(1,n,m),a)}

"genorm" <- function( data )
{
	if( nrow( data ) == 1 ) 
		stop( "Dataset contains a single housekeeping candidate gene." )

    # Transforming data; see notes above for reference
    data <- - sweep( data, 1, apply( data, 1, min ) )

	# auxiliary functions

	calculate.M.score <- function( data ){        # Calculates value M in genorm paper (see page 11)
        cov.data <- cov( t(data) , use="all")	#Eli: se puede añadir use=complete para que calcule cov sin tener en cuenta NAs (pero funciona mal)
        sapply( 1:nrow( cov.data ), function(i) mean( sqrt( cov.data[i,i] + diag( cov.data )[-i] - 2 * cov.data[i,-i] ) ) )
	}

	new.synth.gene <- function( data ){
		colMeans( data ) 
	}
	
	# preparing output
	gene.names <- rownames( data )
	results  <- matrix(NA, 0, length( gene.names ) )
	new.data <- matrix(NA, 0, ncol(data))
	colnames( results ) <- gene.names

	# iterations
	genes.keep <- 1:length(gene.names)
	while( length( genes.keep ) > 1 ){
		data.tmp <- data[genes.keep, ]

		M.score <- calculate.M.score(data.tmp)

		new.row <- rep(NA, ncol(results))
		new.row[genes.keep] <- M.score
		results <- rbind( results, new.row )
		new.data <- rbind( new.data, new.synth.gene(data.tmp) )

		genes.keep <- genes.keep[-which.max( M.score )[1]]
	}

	# transforming output
	rownames( results ) <- paste("iter", 1:nrow(results), sep=".")
	#results <- results[,order( apply( results, 2, function(x) sum(is.na(x)) ) )]
	#data.frame( 
	#		results, 
	#		total.score = apply( results, 1, mean, na.rm = TRUE ),
	#		incrementa.score = c( apply( new.data[-nrow(new.data),, drop=FALSE] - new.data[-1,, drop=FALSE], 1, sd), NA )
	#)
	order( apply( results, 2, function(x) sum(is.na(x)) ) )
}

#############################################################
#
# Funciones para analizar CNV y Genotipo (arrays de SNPs)
#
#############################################################

genotypeEmission <- function(genotypes, states, probHomCall, probMissing, verbose = TRUE) {
	if (!is.numeric(genotypes)) 
		stop("genotypes must be integers (1=AA, 2=AB, 3=BB, 4=missing") 
		
	emissionForGenotypes <- function(probHomGenotype, genotypes) {
	
		isHom <- which(as.vector(genotypes) == 1 | as.vector(genotypes) == 3)
		isHet <- which(as.vector(genotypes) == 2)
		isMissing <- which(as.vector(genotypes) == 4 | is.na(as.vector(genotypes)))
		emission.gt <- rep(NA, length(genotypes))
		emission.gt[isHom] <- probHomGenotype
		emission.gt[isHet] <- 1 - probHomGenotype
		emission.gt[isMissing] <- NA
		
		emission.gt
	}

	emission.gt <- array(NA, dim = c(nrow(GT), ncol(GT), length(states)))
	for (j in 1:ncol(GT)) {
		emission.gt[, j, ] <- sapply(probs, emissionForGenotypes, genotypes = GT[, j])
		if (any(is.na(emission.gt[, j, 1]))) {
			missing <- is.na(emission.gt[, j, 1])
			if (!missing(probMissing)) {
				if (length(probMissing) != length(states))
					stop("probMissing must be a numeric vector equal to the number of states")
				emission.gt[missing, j, ] <- matrix(probMissing, sum(missing), length(states), byrow = TRUE)
			} else {
				if (verbose)
					message("Argument probMissing is not specified. Assume that missing genotype calls are independent of the underling hidden state")
					emission.gt[missing, j, ] <- 1
			}
		}
	}
	dimnames(emission.gt) <- list(rownames(genotypes), colnames(genotypes), states)
	return(suppressWarnings(log(emission.gt)))
}

"cnvHMM" <- function(dataCNV, chr, sampleName, annotData, dataOrd, dataGT, states, cnlocation, phomCall, plotCNV = TRUE) {

	CNVOrd <- dataCNV[dataOrd]
	CT <- CNVOrd[which(annotData$Chromosome == chr)]
	
	# lfilter = 50
	# CTFilter <- filter(CT, rep(1, lfilter), method = "convolution")/lfilter

	CTRaw <- 2*as.matrix(2^(CT))
	colnames(CTRaw) <- sampleName
	sample.sd <- matrix(rowSds(t(log2(CTRaw)), na.rm = T), nrow(CTRaw), ncol(CTRaw))		
	
	ann <- annotData[which(annotData$Chromosome == chr),c(3, 2)]
	rownames(ann) <- annotData[which(annotData$Chromosome == chr), 1]
	colnames(ann) <- c("chromosome", "position")
	ann[, "chromosome"] <- integer2chromosome(ann[, "chromosome"])
	fD <- new("AnnotatedDataFrame", data = ann, varMetadata = data.frame(lableDescription = colnames(ann)))
	pD <- annotatedDataFrameFrom(CTRaw, byrow = FALSE)
	
	GT <- dataGT[dataOrd, sample] 
	GT <- as.matrix(GT[which(annotData$Chromosome == chr)])
	GT <- matrix(as.integer(GT), nrow(GT), ncol(GT))
	dimnames(GT) <- dimnames(CTRaw)
	
	dataChr <- new("oligoSnpSet", copyNumber = CTRaw, calls = GT, featureData = fD, phenoData = pD)
	options <- new("HmmOptions", snpset = dataChr, states = states, copyNumber.location = cnlocation, copyNumber.scale = sample.sd[1], probHomCall = phomCall)
	params <- new("HmmParameter", states = states(options), initialStateProbability = 0.999)
	
	cn.emission <- copyNumber.emission(options)
	gt.emission <- calls.emission(options)
	emission(params) <- cn.emission+gt.emission
	genomicDistance(params) <- exp(-2 * diff(position(dataChr))/(100 * 1e6))
	transitionScale(params) <- matrix(1, length(states), length(states))
	hmmpredict <- hmm(options, params)
	breaks <- findBreaks(predictions(hmmpredict), states = states, position = ann[, "position"], chromosome = chr, sample = sampleName)
	
	if (plotCNV) {
	
		gp <- plot(snpset(options), hmmpredict)
		gp$abline.v <- TRUE
		allParameters <- unlist(snpPar(gp))
		gp$col.predict[2] <- "white"		
		ylim <- 0.5
		if (gp$ylim[1] < 0.5) ylim <-  gp$ylim[1] 
		gp$ylim <- c(ylim, 10)
		gp$cytoband.ycoords <- c(ylim+0.1, ylim+0.05)
		gp$hmm.ycoords <- c(ylim, ylim-0.1)
		show(gp)
		# lines(ann[, "position"], predictions(hmmpredict))
		legend(-0.05, 10, fill = gp$col.predict, legend = states, bty = "n", title = "predicted states")
						
	}
	
	return(breaks)
	
}

#############################################################
#
# Funciones para análisis de Supervivencia
#
#############################################################

"analysisSurvival" <- function(dataMat, dataClinic, geneNames, geneSelection, timeToEvent = "Survival", timeUnits="Months", event = "Died" ,method = "median", q = 0.25, legendx=20, legendy=1, legendUnits="exp",interval_tab=FALSE, pdf = FALSE, pdf_size=8, marks=TRUE, filename = "surv.pdf", maxStatPlot=FALSE) {

	indsel <- which(geneNames %in% geneSelection)
	dataSel <- dataMat[indsel,]
	
	res <- list()
	
	if (pdf) {
		pdf(file = filename, colormode = "rgb", pdf_size, pdf_size)
	}
	
	if (method == "maxstat") {
	
		for (i in 1:nrow(dataSel)) {
			x <- as.numeric(paste(dataSel[i,]))
			time <- as.numeric(paste(dataClinic[,timeToEvent]))
			cens <- as.numeric(paste(dataClinic[,event]))
			mydata <- data.frame(cbind(time, cens, x))
			
			mod1 <- maxstat.test(Surv(time, cens) ~ x, data=mydata, smethod="LogRank", pmethod="exactGauss", abseps=0.01)
			
			x <- (x>mod1$estimate)*1
			mydata <- data.frame(cbind(time, cens, x))
			fit <- survfit(Surv(time, cens) ~ x, data=mydata)
			diff <- survdiff(Surv(time, cens) ~ x, data=mydata)
			pval <- 1-pchisq(diff$chisq, 1)
			
			if (pdf) {
				if(maxStatPlot){
					plot(mod1, xlab="Expression")
					abline(v=mean(as.numeric(paste(dataSel[i,]))))
				}
				titlePS <- paste(rownames(dataMat)[indsel[i]], "-", geneNames[indsel[i]]," (Th = ", round(mod1$estimate, 5), ", p-value=", round(pval,5) ,")")
				plot(fit, xlab=timeUnits, ylab=paste(event,"Rate",sep=" "), main = titlePS, lty=c(2, 1), mark.time=marks)
				legend(legendx, legendy, paste(c("Low", "High"), legendUnits), lty=c(2,1))
			} else{
				if(maxStatPlot){
					plot(mod1, xlab="Expression")
					abline(v=mean(as.numeric(paste(dataSel[i,]))))
				}
				titlePS <- paste(rownames(dataMat)[indsel[i]], "-", geneNames[indsel[i]]," (Th = ", round(mod1$estimate, 5), ", p-value=", round(pval,5) ,")")
				plot(fit, xlab=timeUnits, ylab=paste(event,"Rate",sep=" "), main = titlePS, lty=c(2, 1), mark.time=marks)
				legend(legendx, legendy, paste(c("Low", "High"), legendUnits), lty=c(2,1))
			}
			
			res[[i]] <- pval						
			names(res)[i] <- paste(rownames(dataMat)[indsel[i]], geneNames[indsel[i]], sep = " - ")
			
		}
		
	} else {
		
		for (i in 1:nrow(dataSel)) {
			dataSurvGene.i <- as.numeric(paste(dataSel[i,]))
			fun <- function(data,f,qf) { if(f=="quantile"){ quantile(data, prob=q, na.rm=TRUE) } else { eval(as.name(f))(data, na.rm=TRUE) }}
			dataSurv.i <- cbind(as.numeric(paste(dataClinic[ ,timeToEvent])),as.numeric(paste(dataClinic[ ,event])),
			(dataSurvGene.i > fun(dataSurvGene.i,method,q))*1)
			dataSurv.i <- as.data.frame(dataSurv.i)
			colnames(dataSurv.i) <- c(timeToEvent, event, "Exp")
			fit <- survfit(Surv(dataSurv.i[,timeToEvent], dataSurv.i[,event]) ~ Exp, data=dataSurv.i)
			diff <- survdiff(Surv(dataSurv.i[,timeToEvent], dataSurv.i[,event]) ~ Exp, data=dataSurv.i)
			pval <- 1-pchisq(diff$chisq, 1)
			if (pdf) {
				if(interval_tab)
				{
					risk.data <- data.frame(strata = summary(fit, times = 0, extend = TRUE)$strata, time = summary(fit, times = 0, extend = TRUE)$time, n.risk = summary(fit, times = 0, extend = TRUE)$n.risk)[,-2]
					times <- seq(from=10, to=60, by=10)
					for(j in 1:length(times))
					{
					risk.data <- cbind(risk.data, data.frame(strata = summary(fit, times = times[j], extend = TRUE)$strata, time = summary(fit, times = times[j], extend = TRUE)$time, n.risk = summary(fit, times = times[j], extend = TRUE)$n.risk)[,-c(1:2)])
					}
					colnames(risk.data)[2:ncol(risk.data)] <- c(0,times)
					risk.data[,1] <- c("Low","High")
					par(mar=c(12,8,4,2)+0.1)
					titlePS <- paste(rownames(dataSel)[i], " - ", geneNames[indsel[i]], "(p-value=", round(pval,5) ,")")
					plot(fit, xlab=timeUnits, ylab=paste(event,"Rate",sep=" "), main = titlePS, lty=c(2, 1), xlim=c(0,63), mark.time=marks)
					legend(legendx, legendy, paste(c("Low", "High"),legendUnits), lty=c(2,1))
					mtext(text=paste(risk.data[1,]), side=1, at=c(-5,0,times),line=6)
					mtext(text=paste(risk.data[2,]), side=1, at=c(-5,0,times),line=8)
				}else{
					titlePS <- paste(rownames(dataSel)[i], " - ", geneNames[indsel[i]], "(p-value=", round(pval,5) ,")")
					plot(fit, xlab=timeUnits, ylab=paste(event,"Rate",sep=" "), main = titlePS, lty=c(2, 1), mark.time=marks)
					legend(legendx, legendy, paste(c("Low", "High"),legendUnits), lty=c(2,1))
				}
			} else{
				titlePS <- paste(rownames(dataMat)[indsel[i]], " - ", geneNames[indsel[i]], "(p-value=", round(pval,5) ,")")
				plot(fit, xlab=timeUnits, ylab=paste(event,"Rate",sep=" "), main = titlePS, lty=c(2, 1), mark.time=marks)
				legend(legendx, legendy, paste(c("Low", "High"),legendUnits), lty=c(2,1))
			} 
			res[[i]] <- pval
			names(res)[i] <- paste(rownames(dataMat)[indsel[i]], geneNames[indsel[i]], sep = " - ")
			
		}
		
	}
			
	if (pdf) {
		dev.off()
	} 
	
	return(unlist(res))

}

#’ Create a Kaplan-Meier plot using ggplot2
#’
#’ @param sfit a \code{\link[survival]{survfit}} object
#’ @param table logical: Create a table graphic below the K-M plot, indicating at-risk numbers?
#’ @param returns logical: if \code{TRUE}, return an arrangeGrob object
#’ @param xlabs x-axis label
#’ @param ylabs y-axis label
#’ @param ystratalabs The strata labels. \code{Default = levels(summary(sfit)$strata)}
#’ @param ystrataname The legend name. Default = “Strata”
#’ @param timeby numeric: control the granularity along the time-axis
#’ @param main plot title
#’ @param pval logical: add the pvalue to the plot?
#’ @return a ggplot is made. if return=TRUE, then an arrangeGlob object
#’ is returned
#’ @author Abhijit Dasgupta with contributions by Gil Tomas
#’ \url{http://statbandit.wordpress.com/2011/03/08/an-enhanced-kaplan-meier-plot/}
#’ @export
#’ @examples
#’ \dontrun{
#’ data(colon)
#’  fit <- survfit(Surv(time,status)~rx, data=colon)
#'  ggkm(fit, timeby=500)
#' }

"ggkm" <- function(sfit, table = TRUE, returns = FALSE, xlabs = "Time", ylabs = "survival probability", ystratalabs = NULL, ystrataname = NULL, timeby = 100, main = "Kaplan-Meier Plot", pval = TRUE, legendPos=c(0.85,0.90), strataPos=-0.20) {
	
	require(ggplot2)
	require(survival)
	require(gridExtra)
	require(plyr)
	
	if (is.null(ystratalabs)) {
		ystratalabs <- as.character(levels(summary(sfit)$strata))
	}
	
	m <- max(nchar(ystratalabs))
	
	if (is.null(ystrataname)) ystrataname <- "Strata"
	
	times <- seq(0, max(sfit$time), by = timeby)
	.df <- data.frame(time = sfit$time, n.risk = sfit$n.risk,
    	n.event = sfit$n.event, surv = sfit$surv, strata = summary(sfit, censored = T)$strata,
    	upper = sfit$upper, lower = sfit$lower)
	levels(.df$strata) <- ystratalabs
	zeros <- data.frame(time = 0, surv = 1, strata = factor(ystratalabs, levels=levels(.df$strata)), upper = 1, lower = 1)
	.df <- rbind.fill(zeros, .df)
	d <- length(levels(.df$strata))
	
	p <- ggplot(.df, aes(time, surv, group = strata)) + geom_step(aes(linetype = strata), size = 0.7) +  theme_bw() + scale_x_continuous(xlabs, breaks = times, limits = c(0, max(sfit$time))) + scale_y_continuous(ylabs, limits = c(0, 1)) + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), legend.position = c(legendPos[1], legendPos[2]), legend.background = element_rect(size=0.5, color="black"), legend.key=element_rect(color="white"), plot.margin = unit(c(2, 1, 1, ifelse(m < 10, 2.5,3.5)), "lines"), axis.title.x = element_text(vjust = 0.5)) + labs(linetype = ystrataname) + ggtitle(main)
 
	## Create a blank plot for place-holding
	## .df <- data.frame()

	blank.pic <- ggplot(.df, aes(time, surv)) + geom_blank() + theme_bw() + theme(axis.text.x = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), panel.grid.major = element_blank(), panel.border = element_blank())

	if (pval) {
	
    		sdiff <- survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    		pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
    		pvaltxt <- ifelse(pval < 0.0001, "p < 0.0001", paste("p =", signif(pval, 3)))
    		p <- p + annotate("text", x = 0.6 * max(sfit$time), y = 0.1, label = pvaltxt)
	
	}
	
	if (table) {
    		
		## Create table graphic to include at-risk numbers
    		risk.data <- data.frame(strata = summary(sfit, times = times, extend = TRUE)$strata, time = summary(sfit, times = times, extend = TRUE)$time, n.risk = summary(sfit, times = times, extend = TRUE)$n.risk) 
	
		data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
        	#, color = strata)) +
        	geom_text(size = 3.5) + theme_bw() + scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = ystratalabs) +
        	# scale_y_discrete(#format1ter = abbreviate,
        	# breaks = 1:3,
        	# labels = ystratalabs) +
        	scale_x_continuous("Numbers at risk", limits = c(0, max(sfit$time))) + theme(axis.title.x = element_text(size = 10, vjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(),     axis.text.y = element_text(face = "bold", hjust = 1))
		
    		data.table <- data.table + theme(legend.position = "none") + xlab(NULL) + ylab(NULL)
    		
		data.table <- data.table + theme(plot.margin = unit(c(-1.5, 1, 0.1, ifelse(m < 10, 2.5, 3.5)-strataPos * m), "lines"))
		
		## Plotting the graphs
		## p <- ggplotGrob(p)
		## p <- addGrob(p, textGrob(x = unit(.8, "npc"), y = unit(.25, "npc"), label = pvaltxt,
		## gp = gpar(fontsize = 12)))
		
		grid.arrange(p, blank.pic, data.table, clip = FALSE, nrow = 3, ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))
		
		if(returns) {
		
			a <- arrangeGrob(p, blank.pic, data.table, clip = FALSE, nrow = 3, ncol = 1, heights = unit(c(2, .1, .25),c("null", "null", "null")))
			return(a)
			
		}
		
	} else {
	
    		## p <- ggplotGrob(p)
    		## p <- addGrob(p, textGrob(x = unit(0.5, "npc"), y = unit(0.23, "npc"),
    		## label = pvaltxt, gp = gpar(fontsize = 12)))
    		print(p)
    		if (returns) return(p)
    	}
	
}

 ###########################################
#
# FUNCIONES PARA EL ANÁLISIS DE EXOMAS
#
###########################################

library(VariantAnnotation)
"VCF_ReadANDclean" <- function(vcfparse, org="hg38", filter=TRUE, GT=FALSE)
{
	print("Reading vcf file...")
	vcf_df_list <- list()
	for(i in 1:length(vcfparse))
	{
		vcfread <- readVcf(vcfparse[i],org)
		vcf_df <- cbind(ID=rownames(as.data.frame(vcfread@rowRanges)),as.data.frame(vcfread@rowRanges), vcfread@fixed, vcfread@info)
		if(filter==TRUE)
		{
			vcf_df <- vcf_df[vcf_df$FILTER=="PASS",]
		}
		rownames(vcf_df) <- vcf_df[,1]

		if(GT==TRUE & length(vcfparse)>1)
		{
			print("GT option only available for the joint VCF.")
		}
		
		if(GT==TRUE & length(vcfparse)==1)
		{
			print("Parsing genotype...")
			GT_filter_cod <- VCF_geno(vcfread, vcf_df, filter)
			print("Parsing vcf information...")
			vcf_df <- cleanVCF_eli(vcf_df)
			return(list(vcf_df=vcf_df,GT_df=GT_filter_cod))

		}else{
			print("Parsing vcf information...")
			vcf_df_list[[i]] <- cleanVCF_eli(vcf_df)
		}
	}
	if(GT==FALSE & length(vcfparse)==1)
	{
		return(vcf_df_list[[1]])
	}
	if(length(vcfparse)>1)
	{
		return(vcf_df_list)
	}
}

"cleanVCF_eli" <- function(vcfparse) {
	for(i in 1:ncol(vcfparse))
	{
		if(class(vcfparse[,i])=="CompressedCharacterList")
		{
			vcfparse[,i] <- unlist(lapply(vcfparse[,i], paste, collapse=";"))
		}
		if(class(vcfparse[,i])=="CompressedIntegerList")
		{
			vcfparse[,i] <- unlist(lapply(vcfparse[,i], paste, collapse=";"))
		}
		if(class(vcfparse[,i])=="CompressedNumericList")
		{
			vcfparse[,i] <- unlist(lapply(vcfparse[,i], paste, collapse=";"))
		}
		if(class(vcfparse[,i])=="DNAStringSetList")
		{
			vcfparse[,i] <- unlist(lapply(vcfparse[,i], paste, collapse=";"))
		}
		if(class(vcfparse[,i])=="DNAStringSet")
		{
			vcfparse[,i] <- paste(width(vcfparse[,i]), paste(vcfparse[,i]) ,sep="_")
		}
		if(length(grep("[\\]", vcfparse[,i]))>0)
		{
			vcfparse[,i] <- gsub("[\\]x3[b-d]",";",paste(vcfparse[,i]))
		}
		print(i)
	}
	return(vcfparse)
}
"VCF_geno" <- function(jointvcf, jointvcf_df, filter){
	GT <- geno(jointvcf)$GT
	if(filter==TRUE){
	 GT_filter <- GT[rownames(jointvcf_df),]
    }else{ 
     GT_filter <- GT
    }
	GT_filter_cod <- GT_filter
	GT_filter_cod[GT_filter == "0/0"] <- 1
	GT_filter_cod[GT_filter == "0/1"] <- 2
	GT_filter_cod[GT_filter == "1/1"] <- 3
	GT_filter_cod[GT_filter == "."] <- 0
	GT_filter_cod[GT_filter == "0/2"] <- 0
	GT_filter_cod[GT_filter == "0/3"] <- 0
	GT_filter_cod[GT_filter == "0/4"] <- 0
	GT_filter_cod[GT_filter == "1/2"] <- 0
	GT_filter_cod[GT_filter == "2/2"] <- 0
	GT_filter_cod[GT_filter == "2/3"] <- 0
	GT_filter_cod[GT_filter == "2/4"] <- 0
	return(GT_filter_cod)
}

################################
#
# NUEVAS FUNCIONES 
#
################################

tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
   type <- match.arg(type)
   assign(".type", type, envir=baseenv())
   if(gcFirst) gc(FALSE)
   tic <- proc.time()[type]         
   assign(".tic", tic, envir=baseenv())
   invisible(tic)
}

toc <- function()
{
   type <- get(".type", envir=baseenv())
   toc <- proc.time()[type]
   tic <- get(".tic", envir=baseenv())
   print(toc - tic)
   invisible(toc)
   return(toc - tic)
}

"parseENCODE" <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[seq(2, 20, 2)]
	return(tmp)

}

"parseENCODE_eli" <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[c(grep("gene_id",tmp)+1,grep("gene_type",tmp)+1,grep("gene_name",tmp)+1,grep("level",tmp)+1, ifelse(sum(grep("havana_gene",tmp))>0, grep("havana_gene",tmp)+1, NA))]
	return(tmp)

}

"parseCUFF" <- function(x) {

	tmp <- unlist(strsplit(unlist(strsplit(x, ";")), " "))
	tmp <- tmp[tmp != ""]
	tmp <- tmp[seq(2, 17, 2)]
	if (sum(is.na(tmp))>0) 
		tmp <- c(tmp[1:3], NA, tmp[4], NA, tmp[5:6])
	return(tmp)
}


"analysisGO" <- function(seldatagenes, annotdata, universeRef, gsGO, analysisname, W = 40, H = 40, pval = 0.01, mincat = 2, nfilter = 4, onepage = TRUE) {

	ENS_seldatagenes <- unique(paste(annotdata[paste(annotdata[,2]) %in% paste(seldatagenes$GeneName), 1]))
	GO_seldatagenes <- goGSEnrich2(intersect(ENS_seldatagenes, universeRef), universeRef, gsGO, annotdata, pval, mincat)
	GO_seldatagenes_Filter <- GO_seldatagenes[which(paste(GO_seldatagenes[,"GOID"]) %in% names(which(table(GO_seldatagenes[,"GOID"])>nfilter))),]
	write.table(GO_seldatagenes_Filter, file=paste("GO_", analysisname, ".txt", sep = ""), quote = FALSE, row.names = TRUE, sep="\t")

	GO_seldatagenes_Filter$logp <- (-1)*log10(as.numeric(paste(GO_seldatagenes_Filter$P.Value)))

	sort1 <- unique(GO_seldatagenes_Filter[GO_seldatagenes_Filter[,2] == "BP", c("GOTerm", "logp")])
	sort2 <- unique(GO_seldatagenes_Filter[GO_seldatagenes_Filter[,2] == "MF", c("GOTerm", "logp")])
	sort3 <- unique(GO_seldatagenes_Filter[GO_seldatagenes_Filter[,2] == "CC", c("GOTerm", "logp")])

	sort1 <- transform(sort1, GOTerm = reorder(GOTerm, logp))
	gg1 <- ggplot(sort1, aes(x = GOTerm, y = as.numeric(paste(logp)))) + geom_bar(stat = "identity") + labs(x = "", y = "-log(pvalue)") + coord_flip() + labs(title = "Biological process") + theme(axis.text.x  = element_text(angle = 90))

	sort2 <- transform(sort2, GOTerm = reorder(GOTerm, logp))
	gg2 <- ggplot(sort2, aes(x = GOTerm, y = as.numeric(paste(logp)))) + geom_bar(stat = "identity") + labs(x = "", y = "-log(pvalue)") + coord_flip() + labs(title = "Molecular process") + theme(axis.text.x  = element_text(angle = 90))

	sort3 <- transform(sort3, GOTerm = reorder(GOTerm, logp))
	gg3 <- ggplot(sort3, aes(x = GOTerm, y = as.numeric(paste(logp)))) + geom_bar(stat = "identity") + labs(x = "", y = "-log(pvalue)") + coord_flip() + labs(title = "Cellular component") + theme(axis.text.x  = element_text(angle = 90)) 

	if (onepage) {
	
		pdf(file = paste("GO_Analysis_", analysisname, ".pdf", sep = ""), height = H, width = W, colormodel = "rgb")
		multiplot(gg1, gg2, gg3, layout = matrix(c(1,2,3), ncol=3))
		dev.off()

	} else {

		pdf(file = paste("GO_Analysis_", analysisname, ".pdf", sep = ""), height = H, width = W, colormodel = "rgb")
		gg1
		gg2
		gg3
		dev.off()

	}
}

"goGSEnrich2" <- function(selSet, refSet, gsSet, geneAnnot, hgCutoff, catSize) {

    paramsBP <- GSEAGOHyperGParams(name = "GOBPTest",
    geneSetCollection = gsSet,
    geneIds = selSet,
    universeGeneIds = refSet,
    ontology="BP",
    pvalueCutoff = hgCutoff,
    conditional = TRUE,
    testDirection = "over")

    print("Analizing gene list. GO Biological Process enrichment ...")

    hgOverBP <- hyperGTest(paramsBP)
    x <- sigCategories(hgOverBP)
GO_BP <- data.frame(GOID=rep(paste(x)[1],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), P.Value=rep(pvalues(hgOverBP)[x[1]],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), OR=rep(oddsRatios(hgOverBP)[x[1]],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), ExpCount=rep(expectedCounts(hgOverBP)[x[1]],length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))), Count=paste(unlist(geneIdsByCategory(hgOverBP)[x[1]])), Size= rep(length(unlist(geneIdUniverse(hgOverBP)[x[1]])),length(unlist(geneIdsByCategory(hgOverBP)[x[1]]))))
for(i in 2:length(x))
{
	GO_BP <- rbind(GO_BP, data.frame(GOID=rep(paste(x)[i],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), P.Value=rep(pvalues(hgOverBP)[x[i]],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), OR=rep(oddsRatios(hgOverBP)[x[i]],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), ExpCount=rep(expectedCounts(hgOverBP)[x[i]],length(unlist(geneIdsByCategory(hgOverBP)[x[i]]))), Count=paste(unlist(geneIdsByCategory(hgOverBP)[x[i]])), Size= rep(length(unlist(geneIdUniverse(hgOverBP)[x[i]])),length(unlist(geneIdsByCategory(hgOverBP)[x[i]])))))
}

    paramsMF <- GSEAGOHyperGParams(name = "GOMFTest",
    geneSetCollection = gsSet,
    geneIds = selSet,
    universeGeneIds = refSet,
    ontology="MF",
    pvalueCutoff = hgCutoff,
    conditional = TRUE,
    testDirection = "over")

    print("Analizing gene list. GO Molecular Function enrichment ...")

    hgOverMF <- hyperGTest(paramsMF)
x <- sigCategories(hgOverMF)
GO_MF <- data.frame(GOID=rep(paste(x)[1],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), P.Value=rep(pvalues(hgOverMF)[x[1]],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), OR=rep(oddsRatios(hgOverMF)[x[1]],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), ExpCount=rep(expectedCounts(hgOverMF)[x[1]],length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))), Count=paste(unlist(geneIdsByCategory(hgOverMF)[x[1]])), Size= rep(length(unlist(geneIdUniverse(hgOverMF)[x[1]])),length(unlist(geneIdsByCategory(hgOverMF)[x[1]]))))
for(i in 2:length(x))
{
	GO_MF <- rbind(GO_MF, data.frame(GOID=rep(paste(x)[i],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), P.Value=rep(pvalues(hgOverMF)[x[i]],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), OR=rep(oddsRatios(hgOverMF)[x[i]],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), ExpCount=rep(expectedCounts(hgOverMF)[x[i]],length(unlist(geneIdsByCategory(hgOverMF)[x[i]]))), Count=paste(unlist(geneIdsByCategory(hgOverMF)[x[i]])), Size= rep(length(unlist(geneIdUniverse(hgOverMF)[x[i]])),length(unlist(geneIdsByCategory(hgOverMF)[x[i]])))))
}


    paramsCC <- GSEAGOHyperGParams(name = "GOCCTest",
    geneSetCollection = gsSet,
    geneIds = selSet,
    universeGeneIds = refSet,
    ontology="CC",
    pvalueCutoff = hgCutoff,
    conditional = TRUE,
    testDirection = "over")

    print("Analizing gene list. GO Cellular Component enrichment ...")

    hgOverCC <- hyperGTest(paramsCC)
    x <- sigCategories(hgOverCC)
GO_CC <- data.frame(GOID=rep(paste(x)[1],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), P.Value=rep(pvalues(hgOverCC)[x[1]],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), OR=rep(oddsRatios(hgOverCC)[x[1]],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), ExpCount=rep(expectedCounts(hgOverCC)[x[1]],length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))), Count=paste(unlist(geneIdsByCategory(hgOverCC)[x[1]])), Size= rep(length(unlist(geneIdUniverse(hgOverCC)[x[1]])),length(unlist(geneIdsByCategory(hgOverCC)[x[1]]))))
for(i in 2:length(x))
{
	GO_CC <- rbind(GO_CC, data.frame(GOID=rep(paste(x)[i],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), P.Value=rep(pvalues(hgOverCC)[x[i]],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), OR=rep(oddsRatios(hgOverCC)[x[i]],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), ExpCount=rep(expectedCounts(hgOverCC)[x[i]],length(unlist(geneIdsByCategory(hgOverCC)[x[i]]))), Count=paste(unlist(geneIdsByCategory(hgOverCC)[x[i]])), Size= rep(length(unlist(geneIdUniverse(hgOverCC)[x[i]])),length(unlist(geneIdsByCategory(hgOverCC)[x[i]])))))
}

GO_res <- rbind(cbind(Ontology="BP",GO_BP[GO_BP$Size>=2,]), cbind(Ontology="MF",GO_MF[GO_MF$Size>=2,]), cbind(Ontology="CC",GO_CC[GO_CC$Size>=2,]))
godb_doterm <- Term(GOTERM)
GO_res_annot <- merge(cbind(GO_res[,1:2],GOTerm=godb_doterm[paste(GO_res[,2])], GO_res[3], GeneID=GO_res[,6]), annotEnsembl, by.x=5, by.y=1)
GO_res_annot

}

#############################################################
##
##      ATAC-Seq
##
#############################################################

# Modifico la función de enrichPlot/clusterProfiler para cambiar algunos parámetros por defecto que tiene y personalizar un poco más el plot
# gsInfo, gseaScores y tableGrob2 no las cambio pero necesito que sea visible en el workspace
# Da error TableGrob2
library(grid)
gsInfo <- function(object, geneSetID) {
    geneList <- object@geneList

    if (is.numeric(geneSetID))
        geneSetID <- object@result[geneSetID, "ID"]

    geneSet <- object@geneSets[[geneSetID]]
    exponent <- object@params[["exponent"]]
    df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
    df$ymin <- 0
    df$ymax <- 0
    pos <- df$position == 1
    h <- diff(range(df$runningScore))/20
    df$ymin[pos] <- -h
    df$ymax[pos] <- h
    df$geneList <- geneList

    df$Description <- object@result[geneSetID, "Description"]
    return(df)
}
gseaScores <- getFromNamespace("gseaScores", "DOSE")
ggtable <- function(d, p = NULL) {
    # has_package("ggplotify")
    ggplotify::as.ggplot(tableGrob2(d, p))
}
tableGrob2 <- function(d, p = NULL) {
    # has_package("gridExtra")
    d <- d[order(rownames(d)),]
    tp <- gridExtra::tableGrob(d)
    if (is.null(p)) {
        return(tp)
    }

    # Fix bug: The 'group' order of lines and dots/path is different
    p_data <- ggplot_build(p)$data[[1]]
    # pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
    p_data <- p_data[order(p_data[["group"]]), ]
    pcol <- unique(p_data[["colour"]])
    ## This is fine too
    ## pcol <- unique(p_data[["colour"]])[unique(p_data[["group"]])]  
    j <- which(tp$layout$name == "rowhead-fg")

    for (i in seq_along(pcol)) {
        tp$grobs[j][[i+1]][["gp"]] <- gpar(col = pcol[i])
    }
    return(tp)
}
gseaplot2_eli <- function(x, geneSetID, title = "", color="green", base_size = 11,
                      rel_heights=c(1.5, .5, 1), subplots = 1:3,
                      pvalue_table = FALSE, ES_geom="line", RankedListMetric="Ranked List Metric") {
    ES_geom <- match.arg(ES_geom, c("line", "dot"))

    geneList <- position <- NULL ## to satisfy codetool

    if (length(geneSetID) == 1) {
        gsdata <- gsInfo(x, geneSetID)
    } else {
        gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }

    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
        theme_classic(base_size) +
        theme(panel.grid.major = element_line(colour = "grey92"),
              panel.grid.minor = element_line(colour = "grey92"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()) +
        scale_x_continuous(expand=c(0,0))

    if (ES_geom == "line") {
        es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                              size=1)
    } else {
        es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                               size=1, data = subset(gsdata, position == 1))
    }

    p.res <- p + es_layer +
        theme(legend.position = c(.8, .8), legend.title = element_blank(),
              legend.background = element_rect(fill = "transparent"))

    p.res <- p.res + ylab("Running Enrichment Score") +
        theme(axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x=element_blank(),
              plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))

    i <- 0
    for (term in unique(gsdata$Description)) {
        idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
        gsdata[idx, "ymin"] <- i
        gsdata[idx, "ymax"] <- i + 1
        i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) +
        geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
        xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
        theme(legend.position = "none",
              plot.margin = margin(t=-.1, b=0,unit="cm"),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              axis.line.x = element_blank()) +
        scale_x_continuous(expand=c(0,0)) +
        scale_y_continuous(expand=c(0,0))

    if (length(geneSetID) == 1) {
        ## geneList <- gsdata$geneList
        ## j <- which.min(abs(geneList))
        ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
        ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]

        ## v <- sort(c(v1, v2))
        ## inv <- findInterval(geneList, v)

        v <- seq(1, sum(gsdata$position), length.out=9)
        inv <- findInterval(rev(cumsum(gsdata$position)), v)
        if (min(inv) == 0) inv <- inv + 1

        col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))

        ymin <- min(p2$data$ymin)
        yy <- max(p2$data$ymax - p2$data$ymin) * .3
        xmin <- which(!duplicated(inv))
        xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
        d <- data.frame(ymin = ymin, ymax = yy,
                        xmin = xmin,
                        xmax = xmax,
                        col = col[unique(inv)])
        p2 <- p2 + geom_rect(
                       aes_(xmin=~xmin,
                            xmax=~xmax,
                            ymin=~ymin,
                            ymax=~ymax,
                            fill=~I(col)),
                       data=d,
                       alpha=.9,
                       inherit.aes=FALSE)
    }

    ## p2 <- p2 +
    ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
    ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
    ## theme(legend.position="none") +
    ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))

    df2 <- p$data #data.frame(x = which(p$data$position == 1))
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                              color="grey")
    p.pos <- p.pos + ylab(RankedListMetric) +
        xlab("Rank in Ordered Dataset") +
        theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))

    if (!is.null(title) && !is.na(title) && title != "")
        p.res <- p.res + ggtitle(title)

    if (length(color) == length(geneSetID)) {
        p.res <- p.res + scale_color_manual(values=color)
        if (length(color) == 1) {
            p.res <- p.res + theme(legend.position = "none")
            p2 <- p2 + scale_color_manual(values = "black")
        } else {
            p2 <- p2 + scale_color_manual(values = color)
        }
    }

    if (pvalue_table) {
        pd <- x[geneSetID, c("Description", "NES","pvalue", "p.adjust")]
        # pd <- pd[order(pd[,1], decreasing=FALSE),]
        rownames(pd) <- pd$Description

        pd <- pd[,-1]
        pd <- format(pd, digits=3)

        tp <- tableGrob2(pd, p.res)

        p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp,xmin = quantile(p.res$data$x, .5),xmax = quantile(p.res$data$x, .95),ymin = quantile(p.res$data$runningScore, .75), ymax = quantile(p.res$data$runningScore, .9))
    }


    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] +
        theme(axis.line.x = element_line(),
              axis.ticks.x=element_line(),
              axis.text.x = element_text())

    if (length(subplots) == 1)
        return(plotlist[[1]] + theme(plot.margin=margin(t=.2, r = .2, b=.2,
                                                        l=.2, unit="cm")))


    if (length(rel_heights) > length(subplots))
        rel_heights <- rel_heights[subplots]

    plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
}

gsearank <- function(x, geneSetID, title="") {
    position <- NULL
    gsdata <- gsInfo(x, geneSetID)
    gsdata <- subset(gsdata, position == 1)
    p <- ggplot(gsdata, aes_(x = ~x, y = ~runningScore)) +
        geom_segment(aes_(xend=~x, yend=0)) +
        ggtitle(title) +
        xlab("Position in the Ranked List of Genes") +
        ylab("Running Enrichment Score") +
        theme_minimal()
    return(p)
}
