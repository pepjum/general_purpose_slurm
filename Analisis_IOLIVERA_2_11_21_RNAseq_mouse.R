library(dplyr)
library(edgeR)
library("genefilter")
library(RColorBrewer)
library(viridis)

source("/home/jgonzalezgom/01_NMELERO/funcionesVikv2.R")


files<-list.files("~/04_IOLIVERA/featureCounts", pattern=".txt$")
files<-paste0("~/04_IOLIVERA/featureCounts/", files)


samples<-c()
data_df<-data.frame()
for(i in 1:length(files)){
    cat(paste(files[i]),"\n")
    fileloaded<-read.table(files[i], header=T)
    sample_tmp<-paste(lapply(strsplit(paste(files[i]),"\\/"),"[",4))
    sample_tmp<-paste(lapply(strsplit(paste(sample_tmp),"Aligned"),"[",1))
    colnames(fileloaded)[7]<-paste(sample_tmp)
    samples<-c(samples, sample_tmp)
    if(i ==1){
        data_df<-fileloaded[,c(1,7)]
    }else{
        #colnames(fileloaded)[7]<-paste(sample_tmp)
        #samples<-c(samples,sample_tmp)
        data_df<-cbind(data_df, fileloaded[,7])
    }

}


colnames(data_df)[2:ncol(data_df)]<-paste(samples)
rownames(data_df)<-data_df$Geneid
data_df<-data_df[,-1]
save(data_df, file="/home/jgonzalezgom/04_IOLIVERA/dataframes/RNASEQ_IOLIVERA_raw.Rdata")


# removing low expressed genes in the samples. We want to keep tags that are expressed in at least one
#of wild-type or treated samples. In either case, the tag should be expressed in at least five libraries.
#We seek tags that achieve one count per million for at least four libraries:

# dim (data_df) 55359    10


#keep <- rowSums(cpm(data_df) > 1) >=5 

#data_df_filtered <- data_df[keep,]

#dim(data_df_filtered) 12060    10

#normalizar de otra forma

# log_intensity_threshold <- 5

# numsamples <- 1   #ajustar para que sea el 50% de muestras de cada tipo
# f1 <- kOverA(numsamples, log_intensity_threshold)
# ffun1 <- filterfun(f1)
# whichFilter.dril_18 <- genefilter(data_df[,c(1:2)], ffun1)

# numsamples <- 1   #ajustar para que sea el 50% de muestras de cada tipo
# f1 <- kOverA(numsamples, log_intensity_threshold)
# ffun1 <- filterfun(f1)
# whichFilter.IL12 <- genefilter(data_df[,c(3:4)], ffun1)

# numsamples <- 1   #ajustar para que sea el 50% de muestras de cada tipo
# f1 <- kOverA(numsamples, log_intensity_threshold)
# ffun1 <- filterfun(f1)
# whichFilter.IL12_DR <- genefilter(data_df[,c(5:6)], ffun1)

# numsamples <- 1   #ajustar para que sea el 50% de muestras de cada tipo
# f1 <- kOverA(numsamples, log_intensity_threshold)
# ffun1 <- filterfun(f1)
# whichFilter.MOCK <- genefilter(data_df[,c(7:8)], ffun1)

# numsamples <- 1   #ajustar para que sea el 50% de muestras de cada tipo
# f1 <- kOverA(numsamples, log_intensity_threshold)
# ffun1 <- filterfun(f1)
# whichFilter.NE <- genefilter(data_df[,c(9:10)], ffun1)

# whichFilter <- (whichFilter.dril_18 | whichFilter.IL12 | whichFilter.IL12_DR | whichFilter.MOCK | whichFilter.NE )

# data_df_filtered <- data_df[whichFilter,]

#### drop in another way




# Normalization
d0<-DGEList(data_df)

d0 <- calcNormFactors(d0)

### filter low expression genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
data_df_filtered<-d$counts
#voom transformation and calculation variance weights

snames<-c(rep("DRIL_18", 2),rep("IL12",2), rep("IL_12_DR",2),rep("MOCK",2),rep("NE",2))  
snames<-factor(snames)

design<-model.matrix(~0+snames)
#rownames(design) <- colnames(d)
#colnames(design) <- levels(snames)

data_Normalized <- voom(d, design, plot = T)
k<-data_Normalized$E   #13225 10

save(k, file="/home/jgonzalezgom/04_IOLIVERA/dataframes/data_Normalized_IOLIVERA_NOV21.Rdata")

#fitting linear models

fit <- lmFit(data_Normalized, design)

#make contrasts
contr <- makeContrasts(DRIL_18vsMOCK= snamesDRIL_18 - snamesMOCK, 
                       DRIL_18_vsNE= snamesDRIL_18 - snamesNE,
                       DRIL_18vsIL12= snamesDRIL_18 - snamesIL12,
                       IL12vsMock=snamesIL12 - snamesMOCK,
                       IL12vsNE= snamesIL12 - snamesNE,
                       IL12_DRvsMock=snamesIL_12_DR - snamesMOCK,
                       IL12_DRvsNE= snamesIL_12_DR - snamesNE,
                       IL12_DRvsIL12=snamesIL_12_DR - snamesIL12,
                       IL12_DRvsDRIL_18= snamesIL_12_DR - snamesDRIL_18,
                       MockvsNE=snamesMOCK - snamesNE, 
                       levels = colnames(coef(fit)))


#estimate contrast for each gene

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

#setwd("/home/jgonzalezgom/04_IOLIVERA/dataframes/DEgenes/")
system("mkdir -p /media/inmuno/05_IOLIVERA/DEgenes")
setwd("/media/inmuno/05_IOLIVERA/DEgenes/")

library(xlsx)
DRIL_18vsMOCK <- topTable(tmp, coef = "DRIL_18vsMOCK", n=nrow(tmp), adjust="fdr")
write.xlsx(DRIL_18vsMOCK, file = "DRIL_18vsMOCK.xlsx")

DRIL_18_vsNE <- topTable(tmp, coef = "DRIL_18_vsNE", n=nrow(tmp), adjust="fdr")
write.xlsx(DRIL_18_vsNE, file = "DRIL_18_vsNE.xlsx")

DRIL_18_vsIL12 <- topTable(tmp, coef = "DRIL_18vsIL12", n=nrow(tmp), adjust="fdr")
write.xlsx(DRIL_18_vsIL12, file = "DRIL_18_vsIL12.xlsx")

IL12vsMock <- topTable(tmp, coef = "IL12vsMock", n=nrow(tmp), adjust="fdr")
write.xlsx(IL12vsMock, file = "IL12vsMock.xlsx")

IL12vsNE <- topTable(tmp, coef = "IL12vsNE", n=nrow(tmp), adjust="fdr")
write.xlsx(IL12vsNE, file = "IL12vsNE.xlsx")

IL12_DRvsIL12 <- topTable(tmp, coef = "IL12_DRvsIL12", n=nrow(tmp), adjust="fdr")
write.xlsx(IL12_DRvsIL12, file = "IL12_DRvsIL12.xlsx")

IL12_DRvsMock <- topTable(tmp, coef = "IL12_DRvsMock", n=nrow(tmp), adjust="fdr")
write.xlsx(IL12_DRvsMock, file = "IL12_DRvsMock.xlsx")

IL12_DRvsNE <- topTable(tmp, coef = "IL12_DRvsNE", n=nrow(tmp), adjust="fdr")
write.xlsx(IL12_DRvsNE, file = "IL12_DRvsNE.xlsx")

IL12_DRvsDRIL18 <- topTable(tmp, coef = "IL12_DRvsDRIL_18", n=nrow(tmp), adjust="fdr")
write.xlsx(IL12_DRvsDRIL18, file = "IL12_DRvsDRIL_18.xlsx")

MockvsNE <- topTable(tmp, coef = "MockvsNE", n=nrow(tmp), adjust="fdr")
write.xlsx(MockvsNE, file = "MockvsNE.xlsx")

# QC control


xlim_r = c(min(na.omit(log2(data_df_filtered+1))),max(na.omit(log2(data_df_filtered+1)))) 
xlim = c(min(na.omit(data_Normalized)),max(na.omit(data_Normalized)))
clust.euclid.average <- hclust(dist(t(log2(data_df_filtered+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(data_Normalized)),method="average")

setwd("~/04_IOLIVERA/images/")

pdf(file="RNASEQ_IOLIVERA_QC2.pdf", colormode="rgb", width=20, height=20)

boxplot(log2(data_df_filtered+1),names=colnames(data_df_filtered), cex.axis=0.7, las=2)
multi("density", log2(data_df_filtered+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(data_Normalized, names=colnames(data_Normalized), cex.axis=0.7, las=2)
multi("density", data_Normalized, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)

dev.off()



### PCA

library(ggplot2)
x="PC1"
y="PC2"
type<-snames
fit <- prcomp(t(na.omit(as.matrix(data_Normalized))), scale=T)
pcData <- data.frame(fit$x)

#ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point( size = 6)  + theme_bw() + theme(legend.title=element_blank())

ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=type, shape=type), size = 5) + scale_shape_manual(values=c(15,18,17,16,20)) + theme_bw() + theme(legend.title=element_blank())

#pdf(file = "PCA_RNASeq_MALVAREZ_without_1-3_3-6_samples_AGO21.pdf", width = 15, height = 15, colormodel = "rgb")
pdf(file = "PCA_RNASeq_IOLIVERA_NOV21.pdf", width = 6, height = 6, colormodel = "rgb")

ggp
dev.off()


# heatmap completo
selection<-c("nombres del paper de Irene")

#mirar ensgs
not_in<-setdiff(selection,k_selection$gene_name)

selection<-c("Nipal1","Cd101","C1qc","Trim9","Lyd","Paqr8","Gja5","Krt83","Apol9b","Il17a","Il17f","Cubn","Ppp1r14c","Kln2","Cd86","Serpinb1a",
"Cxcr6","Crispld2","Susd2","Prr29","Cxcr2","Cxcr1","Ly6g","Asprv1","Peli2","Sema3b","Myh11","Serpinb5","Ffar2","Sptbn2","Lgi1","Sox6","Gng11",
"C230014O12Rik","Mecom","1110002E22Rik","Auts2","Pdzd2","Rassf4","Cd40lg","Igsf10","P2rx3","Cnr1","Msi1","Dnah5","Dgkh","B230110C06Rik","Il1rl1",
"Usp18","Barx1","Prss35","Me1","Ifng","Gm49751","Cdkn2a","Ifitm3","Lrrc3b","Arhgef10l","Onecut2","Iigp1","Gm4951","Ido2","Serpinb2",
"Adamts20","Olfr525","Olfr527","Sval1","9130015A21Rik","Pamr1","Il9","Gpr83","Lamc3","Gm32255","Osr2","Pth","Prss2","Otogl","Slc13a3","St8sia6","Csf2")

library("org.Mm.eg.db") # remember to install it if you don't have it already
symbols <- mapIds(org.Mm.eg.db, keys = selection, keytype = "SYMBOL", column="ENSEMBL")

tmp<-data.frame("symbol"=selection,"ENSEMBL"=symbols)


ensgs<-c("ENSMUSG00000036896","ENSMUSG00000021071","ENSMUSG00000019762","ENSMUSG00000057123","ENSMUSG00000067613","ENSMUSG00000068246","ENSMUSG00000025929",
"ENSMUSG00000041872","ENSMUSG00000026726","ENSMUSG00000040653","ENSMUSG00000043932","ENSMUSG00000031825","ENSMUSG00000009210","ENSMUSG00000026180","ENSMUSG00000022582",
"ENSMUSG00000033508","ENSMUSG00000021846","ENSMUSG00000067242","ENSMUSG00000051910","ENSMUSG00000032766","ENSMUSG00000087461","ENSMUSG00000027684","ENSMUSG00000029673",
"ENSMUSG00000036334","ENSMUSG00000044288","ENSMUSG00000054256","ENSMUSG00000022262","ENSMUSG00000097547","ENSMUSG00000021381","ENSMUSG00000111913","ENSMUSG00000045201",
"ENSMUSG00000031549","ENSMUSG00000062345","ENSMUSG00000022449","ENSMUSG00000061489","ENSMUSG00000062782","ENSMUSG00000029865","ENSMUSG00000113985","ENSMUSG00000027188",
"ENSMUSG00000021538","ENSMUSG00000026840","ENSMUSG00000112120","ENSMUSG00000022330","ENSMUSG00000059077","ENSMUSG00000091455","ENSMUSG00000018459","ENSMUSG00000003418")
                        #cnr1
k_selection<-k_merged[which(k_merged$gene_name %in% selection),]  # names del paper

k_merged$ENSG<-paste(lapply(strsplit(paste(k_merged$ENSG),"\\."),"[",1))

k_selection2<-k_merged[which(k_merged$ENSG %in% ensgs),]
# por ensgs


k_all<-rbind(k_selection,k_selection2)

k_sel<-k[which(k$ENSG %in% tmp$ENSEMBL),]

k_all<-k_all[,c(2:11)]

data_Normalized<-k_all

#transformar en zscore

xt<-t(data_Normalized)
xts<-scale(xt)
data_Normalized<-t(xts)


library(pheatmap)

col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(data_Normalized)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("yellow","red","blue","green","black"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("DRIL_18","IL12", "IL_12_DR", "MOCK" , "NE")

colores<-diverge_hsv(6) #blue, white,red
mat_breaks <- seq(min(data_Normalized), max(data_Normalized), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(data_Normalized), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(data_Normalized)/length(colores), max(data_Normalized), length.out=floor(length(colores)/2)))


pdf("GENES_SELECTED_PAPER_EXPRESSED_IOLIVERA_NOV21.pdf")

pheatmap(data_Normalized,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Zscore GENES Expressed (selection paper St.Paul,2020) ",
    cluster_cols = F,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()

#heatmap de solo genes diferencialmente expresados pvalue 0.05 & B > 0

k<-data_Normalized


DRIL_18vsMOCK_de<-DRIL_18vsMOCK[which(DRIL_18vsMOCK$adj.P.Val < 0.05 ),] 
DRIL_18_vsNE_de<-DRIL_18_vsNE[which(DRIL_18_vsNE$adj.P.Val < 0.05),]
DRIL_18_vsIL12_de<-DRIL_18_vsIL12[which(DRIL_18_vsIL12$adj.P.Val < 0.05),] 
IL12vsMock_de<-IL12vsMock[which(IL12vsMock$adj.P.Val < 0.05 ),] 
IL12vsNE_de<-IL12vsNE[which(IL12vsNE$adj.P.Val < 0.05 ),] 
IL12_DRvsIL12_de<-IL12_DRvsIL12[which(IL12_DRvsIL12$adj.P.Val < 0.05),]
IL12_DRvsMock_de<-IL12_DRvsMock[which(IL12_DRvsMock$adj.P.Val < 0.05 ),] 
IL12_DRvsNE_de<-IL12_DRvsNE[which(IL12_DRvsNE$adj.P.Val < 0.05 ),] 
IL12_DRvsDRIL18_de<-IL12_DRvsDRIL18[which(IL12_DRvsDRIL18$adj.P.Val < 0.05 ),] 
#MockvsNE_de<-MockvsNE[which(MockvsNE$adj.P.Val < 0.05 ),] 


save(DRIL_18vsMOCK_de, DRIL_18_vsNE_de, DRIL_18_vsIL12_de,IL12vsMock_de, IL12vsNE_de,IL12_DRvsIL12_de,IL12_DRvsMock_de, IL12_DRvsNE_de, IL12_DRvsDRIL18, MockvsNE_de , file="/media/inmuno/05_IOLIVERA/data/Deg_normalized_genes_pval_0.05_in_all_conditions.Rdata")


DRIL_18vsMOCK_de_0.01<-DRIL_18vsMOCK[which(DRIL_18vsMOCK$adj.P.Val < 0.01 ),] 
DRIL_18_vsNE_de_0.01<-DRIL_18_vsNE[which(DRIL_18_vsNE$adj.P.Val < 0.01),]
DRIL_18_vsIL12_de_0.01<-DRIL_18_vsIL12[which(DRIL_18_vsIL12$adj.P.Val < 0.01),] 
IL12vsMock_de_0.01<-IL12vsMock[which(IL12vsMock$adj.P.Val < 0.01 ),] 
IL12vsNE_de_0.01<-IL12vsNE[which(IL12vsNE$adj.P.Val < 0.01 ),] 
IL12_DRvsIL12_de_0.01<-IL12_DRvsIL12[which(IL12_DRvsIL12$adj.P.Val < 0.01),]
IL12_DRvsMock_de_0.01<-IL12_DRvsMock[which(IL12_DRvsMock$adj.P.Val < 0.01 ),] 
IL12_DRvsNE_de_0.01<-IL12_DRvsNE[which(IL12_DRvsNE$adj.P.Val < 0.01 ),] 
IL12_DRvsDRIL18_de_0.01<-IL12_DRvsDRIL18[which(IL12_DRvsDRIL18$adj.P.Val < 0.01 ),] 
#MockvsNE_de_0.01<-MockvsNE[which(MockvsNE$adj.P.Val < 0.01 ),] 



save(DRIL_18vsMOCK_de_0.01, DRIL_18_vsNE_de_0.01, DRIL_18_vsIL12_de_0.01,IL12vsMock_de_0.01, IL12vsNE_de_0.01,IL12_DRvsIL12_de_0.01,IL12_DRvsMock_de_0.01, IL12_DRvsNE_de_0.01, IL12_DRvsDRIL18, MockvsNE_de_0.01 , file="/media/inmuno/05_IOLIVERA/data/Deg_normalized_genes_pval_0.01_in_all_conditions.Rdata")

library(ggVennDiagram)

pval_005_ttos_VS_Controles<-list("DRIL_18 vs MOCK"= rownames(DRIL_18vsMOCK_de), "DRIL_18 vs NE "=rownames(DRIL_18_vsNE_de),"IL12 vs Mock"=rownames(IL12vsMock_de), "IL12 vs NE"=rownames(IL12vsNE_de) , "IL12_DR vs Mock"=rownames(IL12_DRvsMock_de), "IL12_DR vs NE"=rownames(IL12_DRvsNE_de))


universo<-unique(c(rownames(DRIL_18vsMOCK_de), rownames(DRIL_18_vsNE_de), rownames(IL12vsMock_de),rownames(IL12vsNE_de), rownames(IL12_DRvsMock_de), "IL12_DR vs NE"=rownames(IL12_DRvsNE_de) ) )

matrix_datos<-matrix(, nrow=length(universo), ncol=6 )
rownames(matrix_datos)<-universo
colnames(matrix_datos)<-c("DRIL18_vs_MOCK","DRIL18_vs_NE","IL12_vs_MOCK","IL12_vs_NE", "IL12_DR_vs_Mock", "IL12_DR_vsNE")



pval_005_ttos_VS_ttos<-list("DRIL_18 vs IL12"=rownames(DRIL_18_vsIL12_de),"IL12_DR vs IL12"=rownames(IL12_DRvsIL12_de),"IL12_DR vs DRIL18"=rownames(IL12_DRvsDRIL18_de))

pval_005_signif_ttos_vs_Controles<-list("DRIL_18 vs Mock"= rownames(DRIL_18vsMOCK_de[which(DRIL_18vsMOCK_de$logFC < -1 | DRIL_18vsMOCK_de$logFC > 1 ),]),"DRIL_18 vs NE"=rownames(DRIL_18_vsNE_de[which(DRIL_18_vsNE_de$logFC < -1 | DRIL_18_vsNE_de$logFC > 1 ),]),
"IL12 vs Mock"=rownames(IL12vsMock_de[which(IL12vsMock_de$logFC < -1 | IL12vsMock_de$logFC > 1 ),]),
"IL12 vs NE"=rownames(IL12vsNE_de[which(IL12vsNE_de$logFC < -1 | IL12vsNE_de$logFC > 1 ),]),
"IL12_DR vs Mock"=rownames(IL12_DRvsMock_de[which(IL12_DRvsMock_de$logFC < -1 | IL12_DRvsMock_de$logFC > 1 ),]), 
"IL12_DR vs NE"=rownames(IL12_DRvsNE_de[which(IL12_DRvsNE_de$logFC < -1 | IL12_DRvsNE_de$logFC > 1 ),]))


pval_005_signif_ttos_vs_ttos<-list("DRIL_18 vs IL12"=rownames(DRIL_18_vsIL12_de[which(DRIL_18_vsIL12_de$logFC < -1 | DRIL_18_vsIL12_de$logFC > 1 ),]),
"IL12_DR vs IL12"=rownames(IL12_DRvsIL12_de[which(IL12_DRvsIL12_de$logFC < -1 | IL12_DRvsIL12_de$logFC > 1 ),]),
"IL12_DR vs DRIL18"=rownames(IL12_DRvsDRIL18_de[which(IL12_DRvsDRIL18_de$logFC < -1 | IL12_DRvsDRIL18_de$logFC > 1 ),]))


pval_001_ttos_VS_Controles<-list("DRIL_18 vs MOCK"= rownames(DRIL_18vsMOCK_de_0.01), 
"DRIL_18 vs NE"=rownames(DRIL_18_vsNE_de_0.01),
"IL12 vs Mock"=rownames(IL12vsMock_de_0.01), 
"IL12 vs NE"=rownames(IL12vsNE_de_0.01) , 
"IL12_DR vs Mock"=rownames(IL12_DRvsMock_de_0.01), 
"IL12_DR vs NE"=rownames(IL12_DRvsNE_de_0.01))

pval_001_ttos_VS_ttos<-list("DRIL_18 vs IL12"=rownames(DRIL_18_vsIL12_de_0.01),
"IL12_DR vs IL12"=rownames(IL12_DRvsIL12_de_0.01),
"IL12_DR vs DRIL18"=rownames(IL12_DRvsDRIL18_de_0.01))


pval_001_signif_ttos_vs_Controles<-list("DRIL_18 vs Mock"= rownames(DRIL_18vsMOCK_de_0.01[which(DRIL_18vsMOCK_de_0.01$logFC < -1 | DRIL_18vsMOCK_de_0.01$logFC > 1 ),]),
"DRIL_18 vs NE"=rownames(DRIL_18_vsNE_de_0.01[which(DRIL_18_vsNE_de_0.01$logFC < -1 | DRIL_18_vsNE_de_0.01$logFC > 1 ),]),
"IL12 vs Mock"=rownames(IL12vsMock_de_0.01[which(IL12vsMock_de_0.01$logFC < -1 | IL12vsMock_de_0.01$logFC > 1 ),]),
"IL12 vs NE"=rownames(IL12vsNE_de_0.01[which(IL12vsNE_de_0.01$logFC < -1 | IL12vsNE_de_0.01$logFC > 1 ),]),
"IL12_DR vs Mock"=rownames(IL12_DRvsMock_de_0.01[which(IL12_DRvsMock_de_0.01$logFC < -1 | IL12_DRvsMock_de_0.01$logFC > 1 ),]), 
"IL12_DR vs NE"=rownames(IL12_DRvsNE_de_0.01[which(IL12_DRvsNE_de_0.01$logFC < -1 | IL12_DRvsNE_de_0.01$logFC > 1 ),]))


pval_001_signif_ttos_vs_ttos<-list("DRIL_18 vs IL12"=rownames(DRIL_18_vsIL12_de_0.01[which(DRIL_18_vsIL12_de_0.01$logFC < -1 | DRIL_18_vsIL12_de_0.01$logFC > 1 ),]),
"IL12_DR vs IL12"=rownames(IL12_DRvsIL12_de_0.01[which(IL12_DRvsIL12_de_0.01$logFC < -1 | IL12_DRvsIL12_de_0.01$logFC > 1 ),]),
"IL12_DR vs DRIL18"=rownames(IL12_DRvsDRIL18_de_0.01[which(IL12_DRvsDRIL18_de_0.01$logFC < -1 | IL12_DRvsDRIL18_de_0.01$logFC > 1 ),]))

library(UpSetR)
detach("package:ComplexUpset", unload=TRUE)

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.05_tratamientos_vs_controles.pdf", height=10, width=10)
upset(fromList(pval_005_ttos_VS_Controles),nsets = 6, order.by="freq")
dev.off()

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.05_tratamientos_vs_controles_significativos.pdf", height=10, width=10)
upset(fromList(pval_005_signif_ttos_vs_Controles),nsets = 6, order.by="freq")
dev.off()

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.05_tratamientos_vs_tratamientos.pdf", height=10, width=10)
upset(fromList(pval_005_ttos_VS_ttos),nsets = 3, order.by="freq")
dev.off()

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.01_tratamientos_vs_tratamientos.pdf", height=10, width=10)
upset(fromList(pval_001_ttos_VS_ttos),nsets = 3, order.by="freq")
dev.off()

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.01_tratamientos_vs_tratamientos_significativos.pdf", height=10, width=10)
upset(fromList(pval_001_signif_ttos_vs_ttos),nsets = 3, order.by="freq")
dev.off()

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.05_tratamientos_vs_tratamientos_significativos.pdf", height=10, width=10)
upset(fromList(pval_005_signif_ttos_vs_ttos),nsets = 3, order.by="freq")
dev.off()


pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.01_tratamientos_vs_controles.pdf", height=10, width=10)
upset(fromList(pval_001_ttos_VS_Controles),nsets = 6, order.by="freq")
dev.off()

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_pval0.01_tratamientos_vs_controles_significativos.pdf", height=10, width=10)
upset(fromList(pval_001_signif_ttos_vs_Controles),nsets = 6, order.by="freq")
dev.off()

### heatmap de los DE en IL12_DR vs DRIL18 pval 0.01 
selection<-data_Normalized[which(rownames(data_Normalized) %in% rownames(IL12_DRvsIL12_de_0.01[which(IL12_DRvsIL12_de_0.01$logFC < -1 | IL12_DRvsIL12_de_0.01$logFC > 1 ),])),]

selection<-data_Normalized[which(rownames(data_Normalized) %in% rownames(IL12_DRvsIL12_de[which(IL12_DRvsIL12_de$logFC < -1 | IL12_DRvsIL12_de$logFC > 1 ),])),]


selection2<-data_Normalized[which(rownames(data_Normalized) %in% rownames(IL12_DRvsMock_de_0.01[which(IL12_DRvsMock_de_0.01$logFC < -1 | IL12_DRvsMock_de_0.01$logFC > 1 ),])),]

selection4<-data_Normalized[which(rownames(data_Normalized) %in% rownames(IL12_DRvsMock_de[which(IL12_DRvsMock_de$logFC < -1 | IL12_DRvsMock_de$logFC > 1 ),])),]

selection4<-selection4$E
selection2<-selection2$E
selection<-selection$E

library(pheatmap)
library(colorspace)

# a zscore
data_scaled<-scale(t(data))
data_ordered<-t(data_scaled)


col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(selection4)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("blue","orange","purple","grey","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("DRIL_18", "IL12" , "IL_12_DR","MOCK","NE")

colores<-diverge_hsv(6) #blue, white,red
mat_breaks <- seq(min(selection4), max(selection4), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(selection4), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(selection4)/length(colores), max(selection4), length.out=floor(length(colores)/2)))

#selection2<-as.data.frame(selection2$E)
#selection_ordered<-selection[order(selection[,8], decreasing=T),]
write.xlsx(selection, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.05_logFC_significativo_IL12_DR_vs_IL12.xlsx")
write.xlsx(selection4, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.05_logFC_significativo_IL12_DR_vs_Mock.xlsx")

write.xlsx(selection2, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_MOCK.xlsx")

setwd("/media/inmuno/05_IOLIVERA/Imagenes/")

pdf("Heatmap_of_normalized_samples_ALL_DEG_pval0.05_significativos_logFC_vsIL12_DR_vs_IL12_ZSCORE_IOLIVERA_NOV21.pdf")

pheatmap(selection4,
    color=colores,
    show_rownames     = FALSE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "DEgenes in IL12_DR vs IL12 \n pval0.05 & logFC < -1 | logFC > 1 ",
    cluster_cols = T,
    cluster_rows = T, 
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


########## Anotación de genes
#DRIL_18vsMOCK

library(rtracklayer)

mouse_gtf<-import("/media/inmuno/references/STAR_mouse_mm39/gencode.vM27.annotation.gtf")

gtf_selection<-data.frame("ENSG"=mouse_gtf$gene_id, "SYMBOL"=mouse_gtf$gene_name)
gtf_selection<-unique(gtf_selection)

DRIL_18vsMOCK_de_0.01$ENSG<-rownames(DRIL_18vsMOCK_de_0.01)
DRIL_18vsMOCK_de_0.01_annot<-merge(DRIL_18vsMOCK_de_0.01, gtf_selection, by="ENSG", all.x=T)

DRIL_18_vsNE_de_0.01$ENSG<-rownames(DRIL_18_vsNE_de_0.01)
DRIL_18_vsNE_de_0.01_annot<-merge(DRIL_18_vsNE_de_0.01, gtf_selection, by="ENSG", all.x=T)

DRIL_18_vsIL12_de_0.01$ENSG<-rownames(DRIL_18_vsIL12_de_0.01) 
DRIL_18_vsIL12_de_0.01_annot<-merge(DRIL_18_vsIL12_de_0.01, gtf_selection, by="ENSG", all.x=T)

IL12vsMock_de_0.01$ENSG<-rownames(IL12vsMock_de_0.01) 
IL12vsMock_de_0.01_annot<-merge(IL12vsMock_de_0.01, gtf_selection, by="ENSG", all.x=T)

IL12vsNE_de_0.01$ENSG<-rownames(IL12vsNE_de_0.01)
IL12vsNE_de_0.01_annot<-merge(IL12vsNE_de_0.01, gtf_selection, by="ENSG", all.x=T)

IL12_DRvsIL12_de_0.01$ENSG<-rownames(IL12_DRvsIL12_de_0.01)
IL12_DRvsIL12_de_0.01_annot<-merge(IL12_DRvsIL12_de_0.01, gtf_selection, by="ENSG", all.x=T)

IL12_DRvsMock_de_0.01$ENSG<-rownames(IL12_DRvsMock_de_0.01)
IL12_DRvsMock_de_0.01_annot<-merge(IL12_DRvsMock_de_0.01,gtf_selection, by="ENSG", all.x=T )

IL12_DRvsNE_de_0.01$ENSG<-rownames(IL12_DRvsNE_de_0.01)
IL12_DRvsNE_de_0.01_annot<-merge(IL12_DRvsNE_de_0.01,gtf_selection, by="ENSG", all.x=T )

IL12_DRvsDRIL18_de_0.01$ENSG<-rownames(IL12_DRvsDRIL18_de_0.01)
IL12_DRvsDRIL18_de_0.01_annot<-merge(IL12_DRvsDRIL18_de_0.01,gtf_selection, by="ENSG", all.x=T )

#DRIL_18vsIL_12_de<-DRIL_18vsIL_12[which(DRIL_18vsIL_12$adj.P.Val < 0.01 ),] 
setwd("/media/inmuno/05_IOLIVERA/DEgenes/")

DRIL_18vsMOCK$ENSG<-rownames(DRIL_18vsMOCK)
DRIL_18vsMOCK_annot<-merge(DRIL_18vsMOCK, gtf_selection, by="ENSG", all.x=T)
write.xlsx(DRIL_18vsMOCK_annot, file="DRIL_18vs_MOCK_ALL_annot.xlsx")

DRIL_18_vsNE$ENSG<-rownames(DRIL_18_vsNE)
DRIL_18_vsNE_annot<-merge(DRIL_18_vsNE, gtf_selection, by="ENSG", all.x=T)
write.xlsx(DRIL_18_vsNE_annot, file="DRIL_18vs_NE_ALL_annot.xlsx")

DRIL_18_vsIL12$ENSG<-rownames(DRIL_18_vsIL12) 
DRIL_18_vsIL12_annot<-merge(DRIL_18_vsIL12, gtf_selection, by="ENSG", all.x=T)
write.xlsx(DRIL_18_vsIL12_annot, file="DRIL_18vs_IL12_ALL_annot.xlsx")

IL12vsMock$ENSG<-rownames(IL12vsMock) 
IL12vsMock_annot<-merge(IL12vsMock, gtf_selection, by="ENSG", all.x=T)
write.xlsx(IL12vsMock_annot, file="IL12vsMock_ALL_annot.xlsx")

IL12vsNE$ENSG<-rownames(IL12vsNE)
IL12vsNE_annot<-merge(IL12vsNE, gtf_selection, by="ENSG", all.x=T)
write.xlsx(IL12vsNE_annot, file="IL12vsNE_ALL_annot.xlsx")

IL12_DRvsIL12$ENSG<-rownames(IL12_DRvsIL12)
IL12_DRvsIL12_annot<-merge(IL12_DRvsIL12, gtf_selection, by="ENSG", all.x=T)
write.xlsx(IL12_DRvsIL12_annot, file="IL12_DRvsIL12_ALL_annot.xlsx")

IL12_DRvsMock$ENSG<-rownames(IL12_DRvsMock)
IL12_DRvsMock_annot<-merge(IL12_DRvsMock,gtf_selection, by="ENSG", all.x=T )
write.xlsx(IL12_DRvsMock_annot, file="IL12_DRvsMock_ALL_annot.xlsx")

IL12_DRvsNE$ENSG<-rownames(IL12_DRvsNE)
IL12_DRvsNE_annot<-merge(IL12_DRvsNE,gtf_selection, by="ENSG", all.x=T )
write.xlsx(IL12_DRvsNE_annot, file="IL12_DRvsNE_ALL_annot.xlsx")

IL12_DRvsDRIL18$ENSG<-rownames(IL12_DRvsDRIL18)
IL12_DRvsDRIL18_annot<-merge(IL12_DRvsDRIL18,gtf_selection, by="ENSG", all.x=T )
write.xlsx(IL12_DRvsDRIL18_de_annot, file="IL12_DRvsDRIL18_ALL_annot.xlsx")


library(ggrepel)
library(dplyr)

volcanoPlots_pepe<-function(dataset, pvalth, FCth, topN ,title){

    dataset<-dataset[,1:8]
    dataset<-unique(dataset)

    dataset <- dataset %>% 
    mutate(
        Expression = case_when(logFC >= FCth & adj.P.Val < pvalth  ~ "Up-regulated",
                            logFC <= -FCth & adj.P.Val < pvalth ~ "Down-regulated",
                            TRUE ~ "Unchanged")
        )



    y_int<-(-1)*log10(pvalth)

    p2 <- ggplot(dataset, aes(logFC, -log(adj.P.Val,10))) +
    geom_point(aes(color = Expression), size = 2/5) +
    xlab(expression("log"[2]*"FC")) + 
    ylab(expression("-log"[10]*"P.Value")) +
    scale_color_manual(values = c("#7CAE00", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
    geom_vline(xintercept=FCth) + geom_vline(xintercept=-(FCth))+
    geom_hline(yintercept=y_int)+
    ggtitle(title)+ theme(plot.title=element_text( size=8, face='bold'))



    top <- topN
    top_genes <- bind_rows(
    dataset %>% 
        filter(Expression == 'Up-regulated') %>% 
        arrange(adj.P.Val, desc(abs(logFC))) %>% 
        head(top),
    dataset %>% 
        filter(Expression == 'Down-regulated') %>% 
        arrange(adj.P.Val, desc(abs(logFC))) %>% 
        head(top)
    )

    p3<-p2+ geom_label_repel(data = top_genes,
                    mapping = aes(logFC, -log(adj.P.Val,10), label = SYMBOL),
                    size = 2, max.overlaps = Inf)
    
    return(p3)
}

DRIL_18vsMOCK_plot<-volcanoPlots_pepe(DRIL_18vsMOCK_annot,0.05,1, 30 ,"DRIL_18 vs MOCK")
DRIL_18_vsNE_plot<-volcanoPlots_pepe(DRIL_18_vsNE_annot,0.05,1, 30, "DRIL_18 vs NE")
DRIL_18_vsIL12_plot<-volcanoPlots_pepe(DRIL_18_vsIL12_annot, 0.05,1,30 ,"DRIL_18 vs IL12")
IL12vsMock_annot_plot<-volcanoPlots_pepe(IL12vsMock_annot,0.05,1,30 ,"IL_12 vs MOCK")
IL12vsNE_annot_plot<-volcanoPlots_pepe(IL12vsNE_annot,0.05,1,30 , "IL_12 vs NE")
IL12_DRvsIL12_annot_plot<-volcanoPlots_pepe(IL12_DRvsIL12_annot, 0.05, 1, 30, "IL12_DR vs IL12")
IL12_DRvsMock_annot_plot<-volcanoPlots_pepe(IL12_DRvsMock_annot,0.05,1,30, "IL_12_DR vs MOCK")
IL12_DRvsNE_annot_plot<-volcanoPlots_pepe(IL12_DRvsNE_annot,0.05,1,30, "IL_12_DR vs NE")
IL12_DRvsDRIL_18annot_plot<-volcanoPlots_pepe(IL12_DRvsDRIL18_annot,0.05,1, 30,"IL12_DR vs DRIL_18")

#setwd("/media/inmuno/05_IOLIVERA/Imagenes/")

pdf("VOLCANO_PLOTS_ALL_pval0.05_top30.pdf")
plot(DRIL_18vsMOCK_plot)
plot(IL12vsMock_annot_plot)
plot(IL12_DRvsMock_annot_plot)
plot(IL12_DRvsNE_annot_plot)
plot(DRIL_18_vsNE_plot)
plot(IL12vsNE_annot_plot)
plot(DRIL_18_vsIL12_plot)
plot(IL12_DRvsDRIL_18annot_plot)
plot(IL12_DRvsIL12_annot_plot)

dev.off()


####### GSEA

library(fgsea)
library(gage)

# ENSg to eNTREZID
#CombovsPBS_annot<-CombovsPBS_annot[which(CombovsPBS_annot$adj.P.Val < 0.05),] #345
#Combovs40_de<-Combovs40[which(Combovs40$adj.P.Val < 0.05 ),] 
#Combovs137_de<-Combovs137[which(Combovs137$adj.P.Val < 0.05 ),] 
#PBSvsNaive_annot<-PBSvsNaive_annot[which(PBSvsNaive_annot$adj.P.Val < 0.05),] # 376
#CombovsNaive_annot<-CombovsNaive_annot[which(CombovsNaive_annot$adj.P.Val < 0.05),] #1796

DRIL_18vsMOCK_annot_order<-DRIL_18vsMOCK_annot[order(DRIL_18vsMOCK_annot$t, decreasing=T),]
DRIL_18_vsNE_annot_order<-DRIL_18_vsNE_annot[order(DRIL_18_vsNE_annot$t, decreasing=T),]
DRIL_18_vsIL12_annot_order<-DRIL_18_vsIL12_annot[order(DRIL_18_vsIL12_annot$t, decreasing=T),]
IL12vsMock_annot_order<-IL12vsMock_annot[order(IL12vsMock_annot$t, decreasing=T),]
IL12vsNE_annot_order<-IL12vsNE_annot[order(IL12vsNE_annot$t, decreasing=T),]
IL12_DRvsIL12_annot_order<-IL12_DRvsIL12_annot[order(IL12_DRvsIL12_annot$t, decreasing=T),]
IL12_DRvsMock_annot_order<-IL12_DRvsMock_annot[order(IL12_DRvsMock_annot$t, decreasing=T),]
IL12_DRvsNE_annot_order<-IL12_DRvsNE_annot[order(IL12_DRvsNE_annot$t, decreasing=T),]
IL12_DRvsDRIL18_annot_order<-IL12_DRvsDRIL18_annot[order(IL12_DRvsDRIL18_annot$t, decreasing=T),]

# ordenar por el estadistico 

preparedata_GSEA<-function(dataframe){
    listed_rank<-list()
    for(i in 1:nrow(dataframe)){
        symbol<-dataframe[,"SYMBOL"][i]
        #symbol<-toupper(symbol)
        t_val<-dataframe[,"t"][i]
        listed_rank[[symbol]]<-t_val
    }
    listed_rank2<-unlist(listed_rank)
    return(listed_rank2)
}

DRIL_18vsMOCK_prepared<-preparedata_GSEA(DRIL_18vsMOCK_annot_order)
DRIL_18_vsNE_prepared<-preparedata_GSEA(DRIL_18_vsNE_annot_order)
DRIL_18_vsIL12_prepared<-preparedata_GSEA(DRIL_18_vsIL12_annot_order)
IL12vsMock_prepared<-preparedata_GSEA(IL12vsMock_annot_order)
IL12vsNE_prepared<-preparedata_GSEA(IL12vsNE_annot_order)
IL12_DRvsIL12_prepared<-preparedata_GSEA(IL12_DRvsIL12_annot_order)
IL12_DRvsMock_prepared<-preparedata_GSEA(IL12_DRvsMock_annot_order)
IL12_DRvsNE_prepared<-preparedata_GSEA(IL12_DRvsNE_annot_order)
IL12_DRvsDRIL18_prepared<-preparedata_GSEA(IL12_DRvsDRIL18_annot_order)

save(DRIL_18vsMOCK_prepared, DRIL_18_vsNE_prepared, DRIL_18_vsIL12_prepared, IL12vsMock_prepared, IL12vsNE_prepared,IL12_DRvsIL12_prepared , IL12_DRvsMock_prepared , IL12_DRvsNE_prepared,IL12_DRvsDRIL18_prepared, file="/media/inmuno/05_IOLIVERA/listas_preparadas_GSEA_todas_comparativas.Rdata" )


###### ESTO ESTA MAL. HUMANO pero el resto bien. Voy a coger referencias de ratón 


pathways <- gmtPathways("/media/inmuno/05_IOLIVERA/m5.go.bp.v0.2.symbols.gmt")   # biocarta y demas
pathways<- gmtPathways("/media/inmuno/05_IOLIVERA/m2.all.v0.2.symbols.gmt")
# canonical pathways curated
pathways<- gmtPathways("/media/inmuno/04_ABELLA/DEgenes_NUEVO/c2.cp.v7.4.symbols.gmt")

#examplePathways<- pathways 

fgseaResIL12_DRvsIL12 <- fgsea(pathways = pathways, stats    = IL12_DRvsIL12_prepared , maxSize=500)

### hacer heatmap GLUCOSA METABOLISM m2 category

#############
data_Normalized<-get(load("/media/inmuno/05_IOLIVERA/data/data_Normalized_IOLIVERA_NOV21.Rdata"))
data_Normalized<-as.data.frame(data_Normalized)
data_Normalized$ENSG<-rownames(data_Normalized)
data_Normalized_annot<-merge(data_Normalized,gtf_selection, by.x="ENSG", by.y="gene_id")

path1<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="REACTOME_GLUCOSE_METABOLISM"),]
genes_path1<-unlist(path1$leadingEdge)

selection_glucose_metabolism<-data_Normalized_annot[which(data_Normalized_annot$SYMBOL %in% genes_path1),]

selection_atp<-data_Normalized_annot[which(data_Normalized_annot$gene_id %in% genes_atp),]


rownames(selection_glucose_metabolism)<-selection_glucose_metabolism$SYMBOL
selection_glucose_metabolism$SYMBOL<-NULL
selection_glucose_metabolism$ENSG<-NULL

selection_glucose_metabolism<-selection_glucose_metabolism[,c(5,6,3,4,1,2,7,8,9,10)]

xt<-t(selection_glucose_metabolism)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)


path2<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="REACTOME_GLUCONEOGENESIS"),]
genes_path2<-unlist(path2$leadingEdge)

selection_gluconeogenesis<-data_Normalized_annot[which(data_Normalized_annot$SYMBOL %in% genes_path2),]

rownames(selection_gluconeogenesis)<-selection_gluconeogenesis$SYMBOL
selection_gluconeogenesis$SYMBOL<-NULL
selection_gluconeogenesis$ENSG<-NULL

selection_gluconeogenesis<-selection_gluconeogenesis[,c(5,6,3,4,1,2,7,8,9,10)]


xt<-t(selection_gluconeogenesis)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)

col_groups<-c(rep("IL12_DR",2),rep("IL12",2),rep("DRIL_18",2),rep("MOCK",2),rep("NE",2))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("purple","orange","blue","gray","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("IL12_DR", "IL12","DRIL_18","MOCK","NE")

pdf("Heatmap_zscored_ATP_SYNTHESIS_COUPLED_PROTON_TRANSPORT_IL12DRvsIL12.pdf",4,4)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    show_colnames     =FALSE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 5,
    main              = "IL12DR vs IL12 'ATP_SYNTHESIS_COUPLED_PROTON_TRANSPORT' \n GSEA ANALYSIS",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=5, 
    cellwidth = 10,

)
dev.off()

cmpr<-list(c("IL12","IL12DR"),c("DRIL18","IL12DR"),c("IL12DR","MOCK"))

pdf("Opa1_gene.pdf")
ggplot(opa1_dft, aes(x =group, y = Expression, color=group)) +
  geom_boxplot() + ggtitle("Opa1 quantitative expression") + stat_compare_means(method="t.test",comparisons = cmpr, tip.length=0.01) + theme(legend.position="none")
dev.off()




selection_glycosylation<-data_Normalized_annot[which(data_Normalized_annot$gene_name %in% genes_glycosilation),]

rownames(selection_glycosylation)<-selection_glycosylation$gene_name
selection_glycosylation$gene_name<-NULL
selection_glycosylation$ENSG<-NULL

selection_glycosylation<-selection_glycosylation[,c(5,6,3,4,1,2,7,8,9,10)]

xt<-t(selection_glycosylation)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)


path2<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="REACTOME_GLUCONEOGENESIS"),]
genes_path2<-unlist(path2$leadingEdge)

selection_gluconeogenesis<-data_Normalized_annot[which(data_Normalized_annot$SYMBOL %in% genes_path2),]

rownames(selection_gluconeogenesis)<-selection_gluconeogenesis$SYMBOL
selection_gluconeogenesis$SYMBOL<-NULL
selection_gluconeogenesis$ENSG<-NULL

selection_gluconeogenesis<-selection_gluconeogenesis[,c(5,6,3,4,1,2,7,8,9,10)]


xt<-t(selection_gluconeogenesis)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)

col_groups<-c(rep("IL12_DR",2),rep("IL12",2),rep("DRIL_18",2),rep("MOCK",2),rep("NE",2))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("purple","orange","blue","gray","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("IL12_DR", "IL12","DRIL_18","MOCK","NE")

pdf("Heatmap_zscored_GOBP_GLYCOSILATION_IL12DRvsIL12.pdf",4,4)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    show_colnames     =FALSE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 5,
    main              = "IL12DR vs IL12 'GOBP GLYCOSILATION' \n GSEA ANALYSIS",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=5, 
    cellwidth = 10,

)
dev.off()


path5<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="MIR_155_5P"),]
path6<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="MIR_155_3P"),]

genes_path5<-unlist(path5$leadingEdge)
genes_path6<-unlist(path6$leadingEdge)

genes_path5<-c(genes_path5,"Socs1")
genes_path5<-c(genes_path5,genes_path6)


selection_mir1553p<-data_Normalized_annot[which(data_Normalized_annot$SYMBOL %in% genes_path5),]

rownames(selection_mir1553p)<-selection_mir1553p$SYMBOL
selection_mir1553p$SYMBOL<-NULL
selection_mir1553p$ENSG<-NULL

selection_mir1553p<-selection_mir1553p[,c(5,6,3,4,1,2,7,8,9,10)]


xt<-t(selection_mir1553p)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)

col_groups<-c(rep("IL12_DR",2),rep("IL12",2),rep("DRIL_18",2),rep("MOCK",2),rep("NE",2))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("purple","orange","blue","gray","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("IL12_DR", "IL12","DRIL_18","MOCK","NE")

pdf("Heatmap_zscored_mir155_targets_IL12DRvsIL12.pdf",4,4)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    show_colnames     =FALSE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 5,
    main              = "IL12DR vs IL12 'Mir155 targets' \n GSEA ANALYSIS",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=5, 
    cellwidth = 10,

)
dev.off()



path3<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="WP_GLYCOLYSIS_AND_GLUCONEOGENESIS"),]
genes_path3<-unlist(path3$leadingEdge)

selection_glyco<-data_Normalized_annot[which(data_Normalized_annot$SYMBOL %in% genes_path3),]

rownames(selection_glyco)<-selection_glyco$SYMBOL
selection_glyco$SYMBOL<-NULL
selection_glyco$ENSG<-NULL

selection_glyco<-selection_glyco[,c(5,6,3,4,1,2,7,8,9,10)]
selection_glyco2<-selection_glyco2[,c(5,6,3,4,1,2,7,8,9,10)]


xt<-t(selection_glyco2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)

col_groups<-c(rep("IL12_DR",2),rep("IL12",2),rep("DRIL_18",2),rep("MOCK",2),rep("NE",2))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("purple","orange","blue","gray","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("IL12_DR", "IL12","DRIL_18","MOCK","NE")

pdf("Heatmap_zscored_GOBP_ENRICHMENT_GLYCOSILATION_IL12DRvsIL12.pdf",4,4)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    show_colnames     =FALSE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 5,
    main              = "IL12DR vs IL12 'BP GLYCOSYLATION' \n ENRICHMENT GOBP ANALYSIS",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=5, 
    cellwidth = 10,

)
dev.off()

path4<-fgseaResIL12_DRvsIL12[which(fgseaResIL12_DRvsIL12$pathway=="REACTOME_THE_CITRIC_ACID_TCA_CYCLE_AND_RESPIRATORY_ELECTRON_TRANSPORT"),]
genes_path4<-unlist(path4$leadingEdge)

selection_mito<-data_Normalized_annot[which(data_Normalized_annot$SYMBOL %in% genes_path4),]

rownames(selection_mito)<-selection_mito$SYMBOL
selection_mito$SYMBOL<-NULL
selection_mito$ENSG<-NULL

selection_mito<-selection_mito[,c(5,6,3,4,1,2,7,8,9,10)]


xt<-t(selection_mito)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)

col_groups<-c(rep("IL12_DR",2),rep("IL12",2),rep("DRIL_18",2),rep("MOCK",2),rep("NE",2))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("purple","orange","blue","gray","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("IL12_DR", "IL12","DRIL_18","MOCK","NE")

pdf("Heatmap_zscored_mitochondrion_localization_IL12DRvsIL12.pdf",4,4)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    show_colnames     =FALSE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 5,
    main              = "IL12DR vs IL12 'BP mitochondrion localization' \n ENRICHMENT GOBP ANALYSIS",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=5, 
    cellwidth = 10,

)
dev.off()




topPathwaysUp_IL12_DRvsIL12 <- fgseaResIL12_DRvsIL12[ES > 0][head(order(pval), n=200), pathway]
topPathwaysDown_IL12_DRvsIL12 <- fgseaResIL12_DRvsIL12[ES < 0][head(order(pval), n=200), pathway]
#topPathways_CombovsNaive <- c(topPathwaysUp_CombovsNaive, rev(topPathwaysDown_CombovsNaive))

topPathways_IL12_DRvsIL12 <- c(topPathwaysUp_IL12_DRvsIL12, rev(topPathwaysDown_IL12_DRvsIL12))

write.xlsx(topPathways_IL12_DRvsIL12, file="top_50UP_and_50Down_Pathways_IL12_DRvsIL12_canonical_pathways.xlsx")
write.xlsx(topPathways_IL12_DRvsIL12, file="top_50UP_and_50Down_Pathways_IL12_DRvsIL12_GENE_ONTOLOGY.xlsx")
write.xlsx(topPathways_IL12_DRvsIL12, file="top_50UP_and_50Down_Pathways_IL12_DRvsIL12_inmuneSigDB.xlsx")


fgseaResIL12_DRvsMock <- fgsea(pathways = pathways, stats    = IL12_DRvsMock_prepared , maxSize=500)
topPathwaysUp_IL12_DRvsMock <- fgseaResIL12_DRvsMock[ES > 0][head(order(pval), n=30), pathway]
topPathwaysDown_IL12_DRvsMock <- fgseaResIL12_DRvsMock[ES < 0][head(order(pval), n=30), pathway]


topPathways_DRvsMock <- c(topPathwaysUp_IL12_DRvsMock, rev(topPathwaysDown_IL12_DRvsMock))

write.xlsx(topPathways_DRvsMock, file="top_50UP_and_50Down_Pathways_IL12_DRvsMock_inmuneSigDB.xlsx")
write.xlsx(topPathways_DRvsMock, file="top_50UP_and_50Down_Pathways_IL12_DRvsMock_GeneOntology_BP.xlsx")
write.xlsx(topPathways_DRvsMock, file="top_50UP_and_50Down_Pathways_IL12_DRvsMock_canonical_pathways.xlsx")


#### with gage
plot_GSEA_PEPE<-function(fgseaRes, gene_set_prepared, topPathways_up, topPathways_down, pathways_selected, sizemin, sizemax, n_plot, title_plot){

## Filter FGSEA by using gage results. Must be significant and in same direction to keep 

    gaRes = gage::gage(gene_set_prepared, gsets=pathways_selected, same.dir=TRUE, set.size =c(sizemin,sizemax))
 
    keepups= topPathways_up 
    keepdowns =  topPathways_down

    selection_top<-fgseaRes[which(fgseaRes$pathway %in% keepups),]
    selection_top<-selection_top[order(selection_top$pval),]
    selection_down<-fgseaRes[which(fgseaRes$pathway %in% keepdowns),] 
    selection_down<-selection_down[order(selection_down$pval, decreasing=T),]

    fgseaRes_sel =rbind(selection_top,selection_down)


    fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$NES > 0, "Up-regulated", "Down-regulated")
    filtRes = rbind(head(fgseaRes_sel, n = n_plot),
                  tail(fgseaRes_sel, n = n_plot ))


    upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
    downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
    colos = c(upcols, downcols)
    names(colos) = 1:length(colos)
    filtRes$Index = as.factor(1:nrow(filtRes))

    g_fgseaRes = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
      geom_col( aes(fill = Index )) +
      scale_fill_manual(values = colos ) +
      coord_flip() +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=title_plot) + 
      theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 20)) 

    return(g_fgseaRes)
}

IL12_DRvsMock_plot<-plot_GSEA_PEPE(fgseaResIL12_DRvsMock, IL12_DRvsMock_prepared, topPathwaysUp_IL12_DRvsMock, topPathwaysDown_IL12_DRvsMock, pathways, 1, 500, 10, "GSEA c7 inmuneSigDB IL12_DR vs Mock" )


IL12_DRvsIL12_plot<-plot_GSEA_PEPE(fgseaResIL12_DRvsIL12, IL12_DRvsIL12_prepared, topPathwaysUp_IL12_DRvsIL12, topPathwaysDown_IL12_DRvsIL12, pathways, 1, 500, 10, "GSEA c2 canonical pathways IL12_DR vs IL12" )


pdf("GSEA_IL12_DRvsMock_c7_inmuneSigDB.pdf",30,30)
plot(IL12_DRvsMock_plot)
dev.off()

pdf("GSEA_IL12_DRvs_Mock_c7_inmuneSigDB_plots.pdf")
for(i in 1:length(IL12_DRvsMock_plot$data$pathway)){
  a<- plotEnrichment(pathways[[IL12_DRvsMock_plot$data$pathway[i]]], IL12_DRvsMock_prepared) + labs(title=IL12_DRvsMock_plot$data$pathway[i])
  plot(a)  
}
dev.off()

gene_signature_typeI_interferon<-pathways$GOBP_POSITIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION

gene_signature_typeI_interferon<-tolower(gene_signature_typeI_interferon)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

gene_signature_typeI_interferon2<-paste(lapply(gene_signature_typeI_interferon, FUN=function(x) firstup(x) ))

write.table(sort(gene_signature_typeI_interferon2), file="typeIinterferon_symbols_GO.txt", col.names=F, row.names=F, quote=F, sep="\t")


#heatmap typeI interferon

typeI_interferon<-read.table("/media/inmuno/05_IOLIVERA/data/Mouse....", header=T,sep="\t")

#typeI_interferon_sel<-typeI_interferon[which(typeI_interferon$Input.Type=="current symbol"),]
typei_interferon_sel2<-typeI_interferon[,c("Ensembl.ID","Symbol")]

## a partir del il12_DR_VS_MOCK

#selection4<-selection2
selection4<-as.data.frame(selection4)
selection4$ENSG<-rownames(selection4)
selection4$ENSG<-paste(lapply(strsplit(paste(selection4$ENSG),"\\."),"[",1))

selection_sel<-selection4[which(selection4$ENSG %in% typeI_interferon$Ensembl.ID),]
selection_sel_annot<-merge(selection_sel, typei_interferon_sel2,by.x="ENSG", by.y="Ensembl.ID", all.x=T)

rownames(selection_sel_annot)<-selection_sel_annot$Symbol

selection_sel_annot<-selection_sel_annot[,c(2:11)]
col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(selection_sel_annot)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("blue","orange","purple","grey","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("DRIL_18", "IL12" , "IL_12_DR","MOCK","NE")

colores<-diverge_hsv(10) #blue, white,red
mat_breaks <- seq(min(selection_sel_annot), max(selection_sel_annot), length.out = 10)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(selection_sel_annot), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(selection_sel_annot)/length(colores), max(selection_sel_annot), length.out=floor(length(colores)/2)))

#selection2<-as.data.frame(selection2$E)
#selection_ordered<-selection[order(selection[,8], decreasing=T),]
#write.xlsx(selection, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_IL12.xlsx")
#write.xlsx(selection2, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_MOCK.xlsx")

setwd("/media/inmuno/05_IOLIVERA/Imagenes/")

pdf("Heatmap_of_normalized_samples_ALL_DEG_pval0.05_related_to_typeI_interferon_production_logFC_vsIL12_DR_MOCK_IOLIVERA_NOV21.pdf")

pheatmap(selection_sel_annot,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "DEgenes in IL12_DR vs Mock \n related to TypeI interferon production \n pval0.05 & logFC < -1 | logFC > 1 ",
    cluster_cols = T,
    cluster_rows = T, 
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


### otro dia cargo lista IL12DR vs IL12 para hacer enrich GO
library(clusterProfiler)
library(org.Mm.eg.db)


ids <- bitr(IL12DRvsIL12$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

a<-IL12DRvsIL12
IL12DRvsIL12_m<-merge(a, ids, by.x="ENSG", by.y="ENSEMBL")

enriched_IL12DRvsIL12_sel<-na.omit(unique(IL12DRvsIL12_m[which(IL12DRvsIL12_m$adj.P.Val<0.05 & IL12DRvsIL12_m$logFC >1), "ENTREZID"]))

enriched_GOBP_IL12DRvsIL12 <- enrichGO(enriched_IL12DRvsIL12_sel, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.05)
enriched_GOBP_IL12DRvsIL12_s <- simplify(enriched_GOBP_IL12DRvsIL12, cutoff=0.7, by="p.adjust", select_fun=min)   # UP

enriched_GOBP_IL12DRvsIL12_df<-as.data.frame(enriched_GOBP_IL12DRvsIL12)
enriched_GOBP_IL12DRvsIL12_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_IL12DRvsIL12_df$p.adjust)))

ids2 <- bitr(genes_glyco2, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo


enriched_IL12DRvsIL12_sel_down<-na.omit(unique(IL12DRvsIL12_m[which(IL12DRvsIL12_m$adj.P.Val<0.01 & IL12DRvsIL12_m$logFC < -1), "ENTREZID"]))

enriched_GOBP_IL12DRvsIL12_down <- enrichGO(enriched_IL12DRvsIL12_sel_down, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_IL12DRvsIL12_down_s <- simplify(enriched_GOBP_IL12DRvsIL12_down, cutoff=0.7, by="p.adjust", select_fun=min)   # DOWN

enriched_GOBP_IL12DRvsIL12_down_df<-as.data.frame(enriched_GOBP_IL12DRvsIL12_down)
enriched_GOBP_IL12DRvsIL12_down_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_IL12DRvsIL12_down_df$p.adjust)))


library(viridis)
library(ggplot2)


d <- GOSemSim::godata('org.Mm.eg.db', ont="BP")



####----
#mydotplot <- clusterProfiler::dotplot(enriched_GOBP_IL12DRvsIL12_s, showCategory=10)


# simMatrix <- rrvgo::calculateSimMatrix(enriched_GOBP_IL12DRvsIL12_df$ID,
#                                        orgdb="org.Mm.eg.db",
#                                        ont="BP",
#                                        method="Rel", 
#                                        semdata = d)

# scores <- setNames(-log10(enriched_GOBP_IL12DRvsIL12_df$qvalue), enriched_GOBP_IL12DRvsIL12_df$ID)

# reducedTerms <- rrvgo::reduceSimMatrix(simMatrix_UP,
#                                 scores,
#                                 threshold=0.98, 
#                                 orgdb="org.Mm.eg.db")


# heat_rrvgo <- heatmapPlot(simMatrix,
#             reducedTerms,
#             annotateParent=TRUE,
#             annotationLabel="parentTerm",
#             fontsize=12,
#             color=viridis(15),
#             silent=T)



# ordenar por LOGP y coger el top10 UP
sort1_up <- enriched_GOBP_IL12DRvsIL12_df[order(enriched_GOBP_IL12DRvsIL12_df[, "logp"], decreasing = TRUE), ]
sort1_up_sel<-sort1_up[1:10,]

sort1_down <- enriched_GOBP_IL12DRvsIL12_down_df[order(enriched_GOBP_IL12DRvsIL12_down_df[, "logp"], decreasing = TRUE), ]
sort1_down_sel<-enriched_GOBP_IL12DRvsIL12_down_df[1:10,]

sort1_down_sel$logp<-sort1_down_sel$logp *(-1)


fgseaRes_sel =rbind(sort1_up_sel,sort1_down_sel)


fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$logp > 0, "Up-regulated", "Down-regulated")
filtRes = fgseaRes_sel

upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
colos = c(upcols, downcols)
names(colos) = 1:length(colos)
filtRes$Index = as.factor(1:nrow(filtRes))

title_plot<-"ENRICHMENT GOBP IL12DR vs IL12"

g_fgseaRes = ggplot(filtRes, aes(reorder(Description, logp), logp)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="logP",
        title=title_plot) + 
    theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 20)) 


pdf("enrichment.pdf",40,40)
plot(g_fgseaRes)
dev.off()

blue.bold.italic.0.01.text <- element_text(face = "bold.italic", color = "blue", size = 20)

pdf('GOAnalysis_coherent_DEG_B30_NOV20_BP_best.pdf', 40,40)
 ggplot(data=unique(GO_cohB30_Filter[GO_cohB30_Filter[,1] == "BP",c("GOTerm","logp")]), aes(x = GOTerm, y = logp)) + geom_bar(stat="identity") + labs(x = "", y = "-log(pvalue)") + theme(axis.text.y = blue.bold.italic.0.01.text) + ggtitle("Biological Process")+ coord_flip() + xlim(unique(paste(sort1_sel[,3]))) 

dev.off()


Cytotoxic effector infiltration (Teff) was estimated by using a 6-gene signature (GZMA, GZMB, PRF1, IFNγ, EOMES, and CD8A), as described in previous publications.7,9

https://pubmed.ncbi.nlm.nih.gov/26755520/
https://pubmed.ncbi.nlm.nih.gov/25428504/

GZMA, GZMB, PRF1, IFNγ, EOMES, and CD8A


matriz_genes<-read.xlsx("matriz_genes_pheatmap_pval_0.05_logFC_significativo_IL12_DR_vs_IL12.xlsx", sheetIndex=1)


matriz_genes$ENSG<-paste(lapply(strsplit(paste(matriz_genes$NA.),"\\."),"[",1))

ids <- bitr(matriz_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

matriz_genes_annot<-merge(matriz_genes, ids, by.x="ENSG", by.y="ENSEMBL", all.x=T)

ensgs<-c("ENSMUSG00000023132", "ENSMUSG00000015437","ENSMUSG00000037202","ENSMUSG00000055170", "ENSMUSG00000032446","ENSMUSG00000053977" )

matriz_genes_annot_sel<-matriz_genes_annot[which(matriz_genes_annot$ENSG %in% ensgs),]

#Gzma, Eomes, IFNg

data_Normalized<-get(load("/media/inmuno/05_IOLIVERA/data/data_Normalized_IOLIVERA_NOV21.Rdata"))

rownames(data_Normalized)<-paste(lapply(strsplit(paste(rownames(data_Normalized)),"\\."),"[",1))
data_Normalized_sel<-data_Normalized[which(rownames(data_Normalized) %in% ensgs),]

rownames(data_Normalized_sel)<-c("Cd8a","Eomes","Prf1","Ifng","Gzma","Gzmb")


snames<-c(rep("DRIL-18",2),rep("IL12",2),rep("IL12-DR",2),rep("Mock",2),rep("NE",2))
col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(data_Normalized_sel)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("blue","orange","purple","grey","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("DRIL-18", "IL12" , "IL12-DR","Mock","NE")

colores<-diverge_hsv(12) #blue, white,red
mat_breaks <- seq(min(data_Normalized_sel), max(data_Normalized_sel), length.out = 12)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(data_Normalized_sel), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(data_Normalized_sel)/length(colores), max(data_Normalized_sel), length.out=floor(length(colores)/2)))

#selection2<-as.data.frame(selection2$E)
#selection_ordered<-selection[order(selection[,8], decreasing=T),]
#write.xlsx(selection, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_IL12.xlsx")
#write.xlsx(selection2, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_MOCK.xlsx")

setwd("/media/inmuno/05_IOLIVERA/Imagenes/")

pdf("Heatmap_Teff_IOLIVERA_NOV21.pdf")

pheatmap(data_Normalized_sel,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "DEgenes in IL12_DR vs IL12 in red \n Genes related to Teff infiltration \n pval0.05 ",
    cluster_cols = T,
    cluster_rows = T, 
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


pathways<- gmtPathways("/media/inmuno/04_ABELLA/DEgenes_NUEVO/c8")

### lista seleccionada
TRAVAGLINI_LUNG_CD8_MEMORY_EFFECTOR_T_CELL

a<-paste(pathways$TRAVAGLINI_LUNG_CD8_MEMORY_EFFECTOR_T_CELL)
b<-paste(pathways$TRAVAGLINI_LUNG_CD4_MEMORY_EFFECTOR_T_CELL)

c_list<-c(a,b)
c_list<-unique(c_list)  # 92 genes
c_list<-tolower(c_list)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

c_list<-firstup(c_list)

not<-setdiff(c_list, ids$SYMBOL)

ids <- bitr(c_list, fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo

ids2<-bitr(ids$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")

ids2<-rbind(ids2,faltantes)

output<-data.frame()

for(ensembl in unique(ids2$ENSEMBL)){
    sel<-ids2[ids2$ENSEMBL==ensembl,]
    if(nrow(sel)>1){
        sel<-sel[1,]
        output<-rbind(output,sel)
    }else{
        sel<-sel
        output<-rbind(output,sel)
    }
}


# lista de 90 genes que tienen que ver con Teff, de ellos 86 tienen expresión y 22 están diferencialmente expresados en el tratamiento combinado

#diferencialmente expresados en IL12DRvsIL12

matriz_genes_annot

matriz_genes_annot_sel<-matriz_genes_annot[which(matriz_genes_annot$ENSG %in% output$ENSEMBL),]

data_Normalized_sel<-data_Normalized[which(rownames(data_Normalized) %in% output$ENSEMBL),]
data_Normalized_sel_annot<-merge(data_Normalized_sel, output, by.x="ENSG", by.y="ENSEMBL", all.x=T)

snames<-c(rep("DRIL-18",2),rep("IL12",2),rep("IL12-DR",2),rep("Mock",2),rep("NE",2))
col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(data_Normalized_sel)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("blue","orange","purple","grey","yellow"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("DRIL-18", "IL12" , "IL12-DR","Mock","NE")

colores<-diverge_hsv(8) #blue, white,red
mat_breaks <- seq(min(data_Normalized_sel), max(data_Normalized_sel), length.out = 8)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(data_Normalized_sel), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(data_Normalized_sel)/length(colores), max(data_Normalized_sel), length.out=floor(length(colores)/2)))

#selection2<-as.data.frame(selection2$E)
#selection_ordered<-selection[order(selection[,8], decreasing=T),]
#write.xlsx(selection, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_IL12.xlsx")
#write.xlsx(selection2, file="/media/inmuno/05_IOLIVERA/DEgenes/matriz_genes_pheatmap_pval_0.01_logFC_significativo_IL12_DR_vs_MOCK.xlsx")

setwd("/media/inmuno/05_IOLIVERA/Imagenes/")

pdf("Heatmap_Teff_IOLIVERA_NOV21_TEFF.pdf")

pheatmap(data_Normalized_sel,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 6,
    main              = "Genes signatures Teff \n Genes related to Teff infiltration \n pval0.05 ",
    cluster_cols = T,
    cluster_rows = T, 
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()




