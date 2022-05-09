library(dplyr)
library(edgeR)
library("genefilter")
library(RColorBrewer)
library(viridis)

source("/home/jgonzalezgom/01_NMELERO/funcionesVikv2.R")


files<-list.files("~/05_ABELLA/featureCounts", pattern=".txt$")
files<-paste0("~/05_ABELLA/featureCounts/", files)


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
system("mkdir -p /home/jgonzalezgom/05_ABELLA/dataframes/")
save(data_df, file="/home/jgonzalezgom/05_ABELLA/dataframes/RNASEQ_ABELLA_unordered_raw.Rdata")

samples_137<-which(grepl("-137", samples))
samples_40<-which(grepl("-40", samples))
samples_combo<-which(grepl("Combo", samples))
samples_naive<-which(grepl("Naive", samples))
samples_ova<-which(grepl("OVA", samples))
samples_PBS<-which(grepl("PBS", samples))

data_df_ordered<-data_df[,c(samples_137,samples_40, samples_combo, samples_naive, samples_ova, samples_PBS)]
save(data_df_ordered, file="/home/jgonzalezgom/05_ABELLA/RNA_SEQ_OCT_21/dataframes/RNASEQ_ABELLA_ordered_raw.Rdata")


`%!in%` <- Negate(`%in%`)
unselect<-c("R2BN-OVA","R1BN-OVA","R1BN-137","R3BN-40","R4PBS","R3Naive","R3BN-Combo")
data_df_ordered_sel<-data_df_ordered[,which(colnames(data_df_ordered) %!in% unselect)]

save(data_df_ordered_sel, file="/media/inmuno/04_ABELLA/RNASEQ_ABELLA_ordered_raw_selection.Rdata")


data_df_ordered<-data_df_ordered_sel

# me paso a trabajar en la workstation
source("/home/inmuno/Documents/scripts_JG/funcionesVikv2.R")


data_df<-get(load("/media/inmuno/04_ABELLA/RNASEQ_ABELLA_ordered_raw.Rdata"))
#data_df<-data_df_ordered

selection2<-data_df[,c(9:16, 19:22)]

save(selection2, file="/media/inmuno/04_ABELLA/RNASEQ_ABELLA_ordered_selection_without_others_treatments.Rdata")

selection<-selection2[,c(1:2,4,6:8,9:11)]
save(selection, file="/media/inmuno/04_ABELLA/RNASEQ_ABELLA_ordered_selection_without_others_treatments_without_outliers.Rdata")

data_df<-selection

# Normalization

selection_combo<-colnames(data_df)[7:9]
selection_pbs<-colnames(data_df)[13:15]
selection_naive<-colnames(data_df)[10:12]

selection_cd137<-colnames(data_df)[1:3]
selection_cd40<-colnames(data_df)[4:6]



log_intensity_threshold<-5

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.combo<-genefilter(data_df[,which(colnames(data_df) %in% selection_combo)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.naive<-genefilter(data_df[,which(colnames(data_df) %in% selection_naive)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.pbs<-genefilter(data_df[,which(colnames(data_df) %in% selection_pbs)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.cd137<-genefilter(data_df[,which(colnames(data_df) %in% selection_cd137)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.cd40<-genefilter(data_df[,which(colnames(data_df) %in% selection_cd40)], ffun1)


whichFilter <- (whichFilter.combo | whichFilter.naive | whichFilter.pbs | whichFilter.cd137 | whichFilter.cd40 )   # no sale

selection_data_filter<-data_df[whichFilter,]


group=c(rep("CD137", 3), rep("CD40", 3), rep("COMBO",3), rep("NAIVE",3),rep("PBS",3))

d0<-DGEList(selection_data_filter, group=group)


### filter low expression genes
 
#voom transformation and calculation variance weights

snames<-group

design<-model.matrix(~0+snames)

setwd("/media/inmuno/04_ABELLA/imagenes_3/")
#rownames(design) <- colnames(d)
#colnames(design) <- levels(snames)
pdf("voom.pdf")
data_Normalized <- voom(d0, design, plot = T)
dev.off()

k<-data_Normalized$E   #20498 9

#save(k, file="/media/inmuno/04_ABELLA/data_Normalized_ABELLA_NOV21_selected_samples.Rdata")
save(k, file="/media/inmuno/04_ABELLA/data_Normalized_ABELLA_NOV21_selected_samples_without_treatments.Rdata")

#fitting linear models

fit <- lmFit(data_Normalized, design)

#make contrasts
contr <- makeContrasts(                       
                       CombovsPBS= snamesCombo - snamesPBS,
                       PBSvsNaive= snamesPBS - snamesNaive,
                       CombovsNaive = snamesCombo - snamesNaive,
                       levels = colnames(coef(fit)))


#estimate contrast for each gene
contr <- makeContrasts(                       
                       CombovsPBS= snamesCOMBO - snamesPBS,
                       cd137vsPBS= snamesCD137 - snamesPBS,
                       cd40vsPBS = snamesCD40 - snamesPBS,
                       levels = colnames(coef(fit)))



tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

system("mkdir -p /media/inmuno/04_ABELLA/DEgenes_NUEVO/")
setwd("/media/inmuno/04_ABELLA/DEgenes_NUEVO/")

# PBSvsNaive <- topTable(tmp, coef = "PBSvsNaive", n=nrow(tmp), adjust="fdr")
# write.table(PBSvsNaive, file = "PBSvsNaive.txt", quote = FALSE, row.names = TRUE, sep="\t")

CombovsPBS <- topTable(tmp, coef = "CombovsPBS", n=nrow(tmp), adjust="fdr")
write.table(CombovsPBS, file = "CombovsPBS.txt", quote = FALSE, row.names = TRUE, sep="\t")

PBSvsNaive <- topTable(tmp, coef = "PBSvsNaive", n=nrow(tmp), adjust="fdr")
write.table(PBSvsNaive, file = "PBSvsNaive.txt", quote = FALSE, row.names = TRUE, sep="\t")

CombovsNaive <- topTable(tmp, coef = "CombovsNaive", n=nrow(tmp), adjust="fdr")
write.table(CombovsNaive, file = "CombovsNaive.txt", quote = FALSE, row.names = TRUE, sep="\t")

# Combovs40 <- topTable(tmp, coef = "Combovs40", n=nrow(tmp), adjust="fdr")
# write.table(Combovs40, file = "Combovs40.txt", quote = FALSE, row.names = TRUE, sep="\t")

# Combovs137 <- topTable(tmp, coef = "Combovs137", n=nrow(tmp), adjust="fdr")
# write.table(Combovs137, file = "Combovs137.txt", quote = FALSE, row.names = TRUE, sep="\t")

cd40vsPBS<-topTable(tmp, coef = "cd40vsPBS", n=nrow(tmp), adjust="fdr")
# write.table(k40vsPBS, file = "k40vsPBS.txt", quote = FALSE, row.names = TRUE, sep="\t")

cd137vsPBS<-topTable(tmp, coef = "cd137vsPBS", n=nrow(tmp), adjust="fdr")
# write.table(k137vsPBS, file = "k137vsPBS.txt", quote = FALSE, row.names = TRUE, sep="\t")


# QC control

xlim_r = c(min(na.omit(log2(d0$counts+1))),max(na.omit(log2(d0$counts+1)))) 
xlim = c(min(na.omit(data_Normalized$E)),max(na.omit(data_Normalized$E)))
clust.euclid.average <- hclust(dist(t(log2(d0$counts+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(data_Normalized$E)),method="average")

#system("mkdir -p /media/inmuno/04_ABELLA/images/")
setwd("/media/inmuno/04_ABELLA/imagenes_3/")

pdf(file="RNASEQ_ABELLA_QC.pdf", colormode="rgb", width=20, height=20)

boxplot(log2(d0$counts+1),names=colnames(d0$counts), cex.axis=0.7, las=2)
multi("density", log2(d0$counts+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(data_Normalized$E, names=colnames(data_Normalized$E), cex.axis=0.7, las=2)
multi("density", data_Normalized$E, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)

dev.off()



### PCA

library(ggplot2)
x="PC1"
y="PC2"
type<-snames
fit <- prcomp(t(na.omit(as.matrix(data_Normalized$E))), scale=T)

#type<-snames[c(1:12,17:22)]
#fit2 <- prcomp(t(na.omit(as.matrix(data_Normalized$E))), scale=T)

pcData <- data.frame(fit$x)
#pcData<-data.frame(fit2$x)
#ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point( size = 6)  + theme_bw() + theme(legend.title=element_blank())

ggp <- ggplot(pcData, aes_string(x=x, y=y)) + geom_point(aes(colour=type, shape=type), size = 5) + scale_shape_manual(values=c(15,18,17,16)) + theme_bw() + theme(legend.title=element_blank())

#pdf(file = "PCA_RNASeq_MALVAREZ_without_1-3_3-6_samples_AGO21.pdf", width = 15, height = 15, colormodel = "rgb")
pdf(file = "PCA_RNASeq_ABELLA_LABEL_NOV21_labels.pdf", width = 6, height = 6, colormodel = "rgb")

ggp
dev.off()


# heatmap completo


library(pheatmap)

library(colorspace)

tmp <- detOutliers(data_Normalized$E, type, "DetOutliers_RNA_SEQ_ABELLA_NOV21.pdf", 14, 14)



col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(data_Normalized$E)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue","green"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("Combo","Naive","PBS")

colores<-diverge_hsv(15) #blue, white,red
mat_breaks <- seq(min(data_Normalized$E), max(data_Normalized$E), length.out = 10)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(data_Normalized$E), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(data_Normalized$E)/length(colores), max(data_Normalized$E), length.out=floor(length(colores)/2)))



pdf("Heatmap_ALL_GENES_EXPRESSED_ABELLA_NOV21_selected_samples.pdf")

pheatmap(data_Normalized$E,
    color=colores,
    show_rownames     = FALSE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 14,
    main              = "Heatmap of normalized samples \n (all genes expressed)",
    cluster_cols = T,
    #border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("270"),
    #labels_col=c("Control","BO 112","DMX AA","BO 112 + DMX AA"),

)
dev.off()

#heatmap de solo genes diferencialmente expresados pvalue 0.05 & B > 0

k<-(data_Normalized$E)
k<-as.data.frame(k)
k$ENSG<-rownames(k)
k_m<-merge(k, mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)
rownames(k_m)<-k_m$gene_name

k_m_sel<-k_m[which(k_m$gene_name %in% c("Fasl","Gzmb","Ccl5","Ifng","Tigit","Pdcd1","Lag3","Ctla4","Foxp3","Timd4")),]
rownames(k_m_sel)<-k_m_sel$gene_name
k_m_sel$ENSG<-NULL
k_m_sel$gene_name<-NULL
k_m_sel<-k_m_sel[,c(1:9,13:15)]


xt<-t(k_m_sel)
xts<-scale(xt)
k_sel_s<-t(xts)

col_groups<-c(rep("CD137L",3),rep("CD40L",3),rep("COMBO",3),rep("PBS",3))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue","green","orange"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("CD137L", "CD40L" , "COMBO","PBS")


colfunc <- colorRampPalette(c("green", "black","red"))

colores <- colfunc(5)


mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 5)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(data_Normalized_DEG$E), 0, length.out=ceiling(length(colores)/2) + 1), 
#              seq(max(data_Normalized_DEG$E)/length(colores), max(data_Normalized_DEG$E), length.out=floor(length(colores)/2)))

pdf("Heatmap_suplementario_genes_paper.pdf",5,4)

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    #breaks            = mat_breaks,
    annotation_colors = mat_colors,
    border_color=FALSE,
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 10,
    main              = "Heatmap zscored ",
    cluster_cols = F,
    #border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=10, 
    cellwidth = 12,

)
dev.off()





# gene selection

#PBSvsNaive_de<-PBSvsNaive[which(PBSvsNaive$adj.P.Val < 0.05 ),] 
CombovsPBS_de<-CombovsPBS[which(CombovsPBS$adj.P.Val < 0.05),] #345
#Combovs40_de<-Combovs40[which(Combovs40$adj.P.Val < 0.05 ),] 
#Combovs137_de<-Combovs137[which(Combovs137$adj.P.Val < 0.05 ),] 
PBSvsNaive_de<-PBSvsNaive[which(PBSvsNaive$adj.P.Val < 0.05),] # 376
CombovsNaive_de<-CombovsNaive[which(CombovsNaive$adj.P.Val < 0.05),] #1796

cd137vsPBS_de<-cd137vsPBS[which(cd137vsPBS$adj.P.Val < 0.05),] # 376
cd40vsPBS_de<-cd40vsPBS[which(cd40vsPBS$adj.P.Val < 0.05),] #1796

# significativo el cambio

CombovsPBS_de_signif<-CombovsPBS_de[which(CombovsPBS_de$logFC < -1 | CombovsPBS_de$logFC > 1 ),]
cd137vsPBS_de_signif<-cd137vsPBS_de[which(cd137vsPBS_de$logFC < -1 | cd137vsPBS_de$logFC > 1 ),]
cd40vsPBS_de_signif<-cd40vsPBS_de[which(cd40vsPBS_de$logFC < -1 | cd40vsPBS_de$logFC > 1 ),]


pdf("Genes_DE_3_conditions_vs_PBS_significativos.pdf")
compare3List(paste(unique(rownames(CombovsPBS_de_signif))),paste(unique(rownames(cd137vsPBS_de_signif))), paste(unique(rownames(cd40vsPBS_de_signif))),"Combo vs PBS","CD137-L vs PBS", "CD40-L vs PBS","Venn Diagram Genes with adj.P.value < 0.05 and LogFC <-1 | > 1")
dev.off()


myLists<-list("Combo vs PBS"=paste(unique(rownames(CombovsPBS_de_signif))),"CD137L vs PBS"=paste(unique(rownames(cd137vsPBS_de_signif))), "CD40L vs PBS"=paste(unique(rownames(cd40vsPBS_de_signif))))

Vstem <- Venn(myLists)

pdf("Venn_weighted_ABELLA.pdf")
plot(Vstem, doWeights = TRUE)
dev.off()

########### hacer enriquecimientos de los genes 42 7 50




#guardar para GSEA
library(xlsx)
write.xlsx(CombovsPBS_de, file ="CombovsPBS_de.xlsx")
write.xlsx(PBSvsNaive_de, file ="PBSvsNaive_de.xlsx")
write.xlsx(CombovsNaive_de, file ="CombovsNaive_de.xlsx")


PBSvsNaive_heatmap_matrix<-data_Normalized[which(rownames(data_Normalized) %in% rownames(PBSvsNaive_de)),]

#Genes del COMBO DEG que no están en las otras dos condiciones
pdf("Genes_DE_3_conditions_vs_PBS.pdf")
compare3List(paste(unique(rownames(CombovsPBS_de))),paste(unique(rownames(k137vsPBS_de))), paste(unique(rownames(k40vsPBS_de))),"Combo vs PBS","CD137-L vs PBS", "CD40-L vs PBS","Venn Diagram Genes with P.value < 0.05")
dev.off()


#### ggvenn
library(ggVennDiagram)

x<-list("Combo vs PBS"= rownames(CombovsPBS_de_signif), "CD137-L vs PBS"=rownames(cd137vsPBS_de_signif), "CD40-L vs PBS"=rownames(cd40vsPBS_de_signif))

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_en_las_3_comparativas_vs_PBS_significativos.pdf", height=10, width=10)
plot(ggVennDiagram(x) + scale_fill_gradient(low="white",high = "white") + ggtitle("DEgenes with adj.p.value < 0.05 % and logFC <-1 | > 1")+ theme(legend.position="none"))
plot(ggvenn(x, fill_color = c("#ff0000", "#008000", "#0000ff"))  + ggtitle("DEgenes with adj.p.value < 0.05 % and logFC <-1 | > 1")+ theme(legend.position="none"))

dev.off()



x<-list("Combo vs PBS"= rownames(CombovsPBS_de[which(CombovsPBS_de$logFC >1 | CombovsPBS_de$logFC < -1 ),]), "PBS vs Naive"=rownames(PBSvsNaive_de[which(PBSvsNaive_de$logFC >1 | PBSvsNaive_de$logFC < -1 ),]), "Combo vs Naive"=rownames(CombovsNaive_de[which(CombovsNaive_de$logFC >1 | CombovsNaive_de$logFC < -1 ),]))

pdf("Venn_Diagramm_Genes_diferencialmente_expresados_en_las_3_comparativas_con_logFC significativo.pdf", height=10, width=10)
plot(ggVennDiagram(x) + scale_fill_gradient(low="white",high = "white") + ggtitle("DEgenes with p.value < 0.05 % & logFC >1 | logFC < -1")+ theme(legend.position="none"))
dev.off()

length(x)

CombovsPBS_de_s<-CombovsPBS_de[which(CombovsPBS_de$logFC < -1 | CombovsPBS_de$logFC > 1),]   # 245  DEG y significativos. 

#Comprobar en las otras dos condiciones si tienen un valor mayor a 0.1  # No lo son

# selection<-rownames(CombovsPBS_de_s)

# ls<-PBSvsNaive_de[which(rownames(PBSvsNaive_de) %in% selection),]
# ll<-CombovsNaive_de[which(rownames(CombovsNaive_de) %in% selection),]
data_Normalized_DEG<-data_Normalized[which(rownames(data_Normalized) %in% rownames(CombovsPBS_de_s)),]

# heatmap de 


head(data_Normalized_DEG)


col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(data_Normalized_DEG$E)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue","green"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("Combo","PBS","Naive")

#colores<-diverge_hsv(6) #blue, white,red
colores <- wes_palette("Zissou1", 100, type = "continuous")
mat_breaks <- seq(min(data_Normalized_DEG$E), max(data_Normalized_DEG$E), length.out = 8)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(data_Normalized_DEG$E), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(data_Normalized_DEG$E)/length(colores), max(data_Normalized_DEG$E), length.out=floor(length(colores)/2)))

pdf("Heatmap_of_normalized_samples_GENES_DEG_COMBOvsPBS_pval0.05_logFC_significativo_NOV21.pdf")

pheatmap(data_Normalized_DEG,
    color=colores,
    show_rownames     = FALSE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 10,
    main              = "Heatmap of normalized samples \n GENES DEG COMBO vs PBS \n pval < 0.05 and logFC < -1 | logFC > 1 ",
    cluster_cols = T,
    #border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()

pdf("Heatmap_GENES_DEG_COMBOvsPBS_pval0.05_logFC_significativo_forzando_clustering_grupos_NOV21.pdf")

pheatmap(data_Normalized_DEG,
    color=colores,
    show_rownames     = FALSE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 10,
    main              = "Heatmap of normalized samples \n GENES DEG COMBO vs PBS \n pval < 0.05 and logFC < -1 | logFC > 1 ",
    cluster_cols = F,
    #border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()




 `%!in%` <- Negate(`%in%`)

selection1<- CombovsPBS_de_signif[which(rownames(CombovsPBS_de_signif) %in%  rownames(cd40vsPBS_de_signif)),]
selection1_cd40<-selection1[which(rownames(selection1) %!in% rownames(cd137vsPBS_de_signif)),]

selection2<- CombovsPBS_de_signif[which(rownames(CombovsPBS_de_signif) %in%  rownames(cd137vsPBS_de_signif)),]
selection2_cd137<-selection2[which(rownames(selection2) %!in% rownames(cd40vsPBS_de_signif)),]


##### heatmaps de estos genes 

k<-as.data.frame(data_Normalized$E)

cd40_leading<-k[which(rownames(k) %in% rownames(selection1_cd40)),]
cd40_leading_sel<-cd40_leading[,c(1:9,13:15)]
cd40_leading_sel$ENSG<-rownames(cd40_leading_sel)


cd137_leading<-k[which(rownames(k) %in% rownames(selection2_cd137)),]
cd137_leading_sel<-cd137_leading[,c(1:9,13:15)]
cd137_leading_sel$ENSG<-rownames(cd137_leading_sel)

cd137_leading_sel$ENSG<-paste(lapply(strsplit(paste(cd137_leading_sel$ENSG),"\\."),"[",1))
cd40_leading_sel$ENSG<-paste(lapply(strsplit(paste(cd40_leading_sel$ENSG),"\\."),"[",1))


### hacer enriquecimiento de los genes selection1_cd40 y selection2_cd137

library(fgsea)
library(clusterProfiler)
library("org.Mm.eg.db")

Combo_PBS_137<-CombovsPBS[which(rownames(CombovsPBS) %in% rownames(cd137_leading_sel)),]
Combo_PBS_40<-CombovsPBS[which(rownames(CombovsPBS) %in% rownames(cd40_leading_sel)),]


Combo_PBS_137$ENSG<-paste(lapply(strsplit(paste(rownames(Combo_PBS_137)),"\\."),"[",1))
Combo_PBS_40$ENSG<-paste(lapply(strsplit(paste(rownames(Combo_PBS_40)),"\\."),"[",1))

ids <- bitr(Combo_PBS_137$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

Combo_PBS_137_entrez<-merge(Combo_PBS_137, ids, by.x="ENSG", by.y="ENSEMBL")
Combo_PBS_137_entrez<-unique(Combo_PBS_137_entrez)

ids <- bitr(Combo_PBS_40$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

Combo_PBS_40_entrez<-merge(Combo_PBS_40, ids, by.x="ENSG", by.y="ENSEMBL")
Combo_PBS_40_entrez<-unique(Combo_PBS_40_entrez)



cd137_leading_sel_entrez_selUP<-na.omit(unique(Combo_PBS_137_entrez[which( Combo_PBS_137_entrez$logFC >1), "ENTREZID"]))

enriched_GOBP_cd137_leading_sel_entrez_selUP <- enrichGO(cd137_leading_sel_entrez_selUP, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_cd137_leading_sel_entrez_selUP<-as.data.frame(enriched_GOBP_cd137_leading_sel_entrez_selUP)

#enriched_GOBP_IL12IPvsIL12IV_df<-as.data.frame(enriched_GOBP_IL12IPvsIL12IV)
#enriched_GOBP_IL12IPvsIL12IV_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_IL12IPvsIL12IV$p.adjust)))


enriched_GOBP_cd137_leading_sel_entrez_selDOWN<-na.omit(unique(Combo_PBS_137_entrez[which( Combo_PBS_137_entrez$logFC < (-1)), "ENTREZID"]))

enriched_GOBP_cd137_leading_sel_entrez_selDOWN <- enrichGO(enriched_GOBP_cd137_leading_sel_entrez_selDOWN, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_cd137_leading_sel_entrez_selDOWN<-as.data.frame(enriched_GOBP_cd137_leading_sel_entrez_selDOWN)
#enriched_GOBP_IL12IPvsIL12IV_down_s <- simplify(enriched_GOBP_IL12IPvsIL12IV_down, cutoff=0.7, by="p.adjust", select_fun=min)   # DOWN


selup<-mutate(enriched_GOBP_cd137_leading_sel_entrez_selUP, logp =(-1)*log(p.adjust, base=10))
seldown<-mutate(enriched_GOBP_cd137_leading_sel_entrez_selDOWN, logp =(1)*log(p.adjust, base=10))

fgseaRes_sel =rbind(selup,seldown)


fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$logp > 0, "Up-regulated", "Down-regulated")
filtRes = fgseaRes_sel
filtRes<-na.omit(filtRes)

upcols =  colorRampPalette(colors = c("red4", "red4", "red4"))( sum(filtRes$Enrichment == "Up-regulated"))
downcols =  colorRampPalette(colors = c( "blue4", "blue4", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
downcols_rev<-rev(downcols)
colos = c(upcols, downcols_rev)
names(colos) = 1:length(colos)

filtRes$Index = as.factor(c(1:length(upcols),rev(sum(length(upcols),length(downcols)-1):length(upcols)+1)))
title_plot<-"GOBP Combo vs PBS leading CD137L"

g_fgseaRes = ggplot(filtRes, aes(reorder(Description, logp), logp)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="logP",
        title=title_plot) + 
    theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 10)) 


pdf("enrichment_GOBP_COMBOvsPBS_leading_CD137.pdf",10,5)
plot(g_fgseaRes)
dev.off()

cd40_leading_sel_entrez_selUP<-na.omit(unique(Combo_PBS_40_entrez[which( Combo_PBS_40_entrez$logFC >1), "ENTREZID"]))

enriched_GOBP_cd40_leading_sel_entrez_selUP <- enrichGO(cd40_leading_sel_entrez_selUP, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_cd40_leading_sel_entrez_selUP <- simplify(enriched_GOBP_cd40_leading_sel_entrez_selUP, cutoff=0.6, by="p.adjust", select_fun=min)   # DOWN

enriched_GOBP_cd40_leading_sel_entrez_selUP<-as.data.frame(enriched_GOBP_cd40_leading_sel_entrez_selUP)
#enriched_GOBP_IL12IPvsIL12IV_df<-as.data.frame(enriched_GOBP_IL12IPvsIL12IV)
#enriched_GOBP_IL12IPvsIL12IV_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_IL12IPvsIL12IV$p.adjust)))


enriched_GOBP_cd40_leading_sel_entrez_selDOWN<-na.omit(unique(Combo_PBS_40_entrez[which( Combo_PBS_40_entrez$logFC < (-1)), "ENTREZID"]))

enriched_GOBP_cd40_leading_sel_entrez_selDOWN <- enrichGO(enriched_GOBP_cd40_leading_sel_entrez_selDOWN, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_cd40_leading_sel_entrez_selDOWN<-as.data.frame(enriched_GOBP_cd40_leading_sel_entrez_selDOWN)
#enriched_GOBP_IL12IPvsIL12IV_down_s <- simplify(enriched_GOBP_IL12IPvsIL12IV_down, cutoff=0.7, by="p.adjust", select_fun=min)   # DOWN


selup<-mutate(enriched_GOBP_cd40_leading_sel_entrez_selUP, logp =(-1)*log(p.adjust, base=10))
seldown<-mutate(enriched_GOBP_cd40_leading_sel_entrez_selDOWN, logp =(1)*log(p.adjust, base=10))

fgseaRes_sel =rbind(selup,seldown)


fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$logp > 0, "Up-regulated", "Down-regulated")
filtRes = fgseaRes_sel
filtRes<-na.omit(filtRes)

upcols =  colorRampPalette(colors = c("red4", "red4", "red4"))( sum(filtRes$Enrichment == "Up-regulated"))
#downcols =  colorRampPalette(colors = c( "blue4", "blue4", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
#downcols_rev<-rev(downcols)
colos = c(upcols)
names(colos) = 1:length(colos)

filtRes$Index = as.factor(1:length(upcols))
title_plot<-"GOBP Combo vs PBS leading cd40L"

g_fgseaRes = ggplot(filtRes, aes(reorder(Description, logp), logp)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="logP",
        title=title_plot) + 
    theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 10)) 


pdf("enrichment_GOBP_COMBOvsPBS_leading_cd40.pdf",10,5)
plot(g_fgseaRes)
dev.off()






library(rtracklayer)
setwd("/media/inmuno/04_ABELLA/DEgenes_NUEVO/")

mouse_gtf<-import("/media/inmuno/references/STAR_mouse_mm39/gencode.vM27.annotation.gtf")

mouse_gtf<-as.data.frame(mouse_gtf) 
mouse_gtf_sel<-mouse_gtf[,c("gene_id","gene_name")]
mouse_gtf_sel<-unique(mouse_gtf_sel)

cd40_leading_selec<-merge(cd40_leading_sel, mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)
cd40_leading_selec$ENSG<-NULL
rownames(cd40_leading_selec)<-cd40_leading_selec$gene_name
cd40_leading_selec$gene_name<-NULL

cd137_leading_selec<-merge(cd137_leading_sel, mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)
cd137_leading_selec$ENSG<-NULL
rownames(cd137_leading_selec)<-cd137_leading_selec$gene_name
cd137_leading_selec$gene_name<-NULL




xt<-t(cd40_leading_selec)
xts<-scale(xt)
k_sel_s<-t(xts)

xt<-t(cd137_leading_selec)
xts<-scale(xt)
k_sel_s<-t(xts)



col_groups<-c(rep("CD137L",3),rep("CD40L",3),rep("COMBO",3),rep("PBS",3))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue","green","orange"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("CD137L", "CD40L" , "COMBO","PBS")


colores<-diverge_hsv(3) #blue, white,red

pdf("Heatmap_COMBOvsPBS_pval0.05_logFC_significativo_leadingCD40L_ENE22v4.pdf",6,6)

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
   # breaks            = myBreaks,
    annotation_colors = mat_colors,
    border_color="black",
    #gaps_col=3,
    drop_levels       = FALSE,
    fontsize          = 8,
    main              = "GENES DEG COMBO vs PBS leading by CD40L \n pval < 0.05 and logFC < -1 | logFC > 1 ",
    cluster_cols = F,
    #border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=8, 
    cellwidth = 10,

)
dev.off()

colores<-diverge_hsv(4) #blue, white,red

pdf("Heatmap_COMBOvsPBS_pval0.05_logFC_significativo_leadingCD137L_ENE22v3.pdf",5,7)

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
   # breaks            = myBreaks,
    annotation_colors = mat_colors,
    border_color="black",
    #gaps_col=3,
    drop_levels       = FALSE,
    fontsize          = 8,
    main              = "GENES DEG COMBO vs PBS leading by CD137L \n pval < 0.05 and logFC < -1 | logFC > 1 ",
    cluster_cols = F,
    #border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=8, 
    cellwidth = 10,

)
dev.off()



# list_sel_tmp<-rownames(CombovsPBS_de)[selection]

# Combo_diff<-setdiff(list_sel_tmp, rownames(k40vsPBS_de))

# anotar Combo_diff
library("org.Mm.eg.db")
library(rtracklayer)
setwd("/media/inmuno/04_ABELLA/DEgenes_NUEVO/")

mouse_gtf<-import("/media/inmuno/references/STAR_mouse_mm39/gencode.vM27.annotation.gtf")

mouse_gtf<-as.data.frame(mouse_gtf) 
mouse_gtf_sel<-mouse_gtf[,c("gene_id","gene_name")]
mouse_gtf_sel<-unique(mouse_gtf_sel)

mouse_gtf_selection<-mouse_gtf_sel[which(mouse_gtf_sel$gene_id %in% rownames(data_Normalized_DEG$E)),]

mouse_gtf_selection<-unique(mouse_gtf_selection)

k<-as.data.frame(data_Normalized_DEG$E)
k$ENSG<-rownames(k)

CombovsPBS_de_selection_s<-merge(k, mouse_gtf_selection, by.x="ENSG", by.y="gene_id", all.x=T)

write.xlsx(CombovsPBS_de_selection_s, file="Combo_vs_PBS_metrics_anotada.xlsx")


########## Anotación de genes
#PBSvsNaive

library(rtracklayer)
library(xlsx)

mouse_gtf<-import("/home/jgonzalezgom/references/STAR_mouse_GHRC39/gencode.vM27.annotation.gtf")

gtf_selection<-data.frame("ENSG"=mouse_gtf$gene_id, "SYMBOL"=mouse_gtf$gene_name)
gtf_selection<-unique(gtf_selection)

PBSvsNaive$ENSG<-rownames(PBSvsNaive)
PBSvsNaive_annot<-merge(PBSvsNaive, gtf_selection, by="ENSG", all.x=T)
PBSvsNaive_annot_de<-PBSvsNaive_annot[which(PBSvsNaive_annot$adj.P.Val < 0.05),]
write.xlsx(PBSvsNaive_annot_de, file="PBSvsNaive_annot_de.xlsx")

CombovsPBS$ENSG<-rownames(CombovsPBS)
CombovsPBS_annot<-merge(CombovsPBS, gtf_selection, by="ENSG", all.x=T)
CombovsPBS_annot_de<-CombovsPBS_annot[which(CombovsPBS_annot$adj.P.Val < 0.05),]
write.xlsx(CombovsPBS_annot_de, file="CombovsPBS_annot_de.xlsx")

CombovsNaive$ENSG<-rownames(CombovsNaive)
CombovsNaive_annot<-merge(CombovsNaive, gtf_selection, by="ENSG", all.x=T)
CombovsNaive_annot_de<-CombovsNaive_annot[which(CombovsNaive$adj.P.Val < 0.05),]
write.xlsx(CombovsNaive_annot_de, file="CombovsNaive_annot_de.xlsx")

library(ggrepel)

setwd("/media/inmuno/04_ABELLA/imagenes_3/")

volcanoPlots_pepe<-function(dataset,pvalth, FCth,ntop ,title){

    dataset<-dataset[,1:8]
    dataset<-unique(dataset)

    dataset <- dataset %>% 
    mutate(
        Expression = case_when(logFC >= FCth & adj.P.Val < pvalth ~ "Up-regulated",
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



    top <- ntop
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

PBSvsNaive_plot<-volcanoPlots_pepe(PBSvsNaive_annot,0.05,1,0, "PBS vs Naive")
CombovsPBS_plot<-volcanoPlots_pepe(CombovsPBS_annot,0.05,1,0, "Combo vs PBS")
CombovsNaive_plot<-volcanoPlots_pepe(CombovsNaive_annot,0.05,1,0, "Combo vs Naive")


tiff("VOLCANO_PLOTS_Condiciones_nuevas_pval0.05_without_labels_Combo_vs_Naive")
#plot(PBSvsNaive_plot)
#plot(CombovsPBS_plot)
plot(CombovsNaive_plot)

dev.off()

#data_Normalized_sel<-data_Normalized[which(rownames(data_Normalized) %in% CombovsPBS_de_selection_s$ENSG ),]

####### Heatmap de genes seleccionados por FERNANDO


lista<-read.csv("/media/inmuno/04_ABELLA/lista_fer_2.txt", fileEncoding="UTF-16LE", sep="\t", header=F)


from #gtf_selection

ensgs<-paste(lapply(strsplit(paste(gtf_selection$ENSG),"\\."),"[",1))

gtf_selection2<-cbind(gtf_selection,ensgs)

lista_merged<-merge(lista, gtf_selection2, by.x="V2", by.y="ensgs", all.x=T)

# selection_from data_Normalized 

k_sel<-data_Normalized[which(rownames(data_Normalized) %in% lista_merged$ENSG),]
k_sel<-k_sel$E
k_sel<-as.data.frame(k_sel)

lista_merged_selection<-lista_merged[which(lista_merged$ENSG %in% rownames(k_sel)),]
lista_merged_selection<-unique(lista_merged_selection)
lista_merged_selection<-lista_merged_selection[,c("ENSG","SYMBOL")]

k_sel$ENSG<-rownames(k_sel)

k_sel<-merge(k_sel, lista_merged_selection, by="ENSG", all.x=T)

rownames(k_sel)<-k_sel$SYMBOL

`%!in%` <- Negate(`%in%`)
unselect<-c("ENSG","SYMBOL")

k_sel_s<-k_sel[which(colnames(k_sel) %!in% unselect)]

# scale rows

xt<-t(data_o)
xts<-scale(xt)
k_sel_s<-t(xts)

col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue","green"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("Combo", "Naïve" , "PBS")

colores<-diverge_hsv(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_ABELLA_pval0.05_logFC1_zscore.pdf")

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = FALSE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Zscore GENES DEG COMBO vs PBS \n pval < 0.05 and logFC < −1 | logFC > 1 ",
    cluster_cols = F,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()



CombovsPBS_de_selection_s_signif<-CombovsPBS_de_selection_s[which(CombovsPBS_de_selection_s$logFC < -1 | CombovsPBS_de_selection_s$logFC > 1 ),]


data_Normalized_sel_signif<-data_Normalized_sel[which(rownames(data_Normalized_sel) %in% CombovsPBS_de_selection_s_signif$ENSG),]


ensgs<-rownames(data_Normalized_sel_signif)

df_traza<-data.frame("ENSG"=ensgs)

df_traza<-merge(df_traza, gtf_selection, by="ENSG", all.x=T)

dataN_sel_signif_sort<-data_Normalized_sel_signif[order(rownames(data_Normalized_sel_signif)),]

df_traza_sel<-df_traza[order(df_traza$ENSG),]

rownames(data_Normalized_sel_signif)<-df_traza_sel$SYMBOL

pdf("Heatmap_of_normalized_samples_Combo_not_in_others_VS_PBS_DEG_pval0.05_signif_NOV21.pdf")

pheatmap(data_Normalized_sel_signif,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    annotation_colors = mat_colors,
    drop_levels       = TRUE,
    fontsize          = 6,
    main              = "Heatmap of normalized samples \n DEgenes in Combo not present in others treatments \n pval < 0.05 and LogFC < -1 or LogFC >1",
    cluster_cols = F,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()



#### GSEA ANALISIS DE LA LISTA DE 88 genes

library(fgsea)
library(gage)

# ENSg to eNTREZID
#CombovsPBS_annot<-CombovsPBS_annot[which(CombovsPBS_annot$adj.P.Val < 0.05),] #345
#Combovs40_de<-Combovs40[which(Combovs40$adj.P.Val < 0.05 ),] 
#Combovs137_de<-Combovs137[which(Combovs137$adj.P.Val < 0.05 ),] 
#PBSvsNaive_annot<-PBSvsNaive_annot[which(PBSvsNaive_annot$adj.P.Val < 0.05),] # 376
#CombovsNaive_annot<-CombovsNaive_annot[which(CombovsNaive_annot$adj.P.Val < 0.05),] #1796

CombovsNaive_annot_order<-CombovsNaive_annot[order(CombovsNaive_annot$t, decreasing=T),]
CombovsPBS_annot_order<-CombovsPBS_annot[order(CombovsPBS_annot$t, decreasing=T),]
PBSvsNaive_annot_order<-PBSvsNaive_annot[order(PBSvsNaive_annot$t, decreasing=T),]

# ordenar por el estadistico 

preparedata_GSEA<-function(dataframe){
    listed_rank<-list()
    for(i in 1:nrow(dataframe)){
        symbol<-dataframe[,"SYMBOL"][i]
        symbol<-toupper(symbol)
        t_val<-dataframe[,"t"][i]
        listed_rank[[symbol]]<-t_val
    }
    listed_rank2<-unlist(listed_rank)
    return(listed_rank2)
}

CombovsNaive_prepared<-preparedata_GSEA(CombovsNaive_annot_order)
CombovsPBS_prepared<-preparedata_GSEA(CombovsPBS_annot_order)
PBSvsNaive_prepared<-preparedata_GSEA(PBSvsNaive_annot_order)


pathways <- gmtPathways("/media/inmuno/04_ABELLA/DEgenes_NUEVO/c7.immunesigdb.v7.4.symbols.gmt")
pathways<- gmtPathways("/media/inmuno/04_ABELLA/DEgenes_NUEVO/c5.go.bp.v7.4.symbols.gmt")

#examplePathways<- pathways



fgseaResCombovsNaive <- fgsea(pathways = pathways, stats    = CombovsNaive_prepared , maxSize=500)
#collapsedPathwaysfgseaResCombovsNaive <- collapsePathways(fgseaResCombovsNaive[order(pval)][padj < 0.05], pathways,CombovsNaive_prepared )
#mainPathwaysfgseaResCombovsNaive <- fgseaResCombovsNaive[pathway %in% collapsedPathwaysfgseaResCombovsNaive$mainPathways][order(-NES), pathway]

fgseaResCombovsPBS <- fgsea(pathways = pathways, stats    = CombovsPBS_prepared , maxSize=500)

fgseaResPBSvsNaive <- fgsea(pathways = pathways, stats    = PBSvsNaive_prepared , maxSize=500)



topPathwaysUp_CombovsNaive <- fgseaResCombovsNaive[ES > 0][head(order(pval), n=30), pathway]
topPathwaysDown_CombovsNaive <- fgseaResCombovsNaive[ES < 0][head(order(pval), n=30), pathway]
#topPathways_CombovsNaive <- c(topPathwaysUp_CombovsNaive, rev(topPathwaysDown_CombovsNaive))

topPathwaysUp_CombovsPBS <- fgseaResCombovsPBS[ES > 0][head(order(pval), n=30), pathway]
topPathwaysDown_CombovsPBS <- fgseaResCombovsPBS[ES < 0][head(order(pval), n=30), pathway]


topPathwaysUp_PBSvsNaive <- fgseaResPBSvsNaive[ES > 0][head(order(pval), n=30), pathway]
topPathwaysDown_PBSvsNaive <- fgseaResPBSvsNaive[ES < 0][head(order(pval), n=30), pathway]




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

COMBO_vs_NAIVE_plot<-plot_GSEA_PEPE(fgseaResCombovsNaive, CombovsNaive_prepared, topPathwaysUp_CombovsNaive, topPathwaysDown_CombovsNaive, pathways, 1, 500, 10, "GSEA c5.GeneOntology BP Combo vs Naive" )
COMBO_vs_PBS_plot<-plot_GSEA_PEPE(fgseaResCombovsPBS, CombovsPBS_prepared, topPathwaysUp_CombovsPBS, topPathwaysDown_CombovsPBS, pathways, 1, 500, 10, "GSEA c5.GeneOntology BP Combo vs PBS" )
PBS_vs_NAIVE_plot<-plot_GSEA_PEPE(fgseaResPBSvsNaive, PBSvsNaive_prepared, topPathwaysUp_PBSvsNaive, topPathwaysDown_PBSvsNaive, pathways, 1, 500, 10, "GSEA c5.GeneOntology BP PBS vs Naive" )


pdf("GSEA_ALL_c5_GeneOntology.pdf",30,30)
plot(COMBO_vs_NAIVE_plot)
plot(COMBO_vs_PBS_plot)
plot(PBS_vs_NAIVE_plot)
dev.off()


pdf("fgsea_CombovsNaive_table_top10_up_and_down.pdf",30,30)
plotGseaTable(pathways[topPathways_CombovsNaive], CombovsNaive_prepared ,fgseaResCombovsNaive , gseaParam=0.5)
dev.off()


pdf("fgsea_CombovsPBS_TCELL_vs_TREG_ACT_UP_pathway.pdf")
plotEnrichment(pathways[["GSE7460_CD8_TCELL_VS_TREG_ACT_UP"]], CombovsPBS_prepared) + labs(title="GSE7460_CD8_TCELL_VS_TREG_ACT_UP")
dev.off()

pdf("fgsea_GOBP_ADAPTATIVE_IMMUNE_RESPONSE.pdf")
plotEnrichment(pathways[["GOBP_ADAPTIVE_IMMUNE_RESPONSE"]], CombovsPBS_prepared) + labs(title="GOBP Adaptive Immune response")
dev.off()

pdf("fgsea_GOBP_TGFB_SIGNALING_PATHWAY.pdf")
plotEnrichment(pathways[["GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY"]], CombovsPBS_prepared) + labs(title="GOBP TGFbeta receptor signaling pathway")
dev.off()

pdf("fgsea_GOBP_VASCULATURE_DEVELOPMENT.pdf")
plotEnrichment(pathways[["GOBP_POSITIVE_REGULATION_OF_VASCULATURE_DEVELOPMENT"]], CombovsPBS_prepared) + labs(title="GOBP Positive regulation of vasculature development")
dev.off()


fgseaResPBSvsNaive <- fgsea(pathways = pathways, stats    = PBSvsNaive_prepared , maxSize=500)
collapsedPathwaysfgseaResPBSvsNaive <- collapsePathways(fgseaResPBSvsNaive[order(pval)][padj < 0.05], pathways,PBSvsNaive_prepared )
mainPathwaysfgseaResPBSvsNaive <- fgseaResPBSvsNaive[pathway %in% collapsedPathwaysfgseaResPBSvsNaive$mainPathways][order(-NES), pathway]


topPathwaysUp_PBSvsNaive <- fgseaResPBSvsNaive[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_PBSvsNaive <- fgseaResPBSvsNaive[ES < 0][head(order(pval), n=10), pathway]
topPathways_PBSvsNaive <- c(topPathwaysUp_PBSvsNaive, rev(topPathwaysDown_PBSvsNaive))

# filter FGSEA by using gage results



pdf("fgsea_PBSvsNaive_table_top10_up_and_down.pdf",30,30)
plotGseaTable(pathways[topPathways_PBSvsNaive], PBSvsNaive_prepared ,fgseaResPBSvsNaive , gseaParam=0.5)
dev.off()

pdf("fgsea_PBSvsNaive_first_topUP_pathway.pdf")
plotEnrichment(pathways[["GSE9601_NFKB_INHIBITOR_VS_PI3K_INHIBITOR_TREATED_HCMV_INF_MONOCYTE_DN"]], PBSvsNaive_prepared) + labs(title="Pro-inflamation in monocytes.")
dev.off()


fgseaResCombovsPBS <- fgsea(pathways = pathways, stats    = CombovsPBS_prepared , maxSize=500)
collapsedPathwaysfgseaResCombovsPBS <- collapsePathways(fgseaResCombovsPBS[order(pval)][padj < 0.05], pathways,CombovsPBS_prepared )
mainPathwaysfgseaResCombovsPBS <- fgseaResCombovsPBS[pathway %in% collapsedPathwaysfgseaResCombovsPBS$mainPathways][order(-NES), pathway]


topPathwaysUp_CombovsPBS <- fgseaResCombovsPBS[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown_CombovsPBS <- fgseaResCombovsPBS[ES < 0][head(order(pval), n=10), pathway]
topPathways_CombovsPBS <- c(topPathwaysUp_CombovsPBS, rev(topPathwaysDown_CombovsPBS))

pdf("fgsea_CombovsPBS_table_top10_up_and_down.pdf",30,30)
plotGseaTable(pathways[topPathways_CombovsPBS], CombovsPBS_prepared ,fgseaResCombovsPBS , gseaParam=0.5)
dev.off()

pdf("fgsea_CombovsPBS_first_topUP_pathway.pdf")
plotEnrichment(pathways[["GSE23568_CTRL_VS_ID3_TRANSDUCED_CD8_TCELL_DN"]], CombovsPBS_prepared) + labs(title="Mouse CD8+ T cells affected by ID3")
dev.off()




# fgseaResCombovsNaive_df<-as.data.frame(fgseaResCombovsNaive)
# fgseaResCombovsNaive_pval_less_0.01<-fgseaResCombovsNaive_df[which(fgseaResCombovsNaive_df$pval < 0.01),]


