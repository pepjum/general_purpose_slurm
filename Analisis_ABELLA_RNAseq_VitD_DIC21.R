library(dplyr)
library(edgeR)
library("genefilter")
library(RColorBrewer)
library(viridis)

source("/home/jgonzalezgom/01_NMELERO/funcionesVikv2.R")


files<-list.files("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/counts/", pattern=".txt$")
files<-paste0("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/counts/", files)


samples<-c()
data_df<-data.frame()
for(i in 1:length(files)){
    cat(paste(files[i]),"\n")
    fileloaded<-read.table(files[i], header=T)
    sample_tmp<-paste(lapply(strsplit(paste(files[i]),"\\/"),"[",7))
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
system("mkdir -p /home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/dataframes/")
save(data_df, file="/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/dataframes/RNASEQ_ABELLA_unordered_raw.Rdata")


samples_VitD<-which(grepl("L_", samples))
samples_control<-which(grepl("V_", samples))


data_df_ordered<-data_df[,c(samples_VitD,samples_control)]
save(data_df_ordered, file="/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/dataframes/RNASEQ_ABELLA_ordered_raw.Rdata")

#not in
`%!in%` <- Negate(`%in%`)

not_sel<-c("R3V_S8","R4L_S12")

data_df_selected<-data_df_ordered[,colnames(data_df_ordered) %!in% not_sel ]

selection_VitD<-colnames(data_df_selected)[1:3]
selection_control<-colnames(data_df_selected)[4:6]

data_df_ordered<-data_df_selected

log_intensity_threshold<-5

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.VitD<-genefilter(data_df_ordered[,which(colnames(data_df_ordered) %in% selection_VitD)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.control<-genefilter(data_df_ordered[,which(colnames(data_df_ordered) %in% selection_control)], ffun1)

whichFilter <- (whichFilter.VitD | whichFilter.control )   # no sale

selection_data_filter<-data_df_ordered[whichFilter,]


group=c(rep("VitD", 3), rep("Control", 3))

d0<-DGEList(selection_data_filter, group=group)


### filter low expression genes
 
#voom transformation and calculation variance weights

snames<-group

design<-model.matrix(~0+snames)

setwd("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/imagenes/")
#rownames(design) <- colnames(d)
#colnames(design) <- levels(snames)

pdf("voom2.pdf")
data_Normalized <- voom(d0, design, plot = T)
dev.off()

k<-data_Normalized$E   #19440 6

#save(k, file="/media/inmuno/04_ABELLA/data_Normalized_ABELLA_NOV21_selected_samples.Rdata")
save(k, file="/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/data_Normalized_ABELLA_DIC21_selected_samples.Rdata")

xlim_r = c(min(na.omit(log2(d0$counts+1))),max(na.omit(log2(d0$counts+1)))) 
xlim = c(min(na.omit(data_Normalized$E)),max(na.omit(data_Normalized$E)))
clust.euclid.average <- hclust(dist(t(log2(d0$counts+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(data_Normalized$E)),method="average")


setwd("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/imagenes/")

pdf(file="RNASEQ_ABELLA_VITD_QC2.pdf", colormode="rgb", width=20, height=20)

boxplot(log2(d0$counts+1),names=colnames(d0$counts), cex.axis=0.7, las=2)
multi("density", log2(d0$counts+1), xlim_r, "Density","","")
plot(clust.euclid.average, main="Hierarchical clustering of samples",  hang=-1)
boxplot(data_Normalized$E, names=colnames(data_Normalized$E), cex.axis=0.7, las=2)
multi("density", data_Normalized$E, xlim, "Density","","")
plot(clust.euclid.average2, main="Hierarchical clustering of normalized samples",  hang=-1)

dev.off()


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
pdf(file = "PCA_RNASeq_ABELLA_VitD_DIC21.pdf", width = 6, height = 6, colormodel = "rgb")

ggp
dev.off()




fit <- lmFit(data_Normalized, design)

#make contrasts
contr <- makeContrasts(                       
                       VitDvsControl= snamesVitD - snamesControl,
                       levels = colnames(coef(fit)))


#estimate contrast for each gene

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

system("mkdir -p /home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/DEgenes/")
setwd("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/DEgenes/")


VitDvsControl <- topTable(tmp, coef = "VitDvsControl", n=nrow(tmp), adjust="fdr")
write.table(VitDvsControl, file = "VitDvsControl.txt", quote = FALSE, row.names = TRUE, sep="\t")


library(pheatmap)
library(colorspace)
setwd("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/imagenes/")

tmp <- detOutliers(data_Normalized$E, type, "DetOutliers_RNA_SEQ_ABELLA_DIC21.pdf", 14, 14)


Vitdvscontrol_sel<-VitDvsControl[which(VitDvsControl$P.Value < 0.05),]
Vitdvscontrol_sel_signif<-Vitdvscontrol_sel[which(Vitdvscontrol_sel$logFC < -1 | Vitdvscontrol_sel$logFC > 1 ),]

selection_genes<-data_Normalized$E[which(rownames(data_Normalized$E) %in% rownames(Vitdvscontrol_sel_signif)),]

xt<-t(selection_genes)
xts<-scale(xt)
k_sel_s<-t(xts)

### esto para hacer un heatmap de algunos genes
### load mouse gtf para meter symbols






col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("VitD", "Control")

colores<-diverge_hsv(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_ABELLA_pval0.05_logFC1_zscore2.pdf")

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
    main              = "Zscore GENES DEG VitD vs Control \n pval < 0.05 and signif logFC  ",
    cluster_cols = TRUE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


#selection_genes<-data_Normalized$E[which(rownames(data_Normalized$E) %in% rownames(Vitdvscontrol_sel_signif)),]

xt<-t(data_Normalized$E)
xts<-scale(xt)
k_sel_s<-t(xts)


### load mouse gtf para meter symbols<-
k_temp<-merge(k_sel_s, mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)


rownames(k_temp)<-k_temp$gene_name
ids <- bitr(k_temp$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

k_temp2<-merge(k_temp, ids, by.x="ENSG", by.y="ENSEMBL", all.y=T)


k_temp$ENSG<-NULL
k_temp$gene_name<-NULL

### cargar rutas kegg


#pathways <- fgsea::gmtPathways("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/imagenes/m2.all.v0.2.symbols.gmt")
pathways <- fgsea::gmtPathways("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/imagenes/m5")


#unir todos los genes que aparezcan en las WNT
genes_patways_wnt<-unique(c(pathways$WP_WNT_SIGNALING_PATHWAY_AND_PLURIPOTENCY, pathways$WP_WNT_SIGNALING_IN_KIDNEY_DISEASE, pathways$WP_WNT_SIGNALING_PATHWAY_NETPATH ,pathways$WP_WNT_SIGNALING_PATHWAY, pathways$REACTOME_SIGNALING_BY_WNT, pathways$REACTOME_WNT_LIGAND_BIOGENESIS_AND_TRAFFICKING, pathways$BIOCARTA_MOUSE_WNT_PATHWAY , pathways$BIOCARTA_MOUSE_WNT_LRP6_PATHWAY, pathways$SANSOM_WNT_PATHWAY_REQUIRE_MYC))

genes_patways_hedgehog<-unique(c(pathways$YAUCH_HEDGEHOG_SIGNALING_PARACRINE_UP, pathways$YAUCH_HEDGEHOG_SIGNALING_PARACRINE_DN, pathways$REACTOME_HEDGEHOG_LIGAND_BIOGENESIS ,pathways$REACTOME_SIGNALING_BY_HEDGEHOG, pathways$REACTOME_HEDGEHOG_OFF_STATE, pathways$REACTOME_HEDGEHOG_ON_STATE, pathways$WP_HEDGEHOG_SIGNALING_PATHWAY ))

genes_GOBP_vitD<-unique(c(pathways$GOBP_REGULATION_OF_VITAMIN_METABOLIC_PROCESS, pathways$GOBP_RESPONSE_TO_VITAMIN_D, pathways$GOBP_CELLULAR_RESPONSE_TO_VITAMIN_D, pathways$GOBP_REGULATION_OF_VITAMIN_D_RECEPTOR_SIGNALING_PATHWAY, pathways$GOBP_VITAMIN_D_RECEPTOR_SIGNALING_PATHWAY, pathways$GOBP_VITAMIN_D_METABOLIC_PROCESS, pathways$GOBP_NEGATIVE_REGULATION_OF_VITAMIN_D_BIOSYNTHETIC_PROCESS, pathways$GOBP_VITAMIN_D_METABOLIC_PROCESS, pathways$GOBP_VITAMIN_D_BIOSYNTHETIC_PROCESS, pathways$GOBP_FAT_SOLUBLE_VITAMIN_BIOSYNTHETIC_PROCESS, pathways$GOBP_VITAMIN_TRANSPORT, pathways$GOBP_VITAMIN_D3_METABOLIC_PROCESS, pathways$GOBP_REGULATION_OF_CALCIDIOL_1_MONOOXYGENASE_ACTIVITY, pathways$GOBP_CALCITONIN_FAMILY_RECEPTOR_SIGNALING_PATHWAY))


# solo hay 19 genes que tenemos en el experimento expresados
a<-which(genes_patways_wnt %in% k_temp2$gene_name)
b<-which(genes_patways_hedgehog %in% k_temp2$gene_name)
c<-which(genes_GOBP_vitD %in% k_temp2$gene_name)



for_heatmap_wnt<-genes_patways_wnt[a]
for_heatmap_hedgehog<-genes_patways_hedgehog[b]
for_heatmap_vitD<-genes_GOBP_vitD[c]


k_temp_sel<-k_temp[which(rownames(k_temp) %in% for_heatmap_wnt),]
k_temp_sel<-k_temp[which(rownames(k_temp) %in% for_heatmap_hedgehog),]
k_temp_sel<-k_temp[which(rownames(k_temp) %in% for_heatmap_vitD),]


k_sel_s<-k_temp_sel


col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("VitD", "Control")

colores<-diverge_hsv(2) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 2)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_ABELLA_ALL_EXPRESSED.pdf")

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
    main              = "Zscore of all genes expressed \n VitD vs Control  ",
    cluster_cols = TRUE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


pdf("Heatmap_Genes_ABELLA_VITD_WNT_SIGNALING.pdf")

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Zscored expression of genes related with WNT signaling pathway  \n VitD vs Control  ",
    cluster_cols = TRUE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()

pdf("Heatmap_Genes_ABELLA_VITD_HEDGEHOG_SIGNALING.pdf")

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Zscored expression of genes related with Hedgehog signaling pathway  \n VitD vs Control  ",
    cluster_cols = TRUE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()

pdf("Heatmap_Genes_ABELLA_VITD_VITAMIN_D_RELATED.pdf")

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Zscored expression of genes related with VITD metabolism or transport  \n VitD vs Control  ",
    cluster_cols = TRUE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


library(AnnotationDbi)
library(clusterProfiler)
library(org.Mm.eg.db)
k_sel_s<-as.data.frame(k_sel_s)
#k_sel_s$ENSG<-paste(lapply(strsplit(paste(rownames(k_sel_s)),"\\."),"[",1))


library("org.Mm.eg.db")
library(rtracklayer)


mouse_gtf<-import("/home/jgonzalezgom/references/STAR_mouse_GHRC39/gencode.vM27.annotation.gtf")

mouse_gtf<-as.data.frame(mouse_gtf) 
mouse_gtf_sel<-mouse_gtf[,c("gene_id","gene_name")]
mouse_gtf_sel<-unique(mouse_gtf_sel)

mouse_gtf_selection<-mouse_gtf_sel[which(mouse_gtf_sel$gene_id %in% rownames(k_sel_s)),]

mouse_gtf_selection<-unique(mouse_gtf_selection)

k_sel_s$ENSG<-rownames(k_sel_s)
k_sel_s_annot<-merge(k_sel_s, mouse_gtf_selection, by.x="ENSG", by.y="gene_id", all.x=T)

library(xlsx)
write.xlsx(k_sel_s_annot, file="VitDvsControl_RNAseq_ABELLA_VitD_selección_significativos.xlsx")


VitDvsControl$ENSG<-rownames(VitDvsControl)
VitDvsControl_annot<-merge(VitDvsControl, mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)



### Volcanoplot

library(ggrepel)


volcanoPlots_pepe<-function(dataset,pvalth, FCth,ntop ,title){

    dataset<-dataset[,1:8]
    dataset<-unique(dataset)

    dataset <- dataset %>% 
    mutate(
        Expression = case_when(logFC >= FCth & P.Value < pvalth ~ "Up-regulated",
                            logFC <= -FCth & P.Value < pvalth ~ "Down-regulated",
                            TRUE ~ "Unchanged")
        )



    y_int<-(-1)*log10(pvalth)

    p2 <- ggplot(dataset, aes(logFC, -log(P.Value,10))) +
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
        arrange(P.Value, desc(abs(logFC))) %>% 
        head(top),
    dataset %>% 
        filter(Expression == 'Down-regulated') %>% 
        arrange(P.Value, desc(abs(logFC))) %>% 
        head(top)
    )

    p3<-p2+ geom_label_repel(data = top_genes,
                    mapping = aes(logFC, -log(P.Value,10), label = gene_name),
                    size = 2, max.overlaps = Inf)
    
    return(p3)
}

VitDvsControl_plot<-volcanoPlots_pepe(VitDvsControl_annot,0.05,1,20, "VitD vs Control")

pdf("VOLCANO_PLOTS_VitDvsControl_pval0.05_RNAseq_ABELLA_DIC21_2.pdf")
#plot(PBSvsNaive_plot)
#plot(CombovsPBS_plot)
plot(VitDvsControl_plot)

dev.off()


VitDvsControl_annot_order<-VitDvsControl_annot[order(VitDvsControl_annot$t, decreasing=T),]
VitDvsControl_annot_order$ENSG<-paste(lapply(strsplit(paste(VitDvsControl_annot_order$ENSG),"\\."),"[",1))

ids <- bitr(VitDvsControl_annot_order$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

VitDvsControl_annot_order_s<-merge(VitDvsControl_annot_order, ids, by.x="ENSG", by.y="ENSEMBL", all.y=T)
# ordenar por el estadistico 

VitDvsControl_annot_order_s_0.05<-VitDvsControl_annot_order_s[which(VitDvsControl_annot_order_s$P.Value < 0.05),]

VitDvsControl_annot_order_s_signif<-VitDvsControl_annot_order_s_0.05[which(VitDvsControl_annot_order_s_0.05$logFC < -1 | VitDvsControl_annot_order_s_0.05$logFC > 1),]


enriched_VitDvsControl_selUP<-na.omit(unique(VitDvsControl_annot_order_s_signif[which( VitDvsControl_annot_order_s_signif$logFC >1), "ENTREZID"]))

enriched_GOBP_VitDvsControl <- enrichGO(enriched_VitDvsControl_selUP, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_VitDvsControl_s <- simplify(enriched_GOBP_VitDvsControl, cutoff=0.5, by="p.adjust", select_fun=min)   # UP

enriched_GOBP_VitDvsControl_df<-as.data.frame(enriched_GOBP_VitDvsControl)
enriched_GOBP_VitDvsControl_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_VitDvsControl$p.adjust)))


enriched_VitDvsControl_sel_down<-na.omit(unique(VitDvsControl_annot_order_s_signif[which(VitDvsControl_annot_order_s_signif$logFC < -1), "ENTREZID"]))

enriched_GOBP_VitDvsControl_down <- enrichGO(enriched_VitDvsControl_sel_down, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_VitDvsControl_down_s <- simplify(enriched_GOBP_VitDvsControl_down, cutoff=0.5, by="p.adjust", select_fun=min)   # DOWN

a<-mutate(enriched_GOBP_VitDvsControl_s, logp =(-1)*log(p.adjust, base=10)) %>% 
    dotplot(x="logp", showCategory=30) + ggtitle("GOBP enrichment UPGENES with pval < 0.05")
b<-mutate(enriched_GOBP_VitDvsControl_down_s, logp =(-1)*log(p.adjust, base=10)) %>% 
    dotplot(x="logp", showCategory=30) + ggtitle("GOBP enrichment DOWNGENES with pval < 0.05")

pdf("GOBP_DOTPLOT_ordered_by_logPval_simplified_pathways.pdf",10,10)
plot(a)
plot(b)
dev.off()

pdf("GOBP_DOTPLOT_ordered_by_GeneRatio_simplified_pathways.pdf",10,10)
dotplot(enriched_GOBP_VitDvsControl_s, showCategory=30) + ggtitle("GOBP enrichment UPGENES with pval < 0.05")
dotplot(enriched_GOBP_VitDvsControl_down_s, showCategory=30) + ggtitle("GOBP enrichment DOWNGENES with pval < 0.05")

dev.off()



enriched_GOBP_VitDvsControl_down_df<-as.data.frame(enriched_GOBP_VitDvsControl_down)
enriched_GOBP_VitDvsControl_down_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_VitDvsControl_down_df$p.adjust)))


sort1_up <- enriched_GOBP_VitDvsControl_df[order(enriched_GOBP_VitDvsControl_df[, "logp"], decreasing = TRUE), ]
sort1_up_sel<-sort1_up[1:50,]

sort1_down <- enriched_GOBP_VitDvsControl_down_df[order(enriched_GOBP_VitDvsControl_down_df[, "logp"], decreasing = TRUE), ]
sort1_down_sel<-enriched_GOBP_VitDvsControl_down_df[1:50,]

sort1_down_sel$logp<-sort1_down_sel$logp *(-1)


fgseaRes_sel =rbind(sort1_up_sel,sort1_down_sel)


fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$logp > 0, "Up-regulated", "Down-regulated")
filtRes = fgseaRes_sel

upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
downcols_rev<-rev(downcols)
colos = c(upcols, downcols_rev)
names(colos) = 1:length(colos)
filtRes$Index = as.factor(1:nrow(filtRes))

title_plot<-"ENRICHMENT GOBP VitDvsControl DEgenes"

g_fgseaRes = ggplot(filtRes, aes(reorder(Description, logp), logp)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="logP",
        title=title_plot) + 
    theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 20)) 


pdf("enrichment_GOBP_VITDvsControl_2.pdf",40,40)
plot(g_fgseaRes)
dev.off()


### GOMF

enriched_VitDvsControl_selUP<-na.omit(unique(VitDvsControl_annot_order_s_signif[which( VitDvsControl_annot_order_s_signif$logFC >1), "ENTREZID"]))

enriched_GOMF_VitDvsControl <- enrichGO(enriched_VitDvsControl_selUP, 'org.Mm.eg.db', ont="MF", pvalueCutoff=0.01)
enriched_GOMF_VitDvsControl_s <- simplify(enriched_GOMF_VitDvsControl, cutoff=0.7, by="p.adjust", select_fun=min)   # UP

enriched_GOMF_VitDvsControl_df<-as.data.frame(enriched_GOMF_VitDvsControl)
enriched_GOMF_VitDvsControl_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOMF_VitDvsControl$p.adjust)))



enriched_GOMF_VitDvsControl_down <- enrichGO(enriched_VitDvsControl_sel_down, 'org.Mm.eg.db', ont="MF", pvalueCutoff=0.01)
enriched_GOMF_VitDvsControl_down_s <- simplify(enriched_GOMF_VitDvsControl_down, cutoff=0.7, by="p.adjust", select_fun=min)   # DOWN

enriched_GOMF_VitDvsControl_down_df<-as.data.frame(enriched_GOMF_VitDvsControl_down)
enriched_GOMF_VitDvsControl_down_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOMF_VitDvsControl_down_df$p.adjust)))


sort1_up <- enriched_GOMF_VitDvsControl_df[order(enriched_GOMF_VitDvsControl_df[, "logp"], decreasing = TRUE), ]
sort1_up_sel<-sort1_up[1:10,]

sort1_down <- enriched_GOMF_VitDvsControl_down_df[order(enriched_GOMF_VitDvsControl_down_df[, "logp"], decreasing = TRUE), ]
sort1_down_sel<-enriched_GOMF_VitDvsControl_down_df[1:10,]

sort1_down_sel$logp<-sort1_down_sel$logp *(-1)


fgseaRes_sel =rbind(sort1_up_sel,sort1_down_sel)


fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$logp > 0, "Up-regulated", "Down-regulated")
filtRes = fgseaRes_sel

upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
downcols_rev<-rev(downcols)
colos = c(upcols, downcols_rev)
names(colos) = 1:length(colos)
filtRes$Index = as.factor(1:nrow(filtRes))

title_plot<-"ENRICHMENT GOMF VitDvsControl DEgenes"

g_fgseaRes = ggplot(filtRes, aes(reorder(Description, logp), logp)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="logP",
        title=title_plot) + 
    theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 20)) 


pdf("enrichment_GOMF_VITDvsControl.pdf",40,40)
plot(g_fgseaRes)
dev.off()


## kegg overrepresentation analysis

### GSEA 
VitDvsControl_annot_order_s

ranks<-setNames(VitDvsControl_annot_order_s$t, VitDvsControl_annot_order_s$ENTREZID)

enriched_VitDvsControl_selUP
enriched_VitDvsControl_sel_down

##KEGG pathway overrepresentation analisis (OVA)

kkup <- enrichKEGG(gene         = enriched_VitDvsControl_selUP,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

kkdown <- enrichKEGG(gene         = enriched_VitDvsControl_sel_down,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)

pdf("KEGG_enrichment_pathways_only_genes_DEs_with_pval_0.05.pdf",10,10)
dotplot(kkup, showCategory=30) + ggtitle("KEGG pathways enrichment UPGENES with pval < 0.05")
dotplot(kkdown, showCategory=30) + ggtitle("KEGG pathways enrichment DOWNGENES with pval < 0.05")
dev.off()


      


### Kegg pathway gene set enrichment analysis


VitDvsControl_annot_entrez_ord<-VitDvsControl_annot_order_s[order(VitDvsControl_annot_order_s$t, decreasing=T),]


output<-data.frame()

for(symbol in unique(VitDvsControl_annot_order_s_s$ENTREZID)){
    cat(symbol,"\n")
    selection<-VitDvsControl_annot_order_s_s[which(VitDvsControl_annot_order_s_s$ENTREZID==symbol),]
    selection<-selection[1,]
    output<-rbind(output,selection)
}

# preparedata_Kegg<-function(dataframe){
#     listed_rank<-list()
#     for(i in 1:nrow(dataframe)){
#         symbol<-dataframe[,"ENTREZID"][i]
#         #symbol<-toupper(symbol)
#         t_val<-dataframe[,"t"][i]
#         listed_rank[[symbol]]<-t_val
#     }
#     listed_rank2<-unlist(listed_rank)
#     return(listed_rank2)
# }

# VitDvsControl_prepared_prepared<-preparedata_Kegg(VitDvsControl_annot_entrez_ord)


library(fgsea)

#pathways <- gmtPathways("/home/jgonzalezgom/05_ABELLA/RNAseq_DIC21/imagenes/c2.all.v7.4.entrez.gmt")

ranks <- setNames(output$t, output$ENTREZID)

### GSEA CON OTRO PAQUETE LO MISMO QUE CDITRANI


##### GSEA

####
organism = "org.Mm.eg.db"

gse <- gseGO(geneList=ranks, 
             ont ="BP", 
             keyType = "ENTREZID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             nPermSimple = 10000,
             pAdjustMethod = "fdr")

a<-plotGOgraph(gse, firstSigNodes = 20,
sigForAll = FALSE,
useFullNames = TRUE)

pdf("plotGO_GSEA.pdf")
plot(a$complete.dag)
dev.off()

require(DOSE)

pdf("ENRICHGOBP_DOTPLOT_GSEA.pdf")
dotplot(gse, showCategory=20, split=".sign", font.size=6, label_format=20) + facet_grid(.~.sign)
dev.off()

max.overlaps = getOption("ggrepel.max.overlaps", default = 70)

pdf("ENRICHGOBP_NETPLOT.pdf")
cnetplot(gse, categorySize="pvalue", foldChange=ranks, showCategory = 3)
dev.off()

#GSEA top 10
pdf("GSEA_ENRICHMENTS.pdf")
gseaplot(gse, by = "preranked", title = gse$Description[1], geneSetID = 1)
gseaplot(gse, by = "preranked", title = gse$Description[2], geneSetID = 2)
gseaplot(gse, by = "preranked", title = gse$Description[3], geneSetID = 3)
gseaplot(gse, by = "preranked", title = gse$Description[4], geneSetID = 4)
gseaplot(gse, by = "preranked", title = gse$Description[5], geneSetID = 5)
gseaplot(gse, by = "preranked", title = gse$Description[6], geneSetID = 6)
gseaplot(gse, by = "preranked", title = gse$Description[7], geneSetID = 7)
gseaplot(gse, by = "preranked", title = gse$Description[8], geneSetID = 8)
gseaplot(gse, by = "preranked", title = gse$Description[9], geneSetID = 9)
gseaplot(gse, by = "preranked", title = gse$Description[10], geneSetID = 10)


dev.off()

#### KEGG gene set enrichment



kegg_organism = "mmu"

kk2 <- gseKEGG(geneList     = ranks,
               organism     = kegg_organism,
               minGSSize    = 30,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               eps=0,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

pdf("dotplot_kegg_pathways_nuevo.pdf")
dotplot(kk2, showCategory = 20, title = "Enriched KEGG Pathways" , font.size=6, label_format=50, split=".sign") + facet_grid(.~.sign)
dev.off()


kk2df<-as.data.frame(kk2)
save(kk2df, file="KEGG_ANALYSIS.Rdata")


genes<-kk2df[which(kk2df$Description=="Wnt signaling pathway"),"core_enrichment"]

genes_splitted_entrez<-unlist(strsplit(genes, split = "/"))

genes_entrez_WNT<-data.frame("ENTREZ"=genes_splitted_entrez)

ids<-bitr(genes_splitted_entrez, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Mm.eg.db")
ids2<-bitr(genes_splitted_entrez, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")

ids3<-bitr(ids2$ENSEMBL, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")

#ids <- bitr(k_temp$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

save(ids3, file="Genes_WNT_signaling_pathway.Rdata")

#load(dataNormalized)

rownames(data_Normalized)<-paste(lapply(strsplit(rownames(data_Normalized),"\\."),"[",1))


#WNT heatmap del pathway enriquecido

selection<-data_Normalized[which(rownames(data_Normalized) %in% ids3$ENSEMBL),]
selection<-as.data.frame(selection)
selection$ENSG<-rownames(selection)
selection_s<-merge(selection, ids3, by.x="ENSG", by.y="ENSEMBL")
rownames(selection_s)<-selection_s$SYMBOL
selection_s$ENSG<-NULL
selection_s$SYMBOL<-NULL

#zscoring

xt<-t(selection_s)
xts<-scale(xt)
k_sel_s<-t(xts)

col_groups<-snames
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("VitD", "Control")

colores<-diverge_hsv(4) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 4)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_pathway_WNT.pdf")

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = FALSE,
    fontsize          = 6,
    main              = "Zscore of WNT signaling pathway \n VitD vs Control  ",
    cluster_cols = TRUE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
)
dev.off()


pathways<-gmtPathways("/home/jgonzalezgom/05_ABELLA/m2.all.v0.2.entrez.gmt")

fgseaResVitDvsC <- fgsea(pathways = pathways, stats    = ranks ,minSize=10, maxSize=500, eps=0)

dotplot(fgseaResVitDvsC)

plot_GSEA_PEPE<-function(fgseaRes, gene_set_prepared, pathways_selected, sizemin, sizemax, n_plot, title_plot, pval){

## Filter FGSEA by using gage results. Must be significant and in same direction to keep 

    gaRes = gage::gage(gene_set_prepared, gsets=pathways_selected, same.dir=TRUE, set.size =c(sizemin,sizemax))

    ups = as.data.frame(gaRes$greater) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")

    downs = as.data.frame(gaRes$less) %>% 
    tibble::rownames_to_column("Pathway") %>% 
    dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
    dplyr::select("Pathway")
 
    #keepups= topPathways_up 
    #keepdowns =  topPathways_down
    keepups = fgseaRes[fgseaRes$NES > 0 & !is.na(match(fgseaRes$pathway, ups$Pathway)), ]
    keepdowns = fgseaRes[fgseaRes$NES < 0 & !is.na(match(fgseaRes$pathway, downs$Pathway)), ]

    fgRes = fgseaRes[ !is.na(match(fgseaRes$pathway, 
           c( keepups$pathway, keepdowns$pathway))), ] %>% 
        arrange(desc(NES))
    
    fgRes$pathway = stringr::str_replace(fgRes$pathway, "WP_" , "")
    fgRes$pathway = stringr::str_replace_all(fgRes$pathway, "_" , " ")

    fgRes$Enrichment<-ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
    fgRes<-unique(as.data.frame(fgRes))
    
    filtRes = rbind(head(fgRes, n = n_plot),
                  tail(fgRes, n = n_plot ))
    filtRes = unique(filtRes)

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
      theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 5)) 

    return(g_fgseaRes)
}

VitDvsC_plot<-plot_GSEA_PEPE(fgseaResVitDvsC, ranks, pathways, 10, 500, 30, "GSEA M2 (cp) VitD vs Control", 0.05 )


pdf("GSEA_VitDvsC_plot_m2_pathways.pdf",10,10)
plot(VitDvsC_plot)
dev.off()


# GSEA with dotplot fgseaResVitDvsC

fgseaResVitDvsC$type = "upregulated"
fgseaResVitDvsC$type[fgseaResVitDvsC$NES < 0] = "downregulated"
fgseaResVitDvsC$GeneRatio<-as.numeric(fgseaResVitDvsC$size)/length(ranks)

fgseaResVitDvsC_UP<-fgseaResVitDvsC[which(fgseaResVitDvsC$type=="upregulated"),]
fgseaResVitDvsC_down<-fgseaResVitDvsC[which(fgseaResVitDvsC$type=="downregulated"),]

fgseaResVitDvsC_UP_ord<-fgseaResVitDvsC_UP[order(fgseaResVitDvsC_UP$padj),]
fgseaResVitDvsC_down_ord<-fgseaResVitDvsC_down[order(fgseaResVitDvsC_down$padj),]

top_N<30

fgseaResVitDvsC_UP_ord_s<-fgseaResVitDvsC_UP_ord[1:top_N,]
fgseaResVitDvsC_down_ord_s<-fgseaResVitDvsC_down_ord[1:top_N,]

plotter<-rbind(fgseaResVitDvsC_UP_ord_s,fgseaResVitDvsC_down_ord_s )

p <- ggplot(plotter, aes(x = GeneRatio, y = reorder(pathway, GeneRatio))) + 
               geom_point(aes(size = GeneRatio, color = padj)) +
               theme_bw(base_size = 5) +
        scale_colour_gradient(limits=c(0, 0.10), low="red") +
        ylab(NULL) +
        ggtitle("GSEA (M2) pathway enrichment") + facet_grid(.~type)
pdf("GSEA_M2_GeneRatio.pdf")
print(p)
dev.off()

       