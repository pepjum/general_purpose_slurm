library(dplyr)
library(edgeR)
library("genefilter")
library(RColorBrewer)
library(viridis)

source("/home/jgonzalezgom/01_NMELERO/funcionesVikv2.R")


files<-list.files("/home/jgonzalezgom/06_CDITRANI/counts/", pattern=".txt$")
files<-paste0("/home/jgonzalezgom/06_CDITRANI/counts/", files)


samples<-c()
data_df<-data.frame()
for(i in 1:length(files)){
    cat(paste(files[i]),"\n")
    fileloaded<-read.table(files[i], header=T)
    sample_tmp<-basename(files[i])
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
system("mkdir -p /home/jgonzalezgom/06_CDITRANI/dataframes/")
save(data_df, file="/home/jgonzalezgom/06_CDITRANI/dataframes/RNASEQ_CDITRANI_unordered_raw.Rdata")


samples_PBS<-which(grepl("PBS", samples))
samples_IL12IP<-which(grepl("IP-", samples))
samples_IL12IV<-which(grepl("IV-", samples))


data_df_ordered<-data_df[,c(samples_PBS,samples_IL12IP, samples_IL12IV)]
save(data_df_ordered, file="/home/jgonzalezgom/06_CDITRANI/dataframes/RNASEQ_CDITRANI_ordered_raw.Rdata")

#not in
`%!in%` <- Negate(`%in%`)

not_sel<-c("IP-1_S5","IV-3_S11","PBS-4_S4")

data_df_selected<-data_df[,colnames(data_df) %!in% not_sel ]

samples_PBS<-which(grepl("PBS", colnames(data_df_selected)))
samples_IL12IP<-which(grepl("IP-", colnames(data_df_selected)))
samples_IL12IV<-which(grepl("IV-", colnames(data_df_selected)))


data_df_ordered<-data_df_selected[,c(samples_PBS,samples_IL12IP, samples_IL12IV)]

save(data_df_ordered, file="/home/jgonzalezgom/06_CDITRANI/dataframes/RNASEQ_CDITRANI_ordered_selected_samples_raw.Rdata")


selection_PBS<-colnames(data_df_ordered)[1:3]
selection_IL12IP<-colnames(data_df_ordered)[4:6]
selection_IL12IV<-colnames(data_df_ordered)[7:9]

#data_df_ordered<-data_df_selected

log_intensity_threshold<-5

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.PBS<-genefilter(data_df_ordered[,which(colnames(data_df_ordered) %in% selection_PBS)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.IL12IV<-genefilter(data_df_ordered[,which(colnames(data_df_ordered) %in% selection_IL12IV)], ffun1)

numsamples<-2
f1 <- kOverA(numsamples, log_intensity_threshold)
ffun1 <- filterfun(f1)

whichFilter.IL12IP<-genefilter(data_df_ordered[,which(colnames(data_df_ordered) %in% selection_IL12IP)], ffun1)



whichFilter <- (whichFilter.PBS | whichFilter.IL12IP | whichFilter.IL12IV  )   # no sale

selection_data_filter<-data_df_ordered[whichFilter,]


group=c(rep("PBS", 3), rep("IL12IP", 3), rep("IL12IV", 3))

d0<-DGEList(selection_data_filter, group=group)


### filter low expression genes
 
#voom transformation and calculation variance weights

snames<-group

design<-model.matrix(~0+snames)

system("mkdir -p /home/jgonzalezgom/06_CDITRANI/imagenes/")

setwd("/home/jgonzalezgom/06_CDITRANI/imagenes/")
#rownames(design) <- colnames(d)
#colnames(design) <- levels(snames)

pdf("voom.pdf")
data_Normalized <- voom(d0, design, plot = F)
dev.off()

k<-data_Normalized$E   #20498 9

#save(k, file="/media/inmuno/04_CDITRANI/data_Normalized_CDITRANI_NOV21_selected_samples.Rdata")
save(k, file="/home/jgonzalezgom/06_CDITRANI/dataframes/data_Normalized_CDITRANI_DIC21_selected_samples_quitadas.Rdata")

xlim_r = c(min(na.omit(log2(d0$counts+1))),max(na.omit(log2(d0$counts+1)))) 
xlim = c(min(na.omit(data_Normalized$E)),max(na.omit(data_Normalized$E)))
clust.euclid.average <- hclust(dist(t(log2(d0$counts+1))),method="average")
clust.euclid.average2 <- hclust(dist(t(data_Normalized$E)),method="average")


setwd("/home/jgonzalezgom/06_CDITRANI/imagenes/")

pdf(file="RNASEQ_CDITRANI_IL12IP_QC.pdf", colormode="rgb", width=20, height=20)

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
pdf(file = "PCA_RNASeq_CDITRANI_IL12IP_DIC21.pdf", width = 6, height = 6, colormodel = "rgb")

ggp
dev.off()




fit <- lmFit(data_Normalized, design)

#make contrasts
contr <- makeContrasts(                       
                       IL12IPvsPBS= snamesIL12IP - snamesPBS,
                       IL12IVvsPBS= snamesIL12IV - snamesPBS,
                       IL12IPvsIL12IV= snamesIL12IP - snamesIL12IV,
                       
                       levels = colnames(coef(fit)))


#estimate contrast for each gene

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

system("mkdir -p /home/jgonzalezgom/06_CDITRANI/DEgenes/")
setwd("/home/jgonzalezgom/06_CDITRANI/DEgenes/")


IL12IPvsPBS <- topTable(tmp, coef = "IL12IPvsPBS", n=nrow(tmp), adjust="fdr")
write.table(IL12IPvsPBS, file = "IL12IPvsPBS.txt", quote = FALSE, row.names = TRUE, sep="\t")

IL12IVvsPBS <- topTable(tmp, coef = "IL12IVvsPBS", n=nrow(tmp), adjust="fdr")
write.table(IL12IVvsPBS, file = "IL12IVvsPBS.txt", quote = FALSE, row.names = TRUE, sep="\t")

IL12IPvsIL12IV <- topTable(tmp, coef = "IL12IPvsIL12IV", n=nrow(tmp), adjust="fdr")
write.table(IL12IPvsIL12IV, file = "IL12IPvsIL12IV.txt", quote = FALSE, row.names = TRUE, sep="\t")


IL12IPvsPBS_pval<-IL12IPvsPBS[which(IL12IPvsPBS$adj.P.Val < 0.05 & (IL12IPvsPBS$logFC < (-1) | IL12IPvsPBS$logFC > 1 )), ]
IL12IVvsPBS_pval<-IL12IVvsPBS[which(IL12IVvsPBS$adj.P.Val < 0.05 & (IL12IVvsPBS$logFC < (-1) | IL12IVvsPBS$logFC > 1 )), ]


library("Vennerable")

myLists<-list("IL12-IP vs PBS"=paste(unique(rownames(IL12IPvsPBS_pval))),"IL12-IV vs PBS"=paste(unique(rownames(IL12IVvsPBS_pval))))

Vstem <- Venn(myLists)

pdf("Venn_weighted_CDITRANI.pdf")
plot(Vstem, doWeights = TRUE)
dev.off()



library(pheatmap)
library(colorspace)

setwd("/home/jgonzalezgom/06_CDITRANI/imagenes/")


tmp <- detOutliers(data_Normalized$E, type, "DetOutliers_RNA_SEQ_CDITRANI_DIC21.pdf", 14, 14)


IL12IPvsIL12IV_sel<-IL12IPvsIL12IV[which(IL12IPvsIL12IV$adj.P.Val < 0.05),]
IL12IPvsIL12IV_sel_signif<-IL12IPvsIL12IV_sel[which(IL12IPvsIL12IV_sel$logFC < -1 | IL12IPvsIL12IV_sel$logFC > 1 ),]

selection_genes<-k[which(rownames(k) %in% rownames(IL12IPvsIL12IV_sel_signif)),]

xt<-t(selection_genes)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-as.data.frame(k_sel_s)
#k_sel_s$ENSG<-rownames(k_sel_s)

k_sel_s2<-merge(k_sel_s,mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)
rownames(k_sel_s2)<-k_sel_s2$gene_name

### load mouse gtf para meter symbols<-
k_sel_s<-k_sel_s2[,c(2:4,8:10,5:7)]

#col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12-IV",3),rep("IL12-IP",3))
mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("black","blue","orange"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12-IV","IL12-IP")

#colores<-diverge_hcl(6) #blue, white,red
#mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_pval0.05_logFC1_zscored_nuevo_logFC_1_revisado.pdf",4,52)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    show_colnames  =FALSE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 5,
    main              = "DEg GENES IL12IP vs IL12IV \n pval < 0.05 & logFC <-1 | > 1  ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=4, 
    cellwidth = 10,

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
write.xlsx(k_sel_s_annot, file="RNAseq_CDITRANI_IL12IP_selección_significativos_matrix.xlsx")
write.table(k_sel_s_annot, file="RNAseq_CDITRANI_IL12IP_selección_significativos_matrix.txt", col.names=T, row.names=F, quote=F, sep="\t")



IL12IPvsIL12IV$ENSG<-rownames(IL12IPvsIL12IV)

IL12IPvsIL12IV_annot<-merge(IL12IPvsIL12IV, mouse_gtf_sel, by.x="ENSG", by.y="gene_id", all.x=T)

IL12IPvsIL12IV_annot$ENSG<-paste(lapply(strsplit(paste(IL12IPvsIL12IV_annot$ENSG),"\\."),"[",1))


library("org.Mm.eg.db")
library(clusterProfiler)

ids <- bitr(IL12IPvsIL12IV_annot$ENSG, fromType="ENSEMBL", toType="GENENAME", OrgDb="org.Mm.eg.db")   # universo

IL12IPvsIL12IV_annot_desc<-merge(IL12IPvsIL12IV_annot, ids, by.x="ENSG", by.y="ENSEMBL", all.x=T)


library("org.Mm.eg.db")
library(clusterProfiler)

ids <- bitr(IL12IPvsIL12IV$ENSG, fromType="ENSEMBL", toType="GENENAME", OrgDb="org.Mm.eg.db")   # universo
IL12IPvsIL12IV_annot<-merge(IL12IPvsIL12IV, ids, by.x="ENSG", by.y="ENSEMBL", all.x=T)



write.xlsx(IL12IPvsIL12IV_annot_desc, file="RNAseq_CDITRANI_IL12IPvsIL12IV_DEGs.xlsx")


### Volcanoplot

library(ggrepel)


volcanoPlots_pepe<-function(dataset,pvalth, FCth,ntop ,title, ymin, ymax,xmin,xmax){

    if(ntop==0){
     dataset<-dataset   
    }else{
        dataset<-dataset[,1:8]
    }
    
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
    scale_y_continuous(limits=c(ymin, ymax))+
    scale_x_continuous(limits=c(xmin, xmax))+
    
    scale_color_manual(values = c("blue", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
    geom_vline(xintercept=FCth) + geom_vline(xintercept=-(FCth))+
    geom_hline(yintercept=y_int)+
    ggtitle(title)+ theme(plot.title=element_text( size=8, face='bold'))


    if(ntop >0 ){
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
    }
    if(ntop >0 ){
    p3<-p2+ geom_label_repel(data = top_genes,
                    mapping = aes(logFC, -log(P.Value,10), label = gene_name),
                    size = 2, max.overlaps = Inf)
    
    return(p3)
    }else{
        return(p2)
    }
}


volcanoPlots_pepe2<-function(dataset,pvalth, FCth,title, ymin, ymax,xmin,xmax){

    
    dataset<-dataset   
    
    
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
    ylab(expression("-log"[10]*"(adj.P.Value)")) +
    scale_y_continuous(limits=c(ymin, ymax))+
    scale_x_continuous(limits=c(xmin, xmax))+
    
    scale_color_manual(values = c("blue", "gray50", "firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5))) + theme_classic() +
    geom_vline(xintercept=FCth) + geom_vline(xintercept=-(FCth))+
    geom_hline(yintercept=y_int)+
    ggtitle(title)+ theme(plot.title=element_text( size=8, face='bold'))


    

    top_genes<-dataset[which(dataset$SYMBOL %in% c("Ifng","Gzmb","Cd274","Cxcl9","Cxcl10","Ido1","Cxcl11","Cd209f")),]
    

    p3<-p2+ geom_text_repel(data = top_genes, box.padding=1,
                    mapping = aes(logFC, -log(adj.P.Val,10), label = SYMBOL),
                    size = 5, max.overlaps = Inf)
    
    return(p3)
}





IL12IPvsIL12IV_plot<-volcanoPlots_pepe(IL12IPvsIL12IV,0.05,1,0, "IL12IP vs IL12IV", 0,10,-11,11)
IL12IPvsPBS_plot<-volcanoPlots_pepe2(IL12IPvsPBS_annot_s,0.05,1, "IL12IP vs PBS", 0,5,-10.5,10.5)
IL12IVvsPBS_plot<-volcanoPlots_pepe2(IL12IVvsPBS_annot_s,0.05,1, "IL12IV vs PBS", 0,5,-10.5,10.5)


pdf("VOLCANO_PLOTS_NUEVOS_MISMA_ESCALA_pval0.05_RNAseq_CDITRANI_DIC21_nuevos_w_genesIPPBS.pdf",5,5)
plot(IL12IPvsPBS_plot)
dev.off()

pdf("VOLCANO_PLOTS_NUEVOS_MISMA_ESCALA_pval0.05_RNAseq_CDITRANI_DIC21_nuevos_w_genes_IVPBS.pdf",5,5)
plot(IL12IVvsPBS_plot)
dev.off()

##### GSEA

IL12IPvsIL12IV$ENSG<-paste(lapply(strsplit(rownames(IL12IPvsIL12IV),"\\."),"[",1))
IL12IPvsPBS$ENSG<-paste(lapply(strsplit(rownames(IL12IPvsPBS),"\\."),"[",1))
IL12IVvsPBS$ENSG<-paste(lapply(strsplit(rownames(IL12IVvsPBS),"\\."),"[",1))


ids <- bitr(IL12IPvsIL12IV$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo
IL12IPvsIL12IV_annot_s<-merge(IL12IPvsIL12IV, ids, by.x="ENSG", by.y="ENSEMBL")

 ids <- bitr(IL12IPvsPBS$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo
 IL12IPvsPBS_annot_s<-merge(IL12IPvsPBS, ids, by.x="ENSG", by.y="ENSEMBL")

 ids <- bitr(IL12IVvsPBS$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo
 IL12IVvsPBS_annot_s<-merge(IL12IVvsPBS, ids, by.x="ENSG", by.y="ENSEMBL")


IL12IPvsIL12IV_annot_s_order<-IL12IPvsIL12IV_annot_s[order(IL12IPvsIL12IV_annot_s$t, decreasing=T),]

IL12IPvsIL12IV_annot_s_order<-IL12IPvsIL12IV_annot[order(IL12IPvsIL12IV_annot$t, decreasing=T),]


output<-data.frame()
for(symbol in unique(IL12IPvsIL12IV_annot_s_order$SYMBOL)){
    cat(symbol,"\n")
    selection<-IL12IPvsIL12IV_annot_s_order[which(IL12IPvsIL12IV_annot_s_order$SYMBOL==symbol),]
    selection<-selection[1,]
    output<-rbind(output,selection)
}


#### GSEA

library(fgsea)
library(gage)

ranks <- setNames(output$t, output$SYMBOL)

pathways <- gmtPathways("/home/jgonzalezgom/06_CDITRANI/m2.cp.wikipathways.v0.2.symbols.gmt")

fgseaResIL12IPvsIL12IV <- fgsea(pathways = pathways, stats    = ranks ,minSize=10, maxSize=500, eps=0)
topPathwaysUp_IL12IPvsIL12IV <- fgseaResIL12IPvsIL12IV[ES > 0][head(order(padj), n=200), pathway]
topPathwaysDown_IL12IPvsIL12IV <- fgseaResIL12IPvsIL12IV[ES < 0][head(order(padj), n=200), pathway]
#topPathways_CombovsNaive <- c(topPathwaysUp_CombovsNaive, rev(topPathwaysDown_CombovsNaive))

topPathways_IL12IPvsIL12IV <- c(topPathwaysUp_IL12IPvsIL12IV, rev(topPathwaysDown_IL12IPvsIL12IV))

write.xlsx(topPathways_IL12IPvsIL12IV, file="UP_and_Down_Pathways_IL12_DRvsIL12_m2.cp.KEGG_pathways.xlsx")


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
      theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 20)) 

    return(g_fgseaRes)
}

IL12IPvsIL12IV_plot<-plot_GSEA_PEPE(fgseaResIL12IPvsIL12IV, ranks, pathways, 10, 500, 30, "GSEA M2 (cp) IL12IPvsIL12IV", 0.05 )


pdf("GSEA_IL12IPvsIL12IV_m2_wikipathways.pdf",30,30)
plot(IL12IPvsIL12IV_plot)
dev.off()




# ENRICHGO

IL12IPvsIL12IV_annot_order_s_0.05<-output[which(output$adj.P.Val < 0.05),]

IL12IPvsIL12IV_annot_order_s_signif<-IL12IPvsIL12IV_annot_order_s_0.05[which(IL12IPvsIL12IV_annot_order_s_0.05$logFC < -1 | IL12IPvsIL12IV_annot_order_s_0.05$logFC > 1),]


ids <- bitr(IL12IPvsIL12IV_annot_order_s_signif$ENSG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Mm.eg.db")   # universo

IL12IPvsIL12IV_annot_order_s_signif_entrez<-merge(IL12IPvsIL12IV_annot_order_s_signif, ids, by.x="ENSG", by.y="ENSEMBL")
IL12IPvsIL12IV_annot_order_s_signif_entrez<-unique(IL12IPvsIL12IV_annot_order_s_signif_entrez)


enriched_IL12IPvsIL12IV_selUP<-na.omit(unique(IL12IPvsIL12IV_annot_order_s_signif_entrez[which( IL12IPvsIL12IV_annot_order_s_signif_entrez$logFC >1), "ENTREZID"]))

enriched_GOBP_IL12IPvsIL12IV <- enrichGO(enriched_IL12IPvsIL12IV_selUP, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_IL12IPvsIL12IV_s <- simplify(enriched_GOBP_IL12IPvsIL12IV, cutoff=0.7, by="p.adjust", select_fun=min)   # UP

#enriched_GOBP_IL12IPvsIL12IV_df<-as.data.frame(enriched_GOBP_IL12IPvsIL12IV)
#enriched_GOBP_IL12IPvsIL12IV_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_IL12IPvsIL12IV$p.adjust)))

a<-mutate(enriched_GOBP_IL12IPvsIL12IV_s, logp =(-1)*log(p.adjust, base=10)) %>% dotplot(x="logp", showCategory=20) + ggtitle("GOBP enrichment UPGENES with pval < 0.05")

enriched_IL12IPvsIL12IV_sel_down<-na.omit(unique(IL12IPvsIL12IV_annot_order_s_signif_entrez[which(IL12IPvsIL12IV_annot_order_s_signif_entrez$logFC < -1), "ENTREZID"]))

enriched_GOBP_IL12IPvsIL12IV_down <- enrichGO(enriched_IL12IPvsIL12IV_sel_down, 'org.Mm.eg.db', ont="BP", pvalueCutoff=0.01)
enriched_GOBP_IL12IPvsIL12IV_down_s <- simplify(enriched_GOBP_IL12IPvsIL12IV_down, cutoff=0.7, by="p.adjust", select_fun=min)   # DOWN

#enriched_GOBP_IL12IPvsIL12IV_down_df<-as.data.frame(enriched_GOBP_IL12IPvsIL12IV_down)
#enriched_GOBP_IL12IPvsIL12IV_down_df$logp <- (-1)*log10(as.numeric(paste(enriched_GOBP_IL12IPvsIL12IV_down_df$p.adjust)))

b<-mutate(enriched_GOBP_IL12IPvsIL12IV_down_s, logp =(-1)*log(p.adjust, base=10)) %>% 
    dotplot(x="logp", showCategory=20) + ggtitle("GOBP enrichment DOWNGENES with pval < 0.05")


#3 unir dos enrich objects en 1

enriched_GOBP_IL12IPvsIL12IV_down_s_order<-enriched_GOBP_IL12IPvsIL12IV_down_s[order(enriched_GOBP_IL12IPvsIL12IV_down_s$Count, decreasing=T),]
enriched_GOBP_IL12IPvsIL12IV_down_s_order_sel<-enriched_GOBP_IL12IPvsIL12IV_down_s_order[1:5,]
enriched_GOBP_IL12IPvsIL12IV_down_s_order_sel$group<-"Down-regulated"
enriched_GOBP_IL12IPvsIL12IV_down_s_order_sel$geneRatio<-enriched_GOBP_IL12IPvsIL12IV_down_s_order_sel$Count/316

enriched_GOBP_IL12IPvsIL12IV_s_order<-enriched_GOBP_IL12IPvsIL12IV_s[order(enriched_GOBP_IL12IPvsIL12IV_s$Count, decreasing=T),]
enriched_GOBP_IL12IPvsIL12IV_s_order_sel<-enriched_GOBP_IL12IPvsIL12IV_s_order[1:5,]
enriched_GOBP_IL12IPvsIL12IV_s_order_sel$group<-"Up-regulated"
enriched_GOBP_IL12IPvsIL12IV_s_order_sel$geneRatio<-enriched_GOBP_IL12IPvsIL12IV_s_order_sel$Count/467

list_enr<-list()
list_enr[["a"]]<-enriched_GOBP_IL12IPvsIL12IV_s_order_sel
list_enr[["b"]]<-enriched_GOBP_IL12IPvsIL12IV_down_s_order_sel

plotter<-merge_result(list_enr)
plotter<-as.data.frame(plotter)


pdf("GOBP_DOTPLOT_ordered_by_logPval_simplified_pathways.pdf",10,10)
plot(a)
plot(b)
dev.off()

pdf("GOBP_DOTPLOT_ordered_by_GeneRatio_simplified_pathways_new.pdf",10,10)
dotplot(enriched_GOBP_IL12IPvsIL12IV_s, showCategory=5) + ggtitle("GOBP enrichment UPGENES with pval < 0.05")
dotplot(enriched_GOBP_IL12IPvsIL12IV_down_s, showCategory=5) + ggtitle("GOBP enrichment DOWNGENES with pval < 0.05")

dev.off()

pdf("GOBP_DOTPLOT_ordered_by_GeneRatio_tested.pdf",10,10)
dotplot(binded, showCategory=5, split="group")
dev.off()
pdf("test3.pdf")
ggplot(plotter, aes_string(x="geneRatio2", y="Description", size="Count", color="p.adjust")) +
        geom_point() + scale_size_continuous(range = c(1, 10))+
        scale_color_continuous(low="blue", high="red", name = "logFDR",
                guide=guide_colorbar(reverse=F)) + ylab(NULL) +theme_dose(10) +facet_wrap(~group)
dev.off()

library(forcats)
library(ggplot2)
library(DOSE)

pdf("test4.pdf",10,10)
plotter %>% mutate(Description = fct_reorder(Description, geneRatio)) %>%
ggplot(aes_string(x="geneRatio", y="Description", size="Count", color="group")) +
        geom_point() + scale_size_continuous(range = c(1, 10)) + ylab(NULL) +theme_dose(13) + scale_color_manual(values=c("blue", "red")) + ggtitle("       GOBP enrhichment genes \n adj.pval < 0.05 and logFC < -1 | > 1")
dev.off()


sort1_up <- enriched_GOBP_IL12IPvsIL12IV_df[order(enriched_GOBP_IL12IPvsIL12IV_df[, "logp"], decreasing = TRUE), ]
sort1_up_sel<-sort1_up[1:20,]

sort1_down <- enriched_GOBP_IL12IPvsIL12IV_down_df[order(enriched_GOBP_IL12IPvsIL12IV_down_df[, "logp"], decreasing = TRUE), ]
sort1_down_sel<-enriched_GOBP_IL12IPvsIL12IV_down_df[1:20,]

sort1_down_sel$logp<-sort1_down_sel$logp *(-1)


fgseaRes_sel =rbind(sort1_up_sel,sort1_down_sel)


fgseaRes_sel$Enrichment = ifelse(fgseaRes_sel$logp > 0, "Up-regulated", "Down-regulated")
filtRes = fgseaRes_sel
filtRes<-na.omit(filtRes)

upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
downcols_rev<-rev(downcols)
colos = c(upcols, downcols_rev)
names(colos) = 1:length(colos)

filtRes$Index = as.factor(c(1:length(upcols),rev(sum(length(upcols),length(downcols)-1):length(upcols)+1)))
title_plot<-"Biological Process DEg pval<0.05 IL12IP vs IL12IV"

g_fgseaRes = ggplot(filtRes, aes(reorder(Description, logp), logp)) +
    geom_col( aes(fill = Index )) +
    scale_fill_manual(values = colos ) +
    coord_flip() +
    labs(x="Pathway", y="logP",
        title=title_plot) + 
    theme_minimal() + theme(legend.position = "none") +theme(axis.text = element_text(size = 20)) 


pdf("enrichment_GOBP_IL12IP_vs_IL12IV.pdf",30,30)
plot(g_fgseaRes)
dev.off()


#### barplot top5 up y down
pdf("test.7.pdf")
ggplot(plotterup, aes(geneRatio2, reorder(Description, geneRatio2)) ) +  geom_bar(stat="identity",fill="firebrick3", color="black", width=1) +  labs(x="GeneRatio", y=NULL,
        title=NULL) + 
    theme(aspect.ratio=0.5/1) + theme(legend.position = "none") +theme(axis.text = element_text(size = 10)) 
dev.off()

pdf("test.8.pdf")
ggplot(plotterdown, aes(geneRatio2, reorder(Description, geneRatio2)) ) +  geom_bar(stat="identity",fill="blue", color="black", width=1) +  labs(x="GeneRatio", y=NULL,
        title=NULL) + 
    theme(aspect.ratio=0.5/1) + theme(legend.position = "none") +theme(axis.text = element_text(size = 10)) 
dev.off()

##KEGG pathways solo de los genes diferencialmente expresados

#genes UP

kkup <- enrichKEGG(gene         = enriched_IL12IPvsIL12IV_selUP,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)




write.xlsx(kk, file="KEGG_pathways_Overrepresentation_analysis_VitDvsControl.xlsx")                 

kkdown <- enrichKEGG(gene         = enriched_IL12IPvsIL12IV_sel_down,
                 organism     = 'mmu',
                 pvalueCutoff = 0.05)



pdf("KEGG_enrichment_pathways_only_genes_DEs_with_pval_0.05.pdf",10,10)
dotplot(kkup, showCategory=30) + ggtitle("KEGG pathways enrichment UPGENES with pval < 0.05")
dotplot(kkdown, showCategory=30) + ggtitle("KEGG pathways enrichment DOWNGENES with pval < 0.05")
dev.off()


# HEATMAPS que pide CDITRANI

# leucocyte migration
pathways <- gmtPathways("/home/jgonzalezgom/06_CDITRANI/m2.all.v0.2.symbols.gmt")


leucopath<-pathways$GOBP_POSITIVE_REGULATION_OF_CYTOKINE_PRODUCTION

leuco<-enriched_GOBP_IL12IPvsIL12IV_s[which(enriched_GOBP_IL12IPvsIL12IV_s$Description=="leukocyte migration"),]

leuco_genes<-leuco$geneID

leucogenes<-unlist(strsplit(leuco_genes,"\\/"))

ids1 <- bitr(leucogenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12-IV",3), rep("IL12-IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("black","blue","orange"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12-IV","IL12-IP")

#colores<-diverge_hcl(6) #blue, white,red
#mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_GO_leucocyte_migration_nuevo.pdf",6,6)

pheatmap(k_sel_s,
  #  color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Gene Expression \n 'Leukocyte migration' GO enriched term \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    cellheight=8, 
    cellwidth = 12,
    angle_col=c("45"),
)
dev.off()


# leucocyte migration

cyto<-enriched_GOBP_IL12IPvsIL12IV_s[which(enriched_GOBP_IL12IPvsIL12IV_s$Description=="positive regulation of cytokine production"),]

cytokine_genes<-cyto$geneID

cytokinegenes<-unlist(strsplit(cytokine_genes,"\\/"))

ids1 <- bitr(cytokinegenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


#col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12-IV",3), rep("IL12-IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("black","blue","orange"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12-IV","IL12-IP")

#colores<-diverge_hcl(6) #blue, white,red
#mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_GO_cytokine_production_nuevo.pdf",6,6)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Gene Expression \n 'positive regulation of cytokine production' GO enriched term \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    cellheight=8, 
    cellwidth = 12,

    angle_col=c("45"),
)
dev.off()


reg<-enriched_GOBP_IL12IPvsIL12IV_s[which(enriched_GOBP_IL12IPvsIL12IV_s$Description=="regulation of inflammatory response"),]

reg_genes<-reg$geneID

reggenes<-unlist(strsplit(reg_genes,"\\/"))

ids1 <- bitr(reggenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_GO_regulation_inflammatory_response.pdf")

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
    main              = "Zscored Gene Expression \n 'regulation of inflammatory response' enriched term in GO \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    cellheight=3, 
    cellwidth = 3,

    angle_col=c("45"),
)
dev.off()

tcell<-enriched_GOBP_IL12IPvsIL12IV_s[which(enriched_GOBP_IL12IPvsIL12IV_s$Description=="T cell proliferation"),]

tcell_genes<-tcell$geneID

tcellgenes<-unlist(strsplit(tcell_genes,"\\/"))

ids1 <- bitr(tcellgenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_GO_tcell_proliferation.pdf")

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
    main              = "Zscored Gene Expression \n 'T cell proliferation' enriched term in GO \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    cellheight=3, 
    cellwidth = 3,

    angle_col=c("45"),
)
dev.off()



ifgamma<-enriched_GOBP_IL12IPvsIL12IV_s[which(enriched_GOBP_IL12IPvsIL12IV_s$Description=="response to interferon-gamma"),]

ifgamma_genes<-ifgamma$geneID

ifgammagenes<-unlist(strsplit(ifgamma_genes,"\\/"))

ids1 <- bitr(ifgammagenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_GO_ifgamma_response.pdf")

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
    main              = "Zscored Gene Expression \n 'response to interferon-gamma' enriched term in GO \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=3, 
    cellwidth = 3,

)
dev.off()

fatty<-enriched_GOBP_IL12IPvsIL12IV_down_s[which(enriched_GOBP_IL12IPvsIL12IV_down_s$Description=="fatty acid metabolic process"),]

fatty_genes<-fatty$geneID

fattygenes<-unlist(strsplit(fatty_genes,"\\/"))

ids1 <- bitr(fattygenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-c(rep("PBS",3),rep("IL12-IV",3), rep("IL12-IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("black","blue","orange"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12-IV","IL12-IP")

#colores<-diverge_hcl(6) #blue, white,red
#mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
#              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_GO_fatty_acid_nuevo.pdf",5,4)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "Zscored Gene Expression \n 'fatty acid metabolic process' enriched term in GO \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=8, 
    cellwidth = 10,

)
dev.off()

# heatmaps de KEGG

kkup<-as.data.frame(kkup)


cyto<-kkup[which(kkup$Description=="Cytokine-cytokine receptor interaction"),]

cyto_genes<-cyto$geneID

cytogenes<-unlist(strsplit(cyto_genes,"\\/"))

ids1 <- bitr(cytogenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap
data_Normalized<-data_Normalized$E
rownames(data_Normalized)<-paste(lapply(strsplit(paste(rownames(data_Normalized)),"\\."),"[",1))

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

#colores<-diverge_hcl(6) #blue, white,red
#mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
#myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
#              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_KEGG_cytokine-receptor-interaction.pdf",4,8)

pheatmap(k_sel_s,
    #color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    #breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 8,
    main              = "'Cytokine-cytokine receptor interaction' \n enriched KEGG pathway \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=8, 
    cellwidth = 10,

)
dev.off()


NF<-kkup[which(kkup$Description=="NF-kappa B signaling pathway"),]

NF_genes<-NF$geneID

NFgenes<-unlist(strsplit(NF_genes,"\\/"))

ids1 <- bitr(NFgenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_KEGG_NF-KappaB_signaling_pathway.pdf")

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
    main              = "Zscored Gene Expression \n 'NF-Kappa B signaling pathway' enriched KEGG pathway \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=3, 
    cellwidth = 3,

)
dev.off()

JAK<-kkup[which(kkup$Description=="JAK-STAT signaling pathway"),]

JAK_genes<-JAK$geneID

JAKgenes<-unlist(strsplit(JAK_genes,"\\/"))

ids1 <- bitr(JAKgenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_KEGG_JAK-STAT_signaling_pathway.pdf")

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
    main              = "Zscored Gene Expression \n 'JAK-STAT signaling pathway' enriched KEGG pathway \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    angle_col=c("45"),
    cellheight=3, 
    cellwidth = 3,

)
dev.off()



chemo<-kkup[which(kkup$Description=="Chemokine signaling pathway"),]

chemo_genes<-chemo$geneID

chemogenes<-unlist(strsplit(chemo_genes,"\\/"))

ids1 <- bitr(chemogenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_KEGG_chemokine_signaling_pathway.pdf")

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
    main              = "Zscored Gene Expression \n 'Chemokine signaling pathway' enriched KEGG pathway \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    cellheight=3, 
    cellwidth = 3,

    angle_col=c("45"),
)
dev.off()


TNF<-kkup[which(kkup$Description=="TNF signaling pathway"),]

TNF_genes<-TNF$geneID

TNFgenes<-unlist(strsplit(TNF_genes,"\\/"))

ids1 <- bitr(TNFgenes, fromType="ENTREZID", toType="ENSEMBL", OrgDb="org.Mm.eg.db")   # universo


#cargar df de heatmap

selection_genes<-data_Normalized[which(rownames(data_Normalized) %in% ids1$ENSEMBL),]
selection_genes<-as.data.frame(selection_genes)
selection_genes$ENSG<-rownames(selection_genes)

ids2 <- bitr(selection_genes$ENSG, fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")   # universo

selection_genes2<-merge(selection_genes, ids2, by.x="ENSG", by.y="ENSEMBL", all.x=T)
selection_genes2$ENSG<-NULL
rownames(selection_genes2)<-selection_genes2$SYMBOL
selection_genes2$SYMBOL<-NULL


xt<-t(selection_genes2)
xts<-scale(xt)
k_sel_s<-t(xts)
k_sel_s<-k_sel_s[,c(1:3,7:9,4:6)]

### load mouse gtf para meter symbols<-


col_groups<-snames
col_groups<-c(rep("PBS",3),rep("IL12IV",3), rep("IL12IP",3))

mat_col <- data.frame(group = col_groups)
rownames(mat_col) <- colnames(k_sel_s)
#mat_colors <- list(group = brewer.pal(4, "Set1"))
mat_colors<-list(group=c("red","green","blue"))
#names(mat_colors$group) <- unique(col_groups)
names(mat_colors$group)<-c("PBS", "IL12IV","IL12IP")

colores<-diverge_hcl(6) #blue, white,red
mat_breaks <- seq(min(k_sel_s), max(k_sel_s), length.out = 6)

# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(k_sel_s), 0, length.out=ceiling(length(colores)/2) + 1), 
              seq(max(k_sel_s)/length(colores), max(k_sel_s), length.out=floor(length(colores)/2)))

pdf("Heatmap_Genes_CDITRANI_KEGG_TNF_signaling_pathway.pdf")

pheatmap(k_sel_s,
    color=colores,
    show_rownames     = TRUE,
    annotation_col    = mat_col,
    breaks            = myBreaks,
    annotation_colors = mat_colors,
    #border_color="black",
    #gaps_col=3,
    drop_levels       = TRUE,
    fontsize          = 10,
    main              = "Zscored Gene Expression \n 'TNF signaling pathway' enriched KEGG pathway \n IL12IP vs IL12IV ",
    cluster_cols = FALSE,
    border_color = TRUE,
    clustering_distance_rows="correlation",
    cellheight=12, 
    cellwidth = 12,

    angle_col=c("45"),
)
dev.off()
