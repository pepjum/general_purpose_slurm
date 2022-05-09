## flowSOM


library(FlowSOM)
#library(flowCT)
library(PeacoQC)
library(FlowCT)

metadata<-data.frame("filename"=c("Eomes 15_2bles.fcs", "Eomes 20_2bles.fcs","Eomes 3_2bles.fcs"),sample_id=c(1:3),condition=c("Augmented","depletioncd8","control"))

fcs_ALL <- fcs.SCE(directory = "/home/inmuno/Documents/" , pattern = "fcs", metadata=metadata, events = 50000, transf.cofactor = 500, project.name = "Maite")

# raw and transformed

#metadata<-data.frame("filename"=c("Eomes 15_2bles.fcs", "Eomes 3_2bles.fcs","Eomes 20_2bles.fcs"),sample_id=c(1:3),condition=c("Augmented","control","depletioncd8"))

#fcs_ALL <- fcs.SCE(directory = ".", pattern = "fcs", metadata=metadata, events = 10000, transf.cofactor = 500, project.name = "Maite")

# marker delimitation/selection for further analysis

# marker.names(fcs_ALL)

# surface_markers <- c("FL3-A:B690-PC5.5-A", "FL4-A:Y585-PE-A", "FL5-A:Y610-mCHERRY-A","FL13-A:V525-KrO-A","FL14-A:V610-A","FL15-A:V660-A","FL17-A:UV405-A","FL18-A:UV525-A")

# physical_markers <- c("FSC-A:FSC-A","SSC-A:SSC-A","SSC-H:SSC-H","FSC-H:FSC-H")

# # qc and remove doublets

# fcs_rem <- qc.and.removeDoublets(fcs_ALL, physical.markers = physical_markers, return.fcs = F)  


# # normalizacion marcadores  

# fcs_n <- marker.normalization(fcs.SCE = fcs_rem, marker = c("FL3-A:B690-PC5.5-A", "FL4-A:Y585-PE-A", "FL5-A:Y610-mCHERRY-A","FL13-A:V525-KrO-A","FL14-A:V610-A","FL15-A:V660-A","FL17-A:UV405-A","FL18-A:UV525-A"), method = "gauss")

# # tranformed matrix
# # eliminacion de efecto batch




# fcs <- batch.removal(fcs_n, method = "harmony", batch = "condition", new.matrix.name = "batch")

# pdf("batch.pdf")
# multidensity(fcs.SCE = fcs, assay.i = "raw", show.markers=c("FL3-A:B690-PC5.5-A","FL4-A:Y585-PE-A", "FL5-A:Y610-mCHERRY-A", "FL13-A:V525-KrO-A", "FL14-A:V610-A", "FL15-A:V660-A", "FL17-A:UV405-A","FL18-A:UV525-A"  ), subsampling = 100, ridgeline.lim = 0)
# multidensity(fcs.SCE = fcs, assay.i = "transformed", show.markers=c("FL3-A:B690-PC5.5-A","FL4-A:Y585-PE-A", "FL5-A:Y610-mCHERRY-A", "FL13-A:V525-KrO-A", "FL14-A:V610-A", "FL15-A:V660-A", "FL17-A:UV405-A","FL18-A:UV525-A"  ), subsampling = 100, ridgeline.lim = 0)
# multidensity(fcs.SCE = fcs, assay.i = "batch", show.markers=c("FL3-A:B690-PC5.5-A","FL4-A:Y585-PE-A", "FL5-A:Y610-mCHERRY-A", "FL13-A:V525-KrO-A", "FL14-A:V610-A", "FL15-A:V660-A", "FL17-A:UV405-A","FL18-A:UV525-A"  ), subsampling = 100, ridgeline.lim = 0)

# dev.off()

# ###
# a<-assays(fcs)[["batch"]]
# a_transposed<-t(a)
# a_transposed<-as.data.frame(a_transposed)
# b<-paste(colData(fcs)[["filename"]])
# b_df<-data.frame("sample"=b)

# a_transposed<-cbind(a_transposed,b_df)

# # system("mkdir -p temporary")
# # setwd(paste0(getwd(),"/temporary"))
# selection<-data.frame()

# set.seed(123)

#  for (sample in unique(a_transposed$sample)){
#      df_sub<-a_transposed[which(a_transposed$sample==sample),]
#      name_file<-paste0(df_sub$sample[1],"_temp")
#      df_temp<-df_sub[,c(1:ncol(df_sub)-1)]

#      x <- sample(1:nrow(df_temp), size = 10000)
#      df_temp_s<-df_temp[x,]
#      selection<-rbind(selection,df_temp_s)   
#      #df_temp_flow<-flowAssist::DFtoFF(df_temp)
#      #flowCore::write.FCS(df_temp_flow, name_file, what="numeric", delimiter = "|", endian="big")
#  }



# head(selection)

# selection_sel<-selection[,c("FL3-A:B690-PC5.5-A", "FL4-A:Y585-PE-A","FL5-A:Y610-mCHERRY-A","FL13-A:V525-KrO-A","FL14-A:V610-A","FL15-A:V660-A","FL17-A:UV405-A","FL18-A:UV525-A")]

# selection_ff<-flowAssist::DFtoFF(selection_sel)


# #ff <- flowCore::read.FCS("/home/inmuno/Documents/Eomes 20_2bles.fcs")
# files<-list.files(".", pattern="_temp")
# ff<-flowCore::read.flowSet(files)

# dimensions<-c()
# for(i in 1:length(files)){
#     k<-nrow(ff[[i]])
#     dimensions<-c(dimensions, k)
# }


markers_of_interest <- c("FL3-A","FL4-A" ,"FL5-A", "FL13-A",
 "FL14-A", "FL15-A", "FL17-A" , "FL18-A" )

channels_of_interest <- markers_of_interest

## REmove margins va antes de compensar


files<-list.files("/home/inmuno/Documents/", pattern=".fcs")
files<-paste0("/home/inmuno/Documents/", files)

for (file in files){
    namefilet<-basename(file)
    namefile<-paste(lapply(strsplit(paste(namefilet),"\\."),"[",1))
    cat(namefile,"\n")
    ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)    
    ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest) # more than 76 % is marginal
# # Compensation
    comp <- flowCore::keyword(ff_m)[["SPILL"]]
    ff_c <- flowWorkspace::compensate(ff_m, comp)
#  # Transformation
    transformList <- flowCore::estimateLogicle(ff_c, channels = colnames(comp))
    ff_t <- flowWorkspace::transform(ff_c, transformList)
# #### Remove doublets 
    ff_clean<-RemoveDoublets(ff_t, channel1="FSC-A", channel2="FSC-H", nmad=4, verbose=FALSE, output="frame")
# subsampling
    set.seed(42)
    x <- sample(1:nrow(ff_clean), size = 50000)
    FF<-ff_clean[x,]

    SOM_x <- 12
    SOM_y <- 12
    nClust = 12
    seed <- 2020
    scaling <- FALSE
    #build the model
    fsom <- FlowSOM(input = FF,scale = scaling, nClus=nClust, colsToUse = markers_of_interest, seed = seed, xdim = SOM_x, ydim = SOM_y)
    #output
    pepeFlowSOMmary(fsom = fsom, n_tsne=20000 ,plotFile = paste0("/home/inmuno/Documents/output_fsom/",namefile,"FSOM_output.pdf"))

}

PlotStars(fsom = fsom,
 backgroundValues = fsom$metaclustering)

ggsave(paste0(dir_results, "fsom_tree.pdf"), height = 8.5, width = 11)


#informativo
str(fSOMB$map, max.level = 2)

#third step
fSOM3 <- BuildMST(fSOMB)

#Plot the MST
pdf("MST_sub.pdf")
PlotStars(fSOM3, view="grid")
dev.off()


pdf("Plot_numbers.pdf")
PlotNumbers(fSOM3)
dev.off()



pdf("MST2.pdf")
PlotStars(fSOM3, equalNodeSize = TRUE)
dev.off()


metaClustering <- as.character(metaClustering_consensus(fSOM3$map$codes,k = 6))

pdf("plotLabels.pdf")
PlotLabels(fSOM3, labels=metaClustering)
dev.off()

#extraer informacion de cada nodo

metaClustering_perCell <- GetMetaclusters(fSOM3, metaClustering)
table(metaClustering_perCell)

### Detecting nodes with a specific pattern

# NK query
 query <- c( "FL5-A"="high", "FL15-A"="low", "FL14-A"="low","FL3-A"="low","FL18-A"="low","FL17-A"="low","FL4-A"="high") # fl4-a =high 
 query <- c( "FL5-A"="high", "FL15-A"="low", "FL14-A"="low","FL3-A"="low","FL18-A"="low","FL17-A"="low","FL4-A"="high") # fl4-a =low 

query2<-c("FL3-A"="high","FL17-A"="high") # cd8
query3<-c("FL3-A"="high","FL18-A"="high", "FL17-A"="low","FL15-A"="low","FL14-A"="low","FL5-A"="low")

query4<-c("FL15-A:V660-A"="high","FL3-A:B690-PC5.5-A"="low")

query_res <- QueryStarPlot(fsom, query4, equalNodeSize = TRUE, plot = FALSE)
cellTypes <- factor(rep("Unlabeled", fsom$map$nNodes),levels=c("Unlabeled", "Bcells"))
cellTypes[query_res$selected] <- "Bcells"

p <- PlotStars(fsom, backgroundValues=cellTypes, backgroundColor=c("#FFFFFF00","#0000FF"))

pdf("Bcells.pdf")
plot(p)
dev.off()

pepeFlowSOMmary<-function (fsom, n_tsne, plotFile) 
{
    fsom <- UpdateFlowSOM(fsom)
    metaclustersPresent <- !is.null(fsom$metaclustering)
    if (metaclustersPresent) {
        mfis <- GetMetaclusterMFIs(fsom, colsUsed = TRUE)
        metaclusters <- GetMetaclusters(fsom)
        nMetaclusters <- NMetaclusters(fsom)
    }
    clusters <- GetClusters(fsom)
    nNodes <- seq_len(NClusters(fsom))
    filePresent <- "File" %in% colnames(fsom$data)
    plotList <- list()
    message("Plot FlowSOM trees")
    for (view in c("MST", "grid")) {
        plotList[[paste0("stars_", view)]] <- PlotStars(fsom, 
            view = view, backgroundValues = fsom$metaclustering, 
            title = paste0("FlowSOM ", view))
        if (metaclustersPresent) {
            p2.1 <- PlotFlowSOM(fsom, view = view, title = "FlowSOM Clusters", 
                equalNodeSize = TRUE) %>% AddNodes(values = fsom$metaclustering, 
                showLegend = TRUE, label = "Metaclusters") %>% 
                AddLabels(labels = seq_len(NClusters(fsom)))
            p2.2 <- PlotFlowSOM(fsom, view = view, equalNodeSize = TRUE, 
                title = "FlowSOM Metaclusters") %>% AddNodes(values = fsom$metaclustering, 
                showLegend = TRUE, label = "Metaclusters") %>% 
                AddLabels(labels = as.numeric(fsom$metaclustering))
        }
        else {
            p2.1 <- PlotNumbers(fsom, view = view, title = "FlowSOM Clusters", 
                maxNodeSize = "auto", equalNodeSize = TRUE)
            p2.2 <- NULL
        }
        plotList[[paste0("labels_", view)]] <- ggpubr::ggarrange(p2.1, 
            p2.2, common.legend = TRUE, legend = "right")
    }

    plotList[["p5"]] <- PlotMarker(fsom, marker = fsom$map$colsUsed, 
        refMarkers = fsom$map$colsUsed, equalNodeSize = TRUE)
    if (filePresent) {
        message("Plot file distribution")
        p6 <- PlotPies(fsom, cellTypes = factor(fsom$data[, "File"]), 
            equalNodeSize = TRUE, view = "grid", title = "File distribution per cluster", 
            colorPalette = FlowSOM_colors)
        filesI <- as.character(unique(fsom$data[, "File"]))
        expectedDistr <- rep(1, length(filesI))
        names(expectedDistr) <- filesI
        arcsDf <- ParseArcs(0, 0, expectedDistr, 0.7)
        arcsDf$Markers <- factor(arcsDf$Markers, levels = filesI)
        plotList[["p6"]] <- AddStarsPies(p6, arcsDf, colorPalette = FlowSOM_colors(length(filesI) + 
            1), showLegend = FALSE) + ggplot2::geom_text(ggplot2::aes(x = 0, 
            y = -1, label = "Expected distribution"))
    }
    message("Calculate t-SNE")
    dimred_res <- PlotDimRed(fsom, cTotal = n_tsne, colorBy = fsom$map$colsUsed, 
        seed = 1, returnLayout = TRUE, title = paste0("t-SNE with markers used in FlowSOM ", 
            "call (perplexity = 30, cells =", n_tsne,")"))
    if (metaclustersPresent) {
        plotList[["p7"]] <- PlotDimRed(fsom, dimred = dimred_res$layout, 
            seed = 1, title = paste0("t-SNE with markers used in FlowSOM Global ", 
                "call (perplexity = 30, cells =", n_tsne,")"))    
    }

    query4<-c("FL15-A"="high")

    query_res <- QueryStarPlot(fsom, query4, equalNodeSize = TRUE, plot = FALSE)
    cellTypes <- factor(rep("Unlabeled", fsom$map$nNodes),levels=c("Unlabeled", "Bcells"))
    cellTypes[query_res$selected] <- "Bcells"

    #p <- PlotStars(fsom, backgroundValues=cellTypes, backgroundColor=c("#FFFFFF00","#0000FF"))
    plotList[["p16"]]<-PlotStars(fsom, backgroundValues=cellTypes, backgroundColor=c("#FFFFFF00","#0000FF"))
    
    if (metaclustersPresent) {
        message("Plot cluster per metacluster distibution")
        clusterPerMetacluster <- data.frame(metaclusters, clusters) %>% 
            dplyr::group_by(metaclusters) %>% dplyr::count(clusters) %>% 
            dplyr::mutate(Percentage = .data$n/sum(.data$n) * 
                100) %>% dplyr::mutate(LabelPos = cumsum(.data$Percentage) - 
            1) %>% as.data.frame()
        clusterPerMetacluster$clusters <- factor(clusterPerMetacluster$cluster)
        plotList[["p9"]] <- ggplot2::ggplot(clusterPerMetacluster, 
            ggplot2::aes(x = metaclusters)) + ggplot2::geom_bar(ggplot2::aes(y = .data$Percentage, 
            fill = metaclusters), stat = "identity", col = "black") + 
            ggplot2::geom_text(ggplot2::aes(y = .data$LabelPos, 
                label = clusters)) + ggplot2::theme_minimal() + 
            ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(), 
                panel.grid.minor.x = ggplot2::element_blank(), 
                legend.position = "none") + ggplot2::ggtitle("Percentages clusters per metacluster")
        message("Plot heatmap")
        colnames(mfis) <- fsom$prettyColnames[colnames(mfis)]
        rownames(mfis) <- levels(metaclusters)
        mfis_scaled <- scale(mfis)
        mfis_scaled[is.na(mfis_scaled)] <- 0
        plotList[["empty"]] <- ggplot2::ggplot() + ggplot2::theme_minimal()
        plotList[["p10"]] <- pheatmap::pheatmap(mfis_scaled, 
            scale = "none", display_numbers = round(mfis, 2), 
            main = "Median expression per metacluster", silent = TRUE)
    }
    message("Make tables")
    datatable1 <- data.frame(`Total number of modelated cells` = nrow(fsom$data), 
        check.names = FALSE)
    rownames(datatable1) <- "Summary FlowSOM \n modified by José González Gomariz, PhD"
    if (metaclustersPresent) {
        datatable1[, "Total number of metaclusters"] <- nMetaclusters
    }
    markersInFlowSOM <- split(fsom$prettyColnames[fsom$map$colsUsed], 
        rep(seq_len(ceiling(length(fsom$prettyColnames[fsom$map$colsUsed])/5)), 
            each = 5)[seq_len(length(fsom$prettyColnames[fsom$map$colsUsed]))]) %>% 
        sapply(paste, collapse = ", ") %>% paste(collapse = ",\n")
    datatable1 <- cbind(datatable1, `Total number of clusters` = fsom$map$nNodes, 
        `Markers used for FlowSOM` = markersInFlowSOM)
    datatable1 <- format(datatable1, digits = 2)
    t1 <- ggpubr::ggtexttable(t(datatable1), theme = ggpubr::ttheme("minimal"))
    if (metaclustersPresent) {
        freqMetaclusters <- data.frame(metaclusters = as.character(metaclusters)) %>% 
            dplyr::count(.data$metaclusters) %>% dplyr::mutate(percentage = .data$n/sum(.data$n) * 
            100) %>% as.data.frame()
        datatable2 <- data.frame(Metacluster = levels(metaclusters), 
            `Number of cells` = 0, `Percentage of cells` = 0, 
            `Number of clusters` = 0, Clusters = "", check.names = FALSE)
        rownames(datatable2) <- datatable2$Metacluster
        datatable2[freqMetaclusters[, "metaclusters"], "Number of cells"] <- freqMetaclusters[, 
            "n"]
        datatable2[freqMetaclusters[, "metaclusters"], "Percentage of cells"] <- freqMetaclusters[, 
            "percentage"]
        datatable2[, "Number of clusters"] <- sapply(levels(metaclusters), 
            function(x) {
                which(fsom$metaclustering == x) %>% length()
            })
        datatable2[, "Clusters"] <- sapply(levels(metaclusters), 
            function(x) {
                cl_selected <- which(fsom$metaclustering == x)
                clusters_string <- split(cl_selected, rep(seq_len(ceiling(length(cl_selected)/25)), 
                  each = 25)[seq_len(length(cl_selected))]) %>% 
                  sapply(paste, collapse = ", ") %>% paste(collapse = ",\n")
                return(clusters_string)
            })
        datatable2 <- format(datatable2, digits = 2)
        split_datatable2 <- split(datatable2, rep(seq_len(ceiling(nrow(datatable2)/30)), 
            each = 30, length.out = nrow(datatable2)))
    }
    freqClusters <- as.data.frame(clusters) %>% dplyr::count(.data$clusters) %>% 
        dplyr::mutate(freq = .data$n/sum(.data$n) * 100)
    datatable3 <- data.frame(Cluster = factor(nNodes), `Number of cells` = 0, 
        `Percentage of cells` = 0, check.names = FALSE)
    datatable3[freqClusters[, "clusters"], "Number of cells"] <- freqClusters[, 
        "n"]
    datatable3[freqClusters[, "clusters"], "Percentage of cells"] <- freqClusters[, 
        "freq"]
    if (metaclustersPresent) {
        datatable3 <- cbind(datatable3, `Belongs to metacluster` = fsom$metaclustering, 
            `Percentage in metacluster` = 0)
        datatable3[as.numeric(as.character(clusterPerMetacluster[, 
            "clusters"])), "Percentage in metacluster"] <- clusterPerMetacluster[, 
            "Percentage"]
    }
    datatable3 <- format(datatable3, digits = 2)
    split_datatable3 <- split(datatable3, rep(seq_len(ceiling(nrow(datatable3)/30)), 
        each = 30, length.out = nrow(datatable3)))
    message("Printing")
    if (!is.null(plotFile)) {
        grDevices::pdf(plotFile, width = 17, height = 10)
        print(t1)
        if (metaclustersPresent) {
            for (table2 in split_datatable2) {
                print(ggpubr::ggtexttable(table2, theme = ggpubr::ttheme("minimal"), 
                  rows = NULL))
            }
        }
        for (table3 in split_datatable3) {
            print(ggpubr::ggtexttable(table3, theme = ggpubr::ttheme("minimal"), 
                rows = NULL))
        }
        for (plot in plotList) {
            print(plot)
        }
        dev.off()
    }
    else {
        return(plotList)
    }
}


pepeFlowSOMmary(fsom=fsom, n_tsne=30000, plotfile)

