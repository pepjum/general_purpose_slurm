# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-list.files("/home/jgonzalezgom/01_NMELERO/RIBAS_CCELL_RNASEQ/fastq/", pattern=".fastq")
fileName<-paste(lapply(strsplit(paste(files),"_"),"[", 1))
fileName<-unique(fileName)
files<-paste0("/home/jgonzalezgom/01_NMELERO/RIBAS_CCELL_RNASEQ/fastq/", fileName)



for(i in 1:length(fileName)){
    cat(fileName[i],"\n")
    #y<-paste0("fastqc --outdir /home/jgonzalezgom/01_NMELERO/ERP105482/fastqc/ ", file)
	y<-paste0("STAR --genomeDir /home/jgonzalezgom/references/STAR_human_hg38/ --runThreadN 6 --readFilesIn ", paste0(files[i],"_1.fastq")," ", "  , paste0("/home/jgonzalezgom/01_NMELERO/RIBAS_CCELL_RNASEQ/results_STAR/",fileName[i])," --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard")
    write(y,file="launch_STAR.cmd",append=TRUE)
}

