# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-list.files("~/01_NMELERO/ERP105482/", pattern=".gz")
fileName<-paste(lapply(strsplit(paste(files),"_"),"[", 1))
files<-paste0("/home/jgonzalezgom/01_NMELERO/ERP105482/", fileName)



for(i in 1:length(files)){
    cat(files[i],"\n")
    #y<-paste0("fastqc --outdir /home/jgonzalezgom/01_NMELERO/ERP105482/fastqc/ ", file)
	y<-paste0("STAR --genomeDir /home/jgonzalezgom/references/STAR_hg38_ref2/hg38_index --readFilesCommand zcat --runThreadN 6 --readFilesIn ", paste0(files[i],"_1.fastq.gz")," ",paste0(files[i],"_2.fastq.gz"), " --outFileNamePrefix ", paste0("/home/jgonzalezgom/01_NMELERO/ERP105482/results_STAR/",fileName[i])," --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard")
    write(y,file="launch_STAR.cmd",append=TRUE)
}

