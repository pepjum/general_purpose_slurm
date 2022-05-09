# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-list.files("~/03_MALVAREZ/fastq_trim/", pattern="T.fastq.gz")
fileName<-paste(lapply(strsplit(paste(files),"_"),"[", 1))
fileName<-unique(fileName)
files<-paste0("/home/jgonzalezgom/03_MALVAREZ/fastq_trim/", fileName)



for(i in 1:length(fileName)){
    cat(fileName[i],"\n")
    #y<-paste0("fastqc --outdir /home/jgonzalezgom/01_NMELERO/ERP105482/fastqc/ ", file)
	y<-paste0("STAR --genomeDir /home/jgonzalezgom/references/STAR_mouse_GHRC39/ --readFilesCommand zcat --runThreadN 6 --readFilesIn ", paste0(files[i],"_1_T.fastq.gz")," ",paste0(files[i],"_2_T.fastq.gz"), " --outFileNamePrefix ", paste0("/home/jgonzalezgom/03_MALVAREZ/results_STAR/",fileName[i])," --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard")
    write(y,file="launch_STAR.cmd",append=TRUE)
}

