# ireate unch_fastqc.cmd",append=TRUE)cmd file for download SRA


Bam_files<-list.files("/home/jgonzalezgom/01_NMELERO/ERP105482/results_STAR/", pattern="sorted.bam")
Bam_files<-paste0("/home/jgonzalezgom/01_NMELERO/ERP105482/results_STAR/", Bam_files)


#Bamfiles<-paste(lapply(strsplit(paste(Bam_files),"\\."),"[",1))


for (i in 1:length(Bam_files)){
    y<-paste0("samtools index ", Bam_files[i])
    write(y,file="launch_samtools_index.cmd",append=TRUE)
	cat(y, "\n")
}



