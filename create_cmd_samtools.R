# ireate unch_fastqc.cmd",append=TRUE)cmd file for download SRA


Bam_files<-list.files("/home/jgonzalezgom/01_NMELERO/ERP105482/results_STAR/", pattern=".bam")
Bam_files<-paste0("/home/jgonzalezgom/01_NMELERO/ERP105482/results_STAR/", Bam_files)
Bamfiles<-paste(lapply(strsplit(paste(Bam_files),"\\."),"[",1))


for (i in 1:length(Bam_files)){
    y<-paste0("samtools sort ", Bam_files[i], " -o ", Bamfiles[i],"_sorted.bam")
    write(y,file="launch_samtools_sort.cmd",append=TRUE)
	cat(y, "\n")
}



