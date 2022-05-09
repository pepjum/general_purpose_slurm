# ireate unch_fastqc.cmd",append=TRUE)cmd file for download SRA


Bam_files<-list.files("/home/jgonzalezgom/01_NMELERO/GSM2344976/results_STAR/", pattern="out.bam")
Bamfiles<-paste(lapply(strsplit(paste(Bam_files),"\\."),"[",1))


for (i in 1:length(Bam_files)){
    y<-paste0("featureCounts  -t exon -g gene_id -a /home/jgonzalezgom/references/gencode.v38.annotation.gtf -o ", paste0("/home/jgonzalezgom/01_NMELERO/GSM2344976/results_fC/", Bamfiles[i],".txt")," ", paste0("/home/jgonzalezgom/01_NMELERO/GSM2344976/results_STAR/",Bam_files[i]) )
    write(y,file="launch_fCounts.cmd",append=TRUE)
	cat(y, "\n")
}



