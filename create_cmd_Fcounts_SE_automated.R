# ireate unch_fastqc.cmd",append=TRUE)cmd file for download SRA
args=(commandArgs(TRUE))

pathfiles<-args[1]
gtffile<-args[2]
outputpath<-args[3]

Bam_files<-list.files(pathfiles, pattern="out.bam")
Bamfiles<-paste(lapply(strsplit(paste(Bam_files),"\\."),"[",1))


for (i in 1:length(Bam_files)){
    y<-paste0("featureCounts -t exon -g gene_id -a ", gftfile, " -o ", paste0(outputpath,"/", Bamfiles[i],".txt")," ", paste0(pathfiles,"/",Bam_files[i]) )
    write(y,file="launch_fCounts.cmd",append=TRUE)
	cat(y, "\n")
}



