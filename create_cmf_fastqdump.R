# create cmd file for download SRA


SRA<-read.table("SRR_Acc_List.txt")
SRA<-paste(SRA$V1)

for(file in SRA2){
    cat(file,"\n")
    y<-paste0("fastq-dump --gzip --split-files --outdir /home/jgonzalezgom/01_NMELERO/ERP105482/ ", file)
    write(y,file="launch_fastqdump.cmd",append=TRUE)
}

