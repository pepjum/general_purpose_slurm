# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-list.files("/home/jgonzalezgom/04_IOLIVERA/fastq/", pattern="T.fastq.gz")
files<-paste0("/home/jgonzalezgom/04_IOLIVERA/fastq/", files)

for(file in files){
    cat(file,"\n")
    y<-paste0("fastqc --outdir /home/jgonzalezgom/04_IOLIVERA/trimmed/ ", file)
    write(y,file="launch_fastqc.cmd",append=TRUE)
}

