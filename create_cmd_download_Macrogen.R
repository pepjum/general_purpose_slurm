# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-read.table("/home/jgonzalezgom/04_IOLIVERA/download_macrogen_irene.txt")
files<-paste(files$V1)


for(file in files){
    cat(file,"\n")
    y<-paste0("wget -P /home/jgonzalezgom/04_IOLIVERA/fastq/ ", file)
    write(y,file="launch_download_macrogen.cmd",append=TRUE)
}

