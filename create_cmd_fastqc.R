# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-list.files("~/01_NMELERO/ERP105482/", pattern=".gz")
files<-paste0("/home/jgonzalezgom/01_NMELERO/ERP105482/", files)

for(file in files){
    cat(file,"\n")
    y<-paste0("fastqc --outdir /home/jgonzalezgom/01_NMELERO/ERP105482/fastqc/ ", file)
    write(y,file="launch_fastqc.cmd",append=TRUE)
}

