args=(commandArgs(TRUE))


filequery<-args[1]
outputpath<-args[2]

SRA<-read.table(filequery)
SRA<-paste(SRA$V1)

#SRA<-c("SRR4423266","SRR4423267","SRR4423268","SRR4423269","SRR4423270","SRR4423271","SRR4423272","SRR4423273","SRR4423274","SRR4423275","SRR4423276","SRR4423277")


for(file in SRA){
    cat(file,"\n")
    y<-paste0("fastq-dump --gzip --split-files --outdir ",outputpath," ", file)
    write(y,file="launch_fastqdump.cmd",append=TRUE)
}

