# create cmd file for download SRA
args=(commandArgs(TRUE))

pathfiles<-args[1]
patternfile<-args[2]
outputpath<-args[3]

files<-list.files(pathfiles,  pattern=patternfile)
files<-paste0(pathfiles,"/", files)

for(file in files){
    cat(file,"\n")
    y<-paste0("fastqc --outdir ", outputpath," ", file)
    write(y,file="launch_fastqc.cmd",append=TRUE)
}

