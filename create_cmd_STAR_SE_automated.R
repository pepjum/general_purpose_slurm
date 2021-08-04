# create cmd file for download SRA
args=(commandArgs(TRUE))

pathfile<-args[1]
patternfile<-args[2]
pathgenome<-args[3]
outputpath<-args[4]

files<-list.files(pathfile, pattern=patternfile)
fileName<-paste(lapply(strsplit(paste(files),"_"),"[", 1))
fileName<-unique(fileName)
files<-paste0(pathfile,"/", fileName)



for(i in 1:length(fileName)){

    cat(fileName[i],"\n")
 
	y<-paste0("STAR --genomeDir ", pathgenome," --readFilesCommand zcat --runThreadN 6 --readFilesIn ", paste0(files[i],"_1.fastq.gz")," ", " --outFileNamePrefix ", paste0(outputpath,"/",fileName[i])," --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard")
    write(y,file="launch_STAR.cmd",append=TRUE)
}

