args=(commandArgs(TRUE))

pathfile<-args[1]
patternfile<-args[2]
pathgenome<-args[3]
outputpath<-args[4]

files<-list.files(pathfile, pattern=patternfile)

fileName<-paste(lapply(strsplit(paste(files),"_T"),"[", 1))
fileName<-paste(lapply(strsplit(paste(fileName),"_R"),"[", 1))

fileName<-unique(fileName)
files<-paste0(pathfile,"/", fileName)



for(i in 1:length(fileName)){
    
	cat(fileName[i],"\n")
    
	#y<-paste0("STAR --genomeDir ",pathgenome, "--readFilesCommand zcat --runThreadN 6 --readFilesIn ", paste0(files[i],"_R1_",patternfile)," ",paste0(files[i],"_R2_",patternfile), " --outFileNamePrefix ", paste0(outputpath,"/",fileName[i])," --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard")
	y<-paste0("STAR --genomeDir ",pathgenome, " --readFilesCommand zcat --runThreadN 6 --readFilesIn ", paste0(files[i],"_R1",patternfile)," ",paste0(files[i],"_R2",patternfile), " --outFileNamePrefix ", paste0(outputpath,"/",fileName[i])," --outSAMtype BAM SortedByCoordinate --outSAMattributes Standard")
	cat(y,"\n")

    write(y,file="/home/jgonzalezgom/scripts_general_purpose/launch_STAR.cmd",append=TRUE)
}

