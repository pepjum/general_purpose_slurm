# create cmd file for download SRA
args=(commandArgs(TRUE))

pathfiles<-args[1]
patternfile<-args[2]
outputpath<-args[3]
adapters<-args[4]
slidingwindow<-args[5]
quality<-args[6]
minlen<-args[7]



files<-list.files(pathfile, pattern=patternfile)
files<-paste0(pathfile,"/", files)
SAMPLES_tmp<-basename(files)
SAMPLES_tmp<-basename(files)
SAMPLES_tmp2<-paste(lapply(strsplit(paste(SAMPLES_tmp),".fastq"),"[",1))

SAMPLES<-gsub('.{1}$', '', SAMPLES_tmp2)




for(sample in SAMPLES){
    cat(sample,"\n")
    y<-paste0("java -jar /beegfs/easybuild/common/software/Trimmomatic/0.38-Java-1.8/trimmomatic-0.38.jar PE -phred33 ", paste0(pathfile,"/",sample,"1",patternfile), paste0(pathfile,"/",sample,"2", patternfile), paste0(outputpath,"/",sample,"1_T.fastq.gz "), paste0(outputpath,"/",sample,"1_T_Unpaired.fastq.gz "), paste0(outputpath,"/",sample,"2_T.fastq.gz "), paste0(outputpath,"/",sample,"2_T_Unpaired.fastq.gz "), "ILLUMINACLIP:",adapters,":2:30:10 SLIDINGWINDOW:",slidingwindow,":",quality ," MINLEN:",minlen)
    write(y,file="launch_trimmomatic.cmd",append=TRUE)
}

