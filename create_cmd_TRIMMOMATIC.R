# create cmd file for download SRA


#SRA<-read.table("SRR_Acc_List.txt")
#SRA<-paste(SRA$V1)

files<-list.files("/home/jgonzalezgom/03_MALVAREZ/fastq/", pattern=".gz")
files<-paste0("/home/jgonzalezgom/03_MALVAREZ/fastq/", files)
SAMPLES_tmp<-basename(files)
SAMPLES_tmp<-basename(files)
SAMPLES_tmp2<-paste(lapply(strsplit(paste(SAMPLES_tmp),".fastq"),"[",1))

SAMPLES<-gsub('.{1}$', '', SAMPLES_tmp2)




for(sample in SAMPLES){
    cat(sample,"\n")
    y<-paste0("java -jar /beegfs/easybuild/common/software/Trimmomatic/0.38-Java-1.8/trimmomatic-0.38.jar PE -phred33 ", paste0("/home/jgonzalezgom/03_MALVAREZ/fastq/",sample,"1.fastq.gz "), paste0("/home/jgonzalezgom/03_MALVAREZ/fastq/",sample,"2.fastq.gz "), paste0("/home/jgonzalezgom/03_MALVAREZ/fastq_trim/",sample,"1_T.fastq.gz "), paste0("/home/jgonzalezgom/03_MALVAREZ/fastq_trim/",sample,"1_T_Unpaired.fastq.gz "), paste0("/home/jgonzalezgom/03_MALVAREZ/fastq_trim/",sample,"2_T.fastq.gz "), paste0("/home/jgonzalezgom/03_MALVAREZ/fastq_trim/",sample,"2_T_Unpaired.fastq.gz "), "ILLUMINACLIP:/beegfs/easybuild/common/software/Trimmomatic/0.38-Java-1.8/adapters/TruSeq3-PE-2.fa/:2:30:10 SLIDINGWINDOW:15:25 MINLEN:36")
    write(y,file="launch_trimmomatic.cmd",append=TRUE)
}

