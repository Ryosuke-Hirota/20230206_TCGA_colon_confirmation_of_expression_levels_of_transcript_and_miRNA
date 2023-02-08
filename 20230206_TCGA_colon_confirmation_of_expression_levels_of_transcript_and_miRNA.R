# This script is to draw plot about correlation between expression level of transcript and miRNA 
# 2023/02/06 made

# activate packages
library(stringr)

# import correspondence table between TCGA colon transcriptome bam and TCGA colon miRNA quantification file
# this table is located at "https://github.com/Ryosuke-Hirota/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data"
setwd("C:/Rdata/20221124_make_table_of_TCGA_transcriptome_bam_correspond_to_TCGA_miRNA_expression_data")
cor.table <-read.table("correspondence_table_between_TCGA_colon_transcriptome_bam_and_miRNA_qunatification_file.txt",sep="\t",header = T,stringsAsFactors = F)

# import list of transcripts that intersect with miRNAs in gencode v36
# this list is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata")
primir.list <-read.table("TCGA_hg38_transcript_intersect_with_miRNA.txt",sep="\t",header = F,stringsAsFactors = F)
primir.list[,2] <-primir.list[,2]-1
primir.list <-primir.list[primir.list[,2]!=primir.list[,8]&primir.list[,3]!=primir.list[,9],]

# remove duplicated gene names and make list regarding miRNA name and gene name
mir_transcipt <-primir.list[,c(4,11)]
mir_transcipt <-subset(mir_transcipt,!duplicated(mir_transcipt))
mir_transcipt <-mir_transcipt[order(mir_transcipt[,1]),]
primir.list <-primir.list[,c(4,11,10)]
transcripts <-unique(primir.list[,3])

# list TCGA colon transcript quantification files
# these files are located at "\\fsw-q02\okamura-lab\20221006_TCGA_colon_salmon_quant_transcriptome"
setwd("C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification")
transcript.quant <-list.files(path = "C:/Rdata/20221006_TCGA_colon_salmon_transcriptome_quantification",pattern = ".txt")
t.file.id <-gsub("_quant.txt","",transcript.quant)

# make table summarized TCGA colon transcript quantification
for (i in 1:length(transcript.quant)){
  transcript.file <-read.table(transcript.quant[i],sep = "\t",header = T,stringsAsFactors = F)
  t <-match(transcripts,transcript.file[,1])
  transcript.file <-transcript.file[t,]
  transcript.file <-transcript.file[,c(1,4)]
  colnames(transcript.file)[2] <-transcript.quant[i]
  if(i==1){
    transcript.quant.table <-transcript.file
  }else{
    transcript.quant.table <-merge(transcript.quant.table,transcript.file,by="Name")
  }}

# import table summarized miRNA quantification files
# this table is located at ""
setwd("C:/Rdata/20230105_TCGA_colon_miRNA_quantification")
miRNA.quant.table <-read.table("table_of_TCGA_colon_miRNA_quantifications.txt",sep="\t",header = T,stringsAsFactors = F,check.names = F)

# make new directory
setwd("C:/Rdata")
dir.create("20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")
setwd("C:/Rdata/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")

# make empty summary
sm <-as.data.frame(matrix(nrow = nrow(mir_transcipt),ncol = 10))
colnames(sm) <-c("miRNA","transcript","m_min","m_median","m_max","m_mean","t_min","t_median","t_max","t_mean")

# investigate quantile and mean of each expression level
for (i in 1:nrow(mir_transcipt)){
  # extract expression level of a certain miRNA 
  miRNA.df <-miRNA.quant.table[miRNA.quant.table[,1]==mir_transcipt[i,1],]
  miRNA.df <-as.data.frame(t(miRNA.df),stringsAsFactors = F)
  m.cor <-match(rownames(miRNA.df),cor.table[,5])
  miRNA.df[,2] <-rownames(miRNA.df)
  miRNA.df[,3] <-cor.table[m.cor,2]
  colnames(miRNA.df) <-c(mir_transcipt[i,1],"miRNA_file_name","transcript_file_name")
  rownames(miRNA.df) <-NULL
  miRNA.df <-miRNA.df[-1,]
  
  # extract expression level of a certain transcript
  transcript <-primir.list[primir.list[,1]==mir_transcipt[i,1]&primir.list[,2]==mir_transcipt[i,2],3]
  t <-match(transcript,transcript.quant.table[,1])
  transcript.df <-transcript.quant.table[t,]
  transcript.df <-transcript.df[,-1]
  transcript.df[1,] <-apply(transcript.df, 2, sum)
  transcript.df <-transcript.df[-2,] 
  transcript.df <-as.data.frame(t(transcript.df),stringsAsFactors = F)
  t.cor <-match(t.file.id,cor.table[,1])
  transcript.df[,2] <-cor.table[t.cor,2]
  colnames(transcript.df) <-c(mir_transcipt[i,2],"transcript_file_name")
  rownames(transcript.df) <-NULL
  
  # merge expression levels of miRNA and transcript
  mt.df <-merge(miRNA.df,transcript.df,by="transcript_file_name")
  mt.df <-subset(mt.df,!is.na(mt.df[,1]))
  mt.df <-mt.df[,c(4,2,1,3)]
  mt.df[,2] <-as.numeric(mt.df[,2])
  
  # calculate quantile and mean of miRNA expression
  m.quantile <-as.numeric(quantile(mt.df[,2]))
  m.mean <-mean(mt.df[,2])

  # calculate quantile and mean of transcript expression
  t.quantile <-as.numeric(quantile(mt.df[,1]))
  t.mean <-mean(mt.df[,1])
  
  # write summary
  sm[i,1] <-mir_transcipt[i,1]
  sm[i,2] <-mir_transcipt[i,2]
  sm[i,3] <-m.quantile[1]
  sm[i,4] <-m.quantile[3]
  sm[i,5] <-m.quantile[5]
  sm[i,6] <-m.mean
  sm[i,7] <-t.quantile[1]
  sm[i,8] <-t.quantile[3]
  sm[i,9] <-t.quantile[5]
  sm[i,10] <-t.mean
}

# import result of correlation analysis between expression level of miRNA and transcript
# this result is located at "https://github.com/Ryosuke-Hirota/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA"
setwd("C:/Rdata/20230110_TCGA_colon_transcriptome_bam_correlation_between_transcript_and_miRNA")
cor.result <-read.table("TCGA_colon_transcriptome_summary_of_correlation_between_transcript_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)

# merge a summary and result about correlation analysis
sm <-merge(cor.result,sm,by=c("miRNA","transcript"))

# output merged summary
setwd("C:/Rdata/20230206_TCGA_colon_confirmation_of_expression_levels_of_transcript_and_miRNA")
write.table(sm,"TCGA_colon_confirmation_of_expression_level_of_transcript_and_miRNA.txt",sep="\t",row.names = F,quote = F)

### draw volcano plot to show how many combinations have signigicant positive correlation
# remove NA row and annotate whether combinations have significant positive correlations
vp <-sm
vp <-vp[!is.na(vp[3]),]
vp[,14] <-NA
vp[vp[,3]>0&vp[,4]<0.05,14] <-"red"
vp[vp[,3]<=0|vp[,4]>=0.05,14] <-"grey40"

# draw valcano plot
pdf("valcano_plot_about_TCGA_colon_correlation_between_expression_level_of_transcript_and_miRNA.pdf")
plot(vp[,3],-log10(vp[,4]),xlab="correlation coefficient",ylab="-log10 (p.value)",col=vp[,14],pch=19,
     main = paste0("r>0 p<0.05 ",nrow(vp[vp[,14]=="red",]),"/",nrow(vp)," (sig.posi / total)"))
legend("topleft",legend =c("other","r>0, p<0.05"),col=unique(vp[,14]),pch=19)
abline(h=log10(0.05)*-1,v=0,lty=2)
dev.off()

### draw scatter plot to investigate whether combinations without correlation have small mean of miRNA/transcript expression 
# remove NA row and annotate whether combinations have significant positive correlations
sp <-sm
sp <-sp[!is.na(sp[3]),]
sp[,14] <-NA
sp[sp[,3]>0&sp[,4]<0.05,14] <-"red"
sp[sp[,3]<=0|sp[,4]>=0.05,14] <-"grey40"

# draw scatter plot
pdf("plot_about_mean_of_TCGA_colon_expression_level_of_miRNA_and_transcript.pdf")
plot(log2(sp[,9]),log2(sp[,13]),xlab="log2(mean of miRNA expression)",ylab="log2(mean of transcript expression)",pch=19,col=sp[,14],
     main = paste0("sig.posi.cor = ",nrow(sp[sp[,14]=="red",])," no.sig/sig.nega.cor = ",nrow(sp)-nrow(sp[sp[,14]=="red",])))
legend("topright",legend =c("other","r>0, p<0.05"),col=unique(sp[,14]),pch=19)
abline(h=0,v=0,lty=2)
dev.off()

# cutoff sample < 50
sp1 <-sp[sp[,5]>=50,]
pdf("plot_about_mean_of_TCGA_colon_expression_level_of_miRNA_and_transcript_cutoff_50.pdf")
plot(log2(sp1[,9]),log2(sp1[,13]),xlab="log2(mean of miRNA expression)",ylab="log2(mean of transcript expression)",pch=19,col=sp1[,14],
     main = paste0("sig.posi.cor = ",nrow(sp1[sp1[,14]=="red",])," no.sig/sig.nega.cor = ",nrow(sp1)-nrow(sp1[sp1[,14]=="red",])))
abline(h=0,v=0,lty=2)
dev.off()

sp2 <-sp[sp[,5]>=100,]
pdf("plot_about_mean_of_TCGA_colon_expression_level_of_miRNA_and_transcript_cutoff_100.pdf")
plot(log2(sp2[,9]),log2(sp2[,13]),xlab="log2(mean of miRNA expression)",ylab="log2(mean of transcript expression)",pch=19,col=sp2[,14],
     main = paste0("sig.posi.cor = ",nrow(sp2[sp2[,14]=="red",])," no.sig/sig.nega.cor = ",nrow(sp2)-nrow(sp2[sp2[,14]=="red",])))
abline(h=0,v=0,lty=2)
dev.off()
