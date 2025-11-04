extract_mid_sequences_1 <- function(i,pos1,pos2,seqs){
  midseqs=c()

  s=as.character(seqs[[i]])
  for (j in seq_along(pos1[[i]])) {
    start1 <- start(pos1[[i]])[j]
    end1 <- end(pos1[[i]])[j]
    start2 <- start(pos2[[i]])[j]
    end2 <- end(pos2[[i]])[j]
    if(!is.na(start1[1]) && !is.na(start2[1])){
      midseq=substr(s,end1+1,start2-1)
      if(nchar(midseq)==24){
        midseqs=append(midseqs,midseq)
      }
    }
  }
  return(midseqs)
}
extract_mid_sequences_2 <- function(i,pos1_rc,pos2_rc,seqs){
  midseqs=c()
  
  s=as.character(seqs[[i]])
  for (j in seq_along(pos1_rc[[i]])) {
    start1 <- start(pos1_rc[[i]])[j]
    end1 <- end(pos1_rc[[i]])[j]
    start2 <- start(pos2_rc[[i]])[j]
    end2 <- end(pos2_rc[[i]])[j]
    if(!is.na(start1[1]) && !is.na(start2[1])){
      midseq=substr(s,end2+1,start1-1)
      if(nchar(midseq)==24){
        midseqs=append(midseqs,midseq)
      }
    }
  }
  return(midseqs)
}


setwd("D:/SZY/100genebarcode/downloads/250514-A00599A")


t=Sys.time()
library(ShortRead)
library(Biostrings)
known_seq1 <- DNAString("TGCAAGAGACTTCCATCCAG")
known_seq2 <- DNAString("gcttccggtctggttcgctttgaagctcga")
# 计算已知序列的反向互补序列
known_seq1_rc <- reverseComplement(known_seq1)
known_seq2_rc <- reverseComplement(known_seq2)

#####100gbnrcodePCR-LHE4661-W1_L1_1.fq##+++++++++++
###
# 读取FASTQ文件
file_name="100gbnrcodePCR-LHE4661-W1_L1_1"

fq <- readFastq(paste(c(file_name,".fq"),collapse = ""))
sequences <- sread(fq)

positions1 <- vmatchPattern(known_seq1, sequences,max.mismatch = 2)
positions2 <- vmatchPattern(known_seq2, sequences,max.mismatch = 2)
positions1_rc <- vmatchPattern(known_seq1_rc, sequences,max.mismatch = 2)
positions2_rc <- vmatchPattern(known_seq2_rc, sequences,max.mismatch = 2)
l=length(positions1)
rev_comp<-function(dna){
  return(
    paste(rev(sapply(unlist(strsplit(dna,split="")), function(x) {switch(x, A = "T", T = "A", G = "C", C = "G") })),collapse = ""))
}

results=c()
for (i in :l) {
  
  A=extract_mid_sequences_1(i,positions1,positions2, sequences)
  
  B=extract_mid_sequences_2(i,positions1_rc,positions2_rc, sequences)
  if(is.null(B)==F){B=rev_comp(B)}
  
  
  results=c(results,A,B)
  if(i%%10000==0){cat("Complete:",i/l*100,"%","\n")}
  
}
write.csv(results,paste(c(file_name,".csv"),collapse = ""))



#####
setwd("D:/SZY/100genebarcode")
barcode=read.csv("barcode.txt",header=F)
recg_seq=read.csv("",header = F)
stat_barcode=c()
for (i in seq_along(barcode)) {
  stat_barcode[i]=length(grep(barcode[i],results))+length(grep(rev_comp(barcode[i]),results))
}