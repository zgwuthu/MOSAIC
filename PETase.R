library(ShortRead)
library(Biostrings)
rev_comp<-function(dna){
  return(
    paste(rev(sapply(unlist(strsplit(dna,split="")), function(x) {switch(x, A = "T", T = "A", G = "C", C = "G") })),collapse = ""))
}

#####concatenated_reads#####
setwd("D:/SZY/PETase_lib_FP_scaleup_250802/downloads/combined")
f_reads <- readFastq("PET-LHG19426-w1_L0_1.fq")
r_reads <- readFastq("PET-LHG19426-w1_L0_2.fq")
f_reads=sread(f_reads)
r_reads=sread(r_reads)


#writeFastq(concatenated_reads, "concatenated_reads.fq")

#####recognize sites#####

site1="GTTCTTATGGTGTAGTAAGT"#+1+3
site2="CCAATGCTACGGAGCGTTTT"
site3="GTGCTCTCGATTGGGCGTCA"

site4="TAGCACCTTGGCATACGACT"
site5="TTTTGGGTGGTCAAAATGAT"
site6="ATTGCCCCTGTCTCT"
site7="CACATGCAATTCCAATGTAT"
site8="GCCATAATTTTCCGAATTCG"

site4=rev_comp(site4)
site5=rev_comp(site5)
site6=rev_comp(site6)
site7=rev_comp(site7)
site8=rev_comp(site8)


l=length(f_reads)
var_aa=c()
var_AA=c()
for (i in 1:l) {
  f=as.character( f_reads[i] )
  r=as.character( r_reads[i] )
  pos1=matchPattern(site1,f)
  pos2=matchPattern(site2,f)
  pos3=matchPattern(site3,f)
  pos4=matchPattern(site4,r)
  pos5=matchPattern(site5,r)
  pos6=matchPattern(site6,r)
  pos7=matchPattern(site7,r)
  pos8=matchPattern(site8,r)

  check=!is.na(start(pos1)[1])&& !is.na(start(pos2)[1])&& !is.na(start(pos3)[1])&& !is.na(start(pos4)[1])&& !is.na(start(pos5)[1])&& !is.na(start(pos6)[1])&& !is.na(start(pos7)[1])&& !is.na(start(pos8)[1])
  
  if(check==F){f1=r;  r1=f
  pos1=matchPattern(site1,f1)
  pos2=matchPattern(site2,f1)
  pos3=matchPattern(site3,f1)
  pos4=matchPattern(site4,r1)
  pos5=matchPattern(site5,r1)
  pos6=matchPattern(site6,r1)
  pos7=matchPattern(site7,r1)
  pos8=matchPattern(site8,r1)
  check=!is.na(start(pos1)[1])&& !is.na(start(pos2)[1])&& !is.na(start(pos3)[1])&& !is.na(start(pos4)[1])&& !is.na(start(pos5)[1])&& !is.na(start(pos6)[1])&& !is.na(start(pos7)[1])&& !is.na(start(pos8)[1])
  if(check==T){f=f1;r=r1
  }
  }
  if(check){
    
    aa=paste(c(substr(f,end(pos1)[1]+1,end(pos1)[1]+3),
               substr(f,end(pos2)[1]+1,end(pos2)[1]+3),                substr(f,end(pos3)[1]+1,end(pos3)[1]+3),
rev_comp(paste(c(
  substr(r,start(pos8)[1]-3,start(pos8)[1]-1),
  substr(r,start(pos7)[1]-3,start(pos7)[1]-1),
  substr(r,start(pos6)[1]-3,start(pos6)[1]-1),
  substr(r,start(pos5)[1]-3,start(pos5)[1]-1),
  substr(r,start(pos4)[1]-3,start(pos4)[1]-1)
                                   ),collapse=""))
                  
    ),collapse = "")
    if(grepl("NULL", aa)==T){next}
    
    var_aa=append(var_aa,aa) 
    
    var_AA=append(var_AA,as.character( translate(DNAString(aa)) ))
  }
  if(i%%10000==0){cat("Complete:",i/l*100,"%","\n")}
}


u_var_AA=unique(var_AA)

t_var_AA=table(var_AA)


setwd("D:/SZY/PETase_lib_FP_scaleup_250802/PET_result2")

png("PET_var.jpg")
plot(sort(t_var_AA,decreasing = T)[1:1000],ylab ="count",ylim=c(0,6000))
dev.off()

######save#####c
write.csv(data.frame(var_aa,var_AA),"var_aa.csv")
write.csv(t_var_AA,"var_aa.csv")
