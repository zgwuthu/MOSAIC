library(ShortRead)
library(Biostrings)
rev_comp<-function(dna){c
  return(
    paste(rev(sapply(unlist(strsplit(dna,split="")), function(x) {switch(x, A = "T", T = "A", G = "C", C = "G") })),collapse = ""))
}#get the reverse complementary sequence of DNA
translate_dna <- function(dna_sequence) {c
  # 定义兼并碱基的所有可能对应碱基
  dna_sequence=toupper(dna_sequence)
  dna_sequence=unlist(strsplit(dna_sequence,""))
  degenerate_bases <- list(
    "R" = c("A", "G"),
    "Y" = c("C", "T"),
    "S" = c("G", "C"),
    "W" = c("A", "T"),
    "K" = c("G", "T"),
    "M" = c("A", "C"),
    "B" = c("C", "G", "T"),
    "D" = c("A", "G", "T"),
    "H" = c("A", "C", "T"),
    "V" = c("A", "C", "G"),
    "N" = c("A", "C", "G", "T")
  )
  
  # 检查输入是否为DNA序列
  if (!all(dna_sequence %in% c("A", "C", "G", "T", names(degenerate_bases)))) {
    stop("输入包含无效的碱基")
  }
  
  # 展开兼并碱基的所有可能组合
  all_combinations <- list()
  
  for (i in seq_along(dna_sequence)) {
    base <- dna_sequence[i]
    if (base %in% names(degenerate_bases)) {
      # 如果是兼并碱基，添加所有可能对应碱基
      if (length(all_combinations) == 0) {
        all_combinations <- lapply(degenerate_bases[[base]], function(x) x)
      } else {
        new_combinations <- list()
        for (comb in all_combinations) {
          for (alt in degenerate_bases[[base]]) {
            new_combinations[[length(new_combinations) + 1]] <- paste0(comb, alt)
          }
        }
        all_combinations <- new_combinations
      }
    } else {
      # 如果是普通碱基，保持不变
      if (length(all_combinations) == 0) {
        all_combinations <- list(base)
      } else {
        all_combinations <- lapply(all_combinations, function(x) paste0(x, base))
      }
    }
  }
  
  # 定义密码子表
  codon_table <- c(
    AAA = "K", AAC = "N", AAG = "K", AAT = "N",
    CAA = "Q", CAC = "H", CAG = "Q", CAT = "H",
    GAA = "E", GAC = "D", GAG = "E", GAT = "D",
    TAA = "*", TAC = "Y", TAG = "*", TAT = "Y",
    ACA = "T", ACC = "T", ACG = "T", ACT = "T",
    CCA = "P", CCC = "P", CCG = "P", CCT = "P",
    GCA = "A", GCC = "A", GCG = "A", GCT = "A",
    TCA = "S", TCC = "S", TCG = "S", TCT = "S",
    AGA = "R", AGC = "S", AGG = "R", AGT = "S",
    CGA = "R", CGC = "R", CGG = "R", CGT = "R",
    GGA = "G", GGC = "G", GGG = "G", GGT = "G",
    TGA = "*", TGC = "C", TGG = "W", TGT = "C",
    ATA = "I", ATC = "I", ATG = "M", ATT = "I",
    CTA = "L", CTC = "L", CTG = "L", CTT = "L",
    GTA = "V", GTC = "V", GTG = "V", GTT = "V",
    TTA = "L", TTC = "F", TTG = "L", TTT = "F"
  )
  
  # 翻译所有可能的DNA序列
  all_proteins <- list()
  
  for (dna in all_combinations) {
    protein <- ""
    
    # 确保序列长度是3的倍数
    if (nchar(dna) %% 3 != 0) {
      warning("DNA序列长度不是3的倍数，将截断到最后一个完整的密码子")
      dna <- substr(dna, 1, floor(nchar(dna)/3)*3)
    }
    
    # 翻译密码子
    for (i in seq(1, nchar(dna), by=3)) {
      codon <- substr(dna, i, i+2)
      if (codon %in% names(codon_table)) {
        protein <- paste0(protein, codon_table[[codon]])
      } else {
        warning(paste("未知的密码子:", codon))
        protein <- paste0(protein, "X")  # 用X表示未知氨基酸
      }
    }
    
    all_proteins[[length(all_proteins) + 1]] <- protein
  }
  
  # 返回所有可能的多肽序列
  return(all_proteins)
}#get all posible residue combinations from one degenerated DNA sequence

#####reads#####
f_reads <- readFastq("FP-LHG19425-w1_L0_1.fq")
r_reads <- readFastq("FP-LHG19425-w1_L0_2.fq")
f_reads=sread(f_reads)
r_reads=sread(r_reads)
#####recognized sites#####
site1="GGCGTGCAGTGCTTC"#-9 -1 position of variants
site2="CGCTACCCCGACCACATG"#-3 -1  +1 +3 
site3="ATCCTGGGGCACAAGCTGGAGTACAAC"#+1+6 +13+15
site4="GCCGACAAGCAGAAGAACGGCATCAAG"#-3-1 +1+3 +13+15
site5="CCGACAACCACTACCTGAGC"#+1 +3
site3=reverseComplement(DNAString(site3)) #-6-1  -15-13
site4=reverseComplement(DNAString(site4)) #+1+3 -3-1 -15-13
site5=reverseComplement(DNAString(site5)) #-3-1
l=length(f_reads)
var_aa=c() #variant DNA seq
var_AA=c() #variant protein seq
for (i in 1:l) {
  f=as.character( f_reads[i] )
  r=as.character( r_reads[i] )
  pos1=matchPattern(site1,f)
  pos2=matchPattern(site2,f)
  pos3=matchPattern(site3,r)
  pos4=matchPattern(site4,r)
  pos5=matchPattern(site5,r)
  check=!is.na(start(pos1)[1]) && !is.na(start(pos2)[1]) && !is.na(start(pos3)[1]) && !is.na(start(pos4)[1]) && !is.na(start(pos5)[1])
  if(check==F){f1=r;r1=f
  pos1=matchPattern(site1,f1)
  pos2=matchPattern(site2,f1)
  pos3=matchPattern(site3,r1)
  pos4=matchPattern(site4,r1)
  pos5=matchPattern(site5,r1)
  check=!is.na(start(pos1)[1]) && !is.na(start(pos2)[1]) && !is.na(start(pos3)[1]) && !is.na(start(pos4)[1]) && !is.na(start(pos5)[1])
  if(check==T){f=f1;r=r1
  }
  }
  if(check){
    
    aa=paste( c(  substr(f,start(pos1)[1]-9,start(pos1)[1]-1),
                  substr(f,start(pos2)[1]-3,start(pos2)[1]-1),
                  substr(f,end(pos2)[1]+1,end(pos2)[1]+3),
                  rev_comp(paste(c(substr(r,start(pos5)[1]-3,start(pos5)[1]-1),
                                   substr(r,start(pos4)[1]-15,start(pos4)[1]-13),
                                   substr(r,start(pos4)[1]-3,start(pos4)[1]-1),
                                   substr(r,end(pos4)[1]+1,end(pos4)[1]+3),
                                   substr(r,start(pos3)[1]-15,start(pos3)[1]-13),
                                   substr(r,start(pos3)[1]-6,start(pos3)[1]-1)),collapse=""))
                  
    ),collapse = "")
    if(grepl("NULL", aa)==T){next}
    
    var_aa=append(var_aa,aa) 
    
    var_AA=append(var_AA,as.character( translate(DNAString(aa)) ))
  }
  if(i%%10000==0){cat("Complete:",i/l*100,"%","\n")}
}
#####recognized designed sequences#####
total_AA=unique(translate_dna("TTNDSCYRKKCAARGTWCAWCAANAYGGYGAYCWMC"))#all designed sequences

count=0
for(i in total_AA){
  if(i%in%u_var_AA){count=count+1}
}

print(c("Number of all designed seuqences:",length(total_AA)))
print(c("Number of all recognized seuqences:",count))
print(c("Ratio:",count/length(total_AA)))

#save counts of individual designed sequences
n_total_AA=rep(0,length(total_AA))
for (i in seq_along(total_AA)) {
  n_total_AA[i]=as.numeric( t_var_AA[as.character(total_AA[i])] )
}
write.csv(data.frame(total_AA,n_total_AA),"designed_seq_count.csv")