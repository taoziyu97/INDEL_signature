# 加载数据
library(TCGAmutations)
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
tcga_load("LUAD")
str(tcga_luad_mc3, max.level = 2)
data <- tcga_luad_mc3@data
###### Part1. 挑出mutation数据库中INDEL部分
data_INDEL <- data[Variant_Type %in% c("INS", "DEL")]
# 取第27行的数据进行test
data_sample <- data_INDEL[c(2),]
indel_types = c("1:Del:C:0", "1:Del:C:1", "1:Del:C:2", "1:Del:C:3", "1:Del:C:4", "1:Del:C:5",
                "1:Del:T:0", "1:Del:T:1", "1:Del:T:2", "1:Del:T:3", "1:Del:T:4", "1:Del:T:5",
                "1:Ins:C:0", "1:Ins:C:1", "1:Ins:C:2", "1:Ins:C:3", "1:Ins:C:4", "1:Ins:C:5",
                "1:Ins:T:0", "1:Ins:T:1", "1:Ins:T:2", "1:Ins:T:3", "1:Ins:T:4", "1:Ins:T:5", 
                # >1bp INDELS
                "2:Del:R:0", "2:Del:R:1", "2:Del:R:2", "2:Del:R:3", "2:Del:R:4", "2:Del:R:5",
                "3:Del:R:0", "3:Del:R:1", "3:Del:R:2", "3:Del:R:3", "3:Del:R:4", "3:Del:R:5",
                "4:Del:R:0", "4:Del:R:1", "4:Del:R:2", "4:Del:R:3", "4:Del:R:4", "4:Del:R:5",
                "5:Del:R:0", "5:Del:R:1", "5:Del:R:2", "5:Del:R:3", "5:Del:R:4", "5:Del:R:5",
                "2:Ins:R:0", "2:Ins:R:1", "2:Ins:R:2", "2:Ins:R:3", "2:Ins:R:4", "2:Ins:R:5", 
                "3:Ins:R:0", "3:Ins:R:1", "3:Ins:R:2", "3:Ins:R:3", "3:Ins:R:4", "3:Ins:R:5", 
                "4:Ins:R:0", "4:Ins:R:1", "4:Ins:R:2", "4:Ins:R:3", "4:Ins:R:4", "4:Ins:R:5",
                "5:Ins:R:0", "5:Ins:R:1", "5:Ins:R:2", "5:Ins:R:3", "5:Ins:R:4", "5:Ins:R:5",
                #MicroHomology INDELS
                "2:Del:M:1", "3:Del:M:1", "3:Del:M:2", "4:Del:M:1", "4:Del:M:2", "4:Del:M:3",
                "5:Del:M:1", "5:Del:M:2", "5:Del:M:3", "5:Del:M:4", "5:Del:M:5", "2:Ins:M:1", 
                "3:Ins:M:1", "3:Ins:M:2", "4:Ins:M:1", "4:Ins:M:2", "4:Ins:M:3", "5:Ins:M:1", 
                "5:Ins:M:2", "5:Ins:M:3", "5:Ins:M:4", "5:Ins:M:5", "complex", "non_matching")
ref <- data_sample$Reference_Allele
mut <- data_sample$Tumor_Seq_Allele2
indel_type_count <- 0
mutation_ID <- data.table(indel_types,indel_type_count)
###### Tips:为了减少条件语句的使用，将该样本突变处有序列的部分赋值给ID_type
if (ref != '-' & mut == '-'){
  ID_type <- ref
}else if(mut != '-' & ref == '-'){
  ID_type <- mut
}
###### 先将complex分出来
if(ref != '-' & mut != '-'){
  real_type <- 'complex'
  for (i in 1:96){
    if(real_type == indel_types[i]){
      mutation_ID$indel_type_count[i] <- mutation_ID$indel_type_count[i] + 1
      break
    }
  }
}else{
###### Part2.1 记录mutation类型是ins还是del
  if (mut == '-'){
    mut_type = 'Del'
    mut_base = substr(ref,1,1)
  }else{
    mut_type = 'Ins'
    mut_base = substr(mut,1,1)
  }
  ###### Part3. 分类并计数indel_type
  ###### Part3.1 取上下游序列
  data_loc = data_sample[, list(Chromosome=paste0("chr",Chromosome),Start=Start_Position-50,End = End_Position+250)]
  get_seq <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg19,
                              names = data_loc$Chromosome,
                              start = data_loc$Start,
                              end = data_loc$End)
  upstream = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,
                              names = data_loc$Chromosome,
                              start = data_loc$Start,
                              end = data_loc$Start+49)
  downstream = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,
                                names = data_loc$Chromosome,
                                start = data_loc$End-249,
                                end = data_loc$End)
  ###### Part3.2 计数突变处序列长度:length_INDEL
  length_INDEL <- 0
  if (nchar(ID_type) < 5){
    length_INDEL <- nchar(ID_type)
  }else if(nchar(ID_type) >= 5){
    length_INDEL <- 5
  }
  ###### Part3.3 记录和突变位点匹配的次数
  # a = "AA"
  # # b = "AAAACCTT"
  # paste0()
  # sub("(^AA+).*", "\\1", b) 
  # sub("^((AA)+).*$", "\\1", "AAAACC") 
  
  count_repeat = function(x, sequence) {
    if (startsWith(sequence, x)) {
      pattern = paste0(
        "^((",
        x,
        ")+).*$"
      )
      nchar(sub(pattern,"\\1", sequence)) / nchar(x)
    } else {
      return(0L)
    }
  }
  count_repeat <- count_repeat(ID_type,as.character(downstream))

  ###### Part3.4 记录homology次数
  count_homology_size = function(x, upstream, downstream) {
    size = nchar(x)
    if (size < 2) {
      stop("Bad input!")
    }
    x_down = substr(x, 1, size - 1)
    
    size_down = 0
    while(nchar(x_down) > 0) {
      if (startsWith(downstream, x_down)) {
        size_down = nchar(x_down)
        break()
      } else {
        x_down = substr(x_down, 1, nchar(x_down) -1)
      }
    }
    
    # Upstream
    # ACAAC|TC|AAGCGGC 
    x_up = substring(x, 2)
    size_up = 0
    while(nchar(x_up) > 0) {
      if (endsWith(upstream, x_up)) {
        size_up = nchar(x_up)
        break()
      } else {
        x_up = substring(x_up, 2)
      }
    }
    
    max(size_up, size_down)
  }
  count_homology_size <- count_homology_size(ID_type,as.character(upstream), as.character(downstream))
  ###### Part3.4 判断单位点还是多位点，区别计数
  ###### repeat_num是重复数的计数（单个位点和多个位点的homology计数方法不同）
  repeat_num <- 0
  ###### Part3.4.1 单个位点
  # 转化整个链（针对1bp indel）
  if (nchar(ID_type) == 1){
    if (mut_base == 'G'){
      mut_base = 'C'
      strand <- '-1'
    }else if(mut_base == 'A'){
      mut_base = 'T'
      strand <- '-1'
    }
    if(count_repeat > 5){
      count_repeat <- 5
    }
    ###### Part3.4.2 多个位点
  }else if(nchar(ID_type) > 1){
    mut_base <- 'R'
    if(count_repeat > 5){
      count_repeat <- 5
      }
    }else if(count_homology_size > 0){
      mut_base <- 'M'
      if(count_homology_size > 5){
        count_homology_size <- 5
      }
    }
  ###### Part4 整合indel_type为标准格式，记录为real_types
  if(mut_base == 'R' & mut_type == 'Del'){
    real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_repeat+1,sep="")
  }else if(mut_base == 'R' & mut_type == 'Ins'){
    real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_repeat,sep="")
  }else{
    real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",count_homology_size,sep="")
  }
  ###### Part5 对indel_type进行计数
  indel_type_count <- 0
  mutation_ID <- data.table(indel_types,indel_type_count)
  for (i in 1:96){
    if(real_type == indel_types[i]){
      mutation_ID$indel_type_count[i] <- mutation_ID$indel_type_count[i] + 1
    }
  }
}

mutation_ID[indel_type_count > 0]


