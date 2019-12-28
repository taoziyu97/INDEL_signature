library(TCGAmutations)
tcga_load("LUAD")
library(BSgenome.Hsapiens.UCSC.hg19)
library(data.table)
str(tcga_luad_mc3, max.level = 2)
#不需要挑出突变位点为TC之外的序列，因为INDEl处不止有1个变化，无法进行简化并做reverse这一步骤（参考SNV的96种）
# 把MAF文件格式进行标准化使得原作者python脚本能够识别
Entrez_Gene_Id <- c("-")
Center <- c("-")
NCBI_Build <- c("GRCh37")
Strand <- c("-")
Tumor_Seq_Allele1 <- c("-")
dbSNP_RS <- c("-")
dbSNP_Val_Status <- c("-")
data<-cbind(tcga_luad_mc3@data[,1],Entrez_Gene_Id,
            Center,NCBI_Build,tcga_luad_mc3@data[,2:4],
            Strand,tcga_luad_mc3@data[,7:8],
            tcga_luad_mc3@data[,5],Tumor_Seq_Allele1,
            tcga_luad_mc3@data[,6],dbSNP_RS,dbSNP_Val_Status,
            tcga_luad_mc3@data[,9:ncol(tcga_luad_mc3@data)])
###### Part1. 挑出mutation数据库中INDEL部分
data_INDEL <- data[Variant_Type %in% c("INS", "DEL")]
# 取第5行的数据进行test
data_sample <- data_INDEL[c(32),]
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
sample <- data_sample$Tumor_Sample_Barcode
chrom <- data_sample$Chromosome
start <- data_sample$Start_Position
ref <- data_sample$Reference_Allele
mut <- data_sample$Tumor_Seq_Allele2
variant_type <- data_sample$Variant_Type
###### Tips:为了减少条件语句的使用，将该样本突变处有序列的部分赋值给ID_type
ID_type <- ''
if (ref != '-' & mut == '-'){
  ID_type <- ref
}else if(mut != '-' & ref == '-'){
  ID_type <- mut
}
###### 先将complex分出来
if(ref != '-' & mut != '-' & (variant_type == 'DEL' | variant_type == 'INS')){
  real_type <- 'complex'
  indel_type_count <- 0
  mutation_ID <- data.table(indel_types,indel_type_count)
  for (i in 1:96){
    if(real_type == indel_types[i]){
      mutation_ID$indel_type_count[i] <- mutation_ID$indel_type_count[i] + 1
    }
  }
}else{
  ###### Part2.1 记录mutation类型是ins还是del
  ###### Tips:加上判断mut或ref等于‘-’与否，因为单个位点的indel字符数也是1，需要加上这个额外条件
  if (length(ref) - length(mut) == length(ref)-1 & mut == '-'){
    mut_type = 'Del'
    mut_base = substr(ref,1,1)
  }else if(length(mut) - length(ref) == length(mut)-1 & ref == '-'){
    mut_type = 'Ins'
    mut_base = substr(mut,1,1)
  }
  ######Part3. 分类并计数indel_type
  ###### Part3.1 取上下游序列
  str(tcga_luad_mc3, max.level = 2)
  data = tcga_luad_mc3@data
  data = data[Variant_Type %in% c("INS", "DEL")]
  data_loc = data_sample[, list(Chromosome=paste0("chr",Chromosome),Start=Start_Position-25,End = End_Position+25)]
  get_seq <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg19,
                              names = data_loc$Chromosome,
                              start = data_loc$Start,
                              end = data_loc$End)
  upstream = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,
                              names = data_loc$Chromosome,
                              start = data_loc$Start,
                              end = data_loc$Start+24)
  downstream = BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg19,
                                names = data_loc$Chromosome,
                                start = data_loc$End-24,
                                end = data_loc$End)
  ###### Part3.2 计数突变处序列长度:length_INDEL
  length_INDEL <- 0
  if (nchar(ID_type) < 5){
    length_INDEL <- nchar(ID_type)
  }else if(nchar(ID_type) >= 5){
    length_INDEL <- 5
  }
  ###### Part3.3 记录和突变位点匹配的序列:match_sequence
  ###### 先计数多少个序列匹配：number_down，再按照具体分类区分：real_number_down
  number_down <- 0
  match_sequence <- ''
  for (i in seq(1,50,1)){
    repeat_ID_type <- paste0(replicate(round(50/(nchar(ID_type)))+1, ID_type),collapse = "")
    if (substr(repeat_ID_type,i,i) == as.character(substr(downstream,i,i))){
      number_down <- number_down+1
      match_sequence <- paste(match_sequence,substr(repeat_ID_type,i,i),sep = "")
      next
    }else{
      break
    }
  }
  for (i in seq(1,50,1)){
    if (substr(repeat_ID_type,round(50/(nchar(ID_type)))+2-i,round(50/(nchar(ID_type)))+2-i) == as.character(substr(upstream,round(50/(nchar(ID_type)))+2-i,round(50/(nchar(ID_type)))+2-i))){
      up_number <- up_number+1
      match_sequence <- paste(match_sequence,substr(repeat_ID_type,round(50/(nchar(ID_type)))+2-i,round(50/(nchar(ID_type)))+2-i),sep = "")
      next
    }else{
      break
    }
  }
  number <- 0
  up_number <- 0
  if (up_number > number_down){
    number <- up_number
  }else if(up_number <= number_down){
    number <- number_down
  }
  ###### 计数以突变位点为一个单位的重复次数：real_num
  ###### 对于单个位点，虽然进行这步，但是该步骤得到的real_num不计入mutation_ID表格中，不用在意
  real_num <- 0
  count_num <- nchar(ID_type)
  if (number_down > 0){
    for (i in seq(1,nchar(match_sequence),count_num)){
      if(ID_type == substr(match_sequence,i,i+count_num-1)){
        real_num <- real_num+1
      }
    }
  }
  ###### Part3.4 判断单位点还是多位点，区别计数
  ###### repeat_num是重复数的计数（单个位点和多个位点的homology计数方法不同）
  repeat_num <- 0
  ###### Part3.4.1 单个位点
  if (nchar(ID_type) == 1){
    if (mut_base == 'G'){
      mut_base = 'C'
      strand <- '-1'
    }else if(mut_base == 'A'){
      mut_base = 'T'
      strand <- '-1'
    }
    if (number_down < 5){
      repeat_num <- number_down
    }else if(number_down >= 5){
      repeat_num <- 5
    }
    ###### Part3.4.2 多个位点
  }else if(nchar(ID_type) > 1){
    if(nchar(match_sequence) >= nchar(ID_type)){
      mut_base <- 'R'
      if (number_down < 5){
        repeat_num <- number_down
      }else if(number_down >= 5){
        repeat_num <- 5
      }
    }else if(nchar(match_sequence) < nchar(ID_type)){
      mut_base <- 'M'
      if (number_down < 5){
        repeat_num <- number_down
      }else if(number_down >= 5){
        repeat_num <- 5
      }
    }
  }
  
  ###### Part4 整合indel_type为标准格式，记录为real_types
  if(mut_base == 'R' & mut_type == 'Del'){
    real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",real_num+1,sep="")
  }else if(mut_base == 'R' & mut_type == 'Ins'){
    real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",real_num,sep="")
  }else{
    real_type <- paste(length_INDEL,":",mut_type,":",mut_base,":",number_down,sep="")
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
