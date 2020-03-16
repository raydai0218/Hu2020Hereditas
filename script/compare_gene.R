#!/usr/bin/env Rscript
# 
# Copyright 2016-2018 Yong-Xin Liu <metagenome@126.com>

# If used this script, please cited:
#   Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
# Root microbiota shift in rice correlates with resident time in the field and developmental stage. Sci China Life Sci 61, 
# https://doi.org/10.1007/s11427-018-9284-4

# 手动运行脚本请，需要设置工作目录，使用 Ctrl+Shift+H 或 Session - Set Work Directory - Choose Directory / To Source File Location 设置工作目录



# 1.1 程序功能描述和主要步骤

# 程序功能：高通量测序reads counts值的组间比较并筛选
# Functions: Calculate pvalue and FDR for each OTUs by edgeR or wilcon
# Main steps: 
# - Reads data matrix and design
# - Calculate pvalue and FDR
# - Save result table in *_all/sig.txt

# 清空工作环境 Clean enviroment object
rm(list=ls()) 


# 2.1 安装CRAN来源常用包
site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# 依赖包列表：参数解析、数据变换、绘图和开发包安装、安装依赖、ggplot主题
package_list = c("limma","ggplot2","pheatmap","dplyr","devtools")
# 判断R包加载是否成功来决定是否安装后再加载
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.2 安装bioconductor常用包
package_list = c("edgeR")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    source("https://bioconductor.org/biocLite.R")
    biocLite(p)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

# 2.3 安装Github常用包
# 参数解析、数据变换、绘图和开发包安装
package_list = c("kassambara/ggpubr")
for(p in package_list){
  q=unlist(strsplit(p,split = "/"))[2]
  if(!suppressWarnings(suppressMessages(require(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install_github(p)
    suppressWarnings(suppressMessages(library(q, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# 3. 读取输入文件
write.table("ID\ttype", file=paste("diff.list", sep=""), append = F, sep="\t", quote=F, row.names=F, col.names=F)

# 读取实验设计
design = read.table("../data/design.txt", header=T, row.names= 1, sep="\t") 
# 统一改实验列为group
design$group = design$groupID

# 按实验组筛选 Select by manual set group,A50的
if (TRUE){
  design = design
  # 调置组排序 Set group order
  #design$group  = factor(design$group, levels=c(""))
}

# 读取OTU表
otutab = read.table("../data/otutab_normed.txt", header=T, row.names=1, quote = "", sep="\t", comment.char="") 
# Show features (include OTU/ASV/Taxonomy) numbers
print("Total features number")
print(dim(otutab)[1])

# 实验设计与输入文件交叉筛选
idx = rownames(design) %in% colnames(otutab)
design = design[idx, , drop = F]
otutab = otutab[,rownames(design)]

# 按丰度值按组中位数筛选OTU
# 标准化为百分比例，并转置
if (TRUE){
  norm = t(otutab)/colSums(otutab,na=T)*100
}else{
  # 非标准化时为，默认抽样10000，除以100标准为百分比
  norm=t(otutab)*1
}
# 检查样本标准化后是否为100
# rowSums(norm)

# 筛选组，按组求中位数
# grp = design[, "GroupID", drop=F] # 需要按第二条件筛选时使用
grp = design[, "groupID", drop=F]
# 按行名合并
mat_t2 = merge(grp, norm, by="row.names")
mat_t2 = mat_t2[,-1]
# 按组求中位数，中位数筛选更有效去除异常值
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=median) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno
# 按丰度按组中位数筛选,万分之一
filtered = mat_mean_final[apply(mat_mean_final,1,max) >= 0.01, ] # select OTU at least one sample > 0.01%
otutab = otutab[rownames(filtered),]

# 按均值输出和保存相对丰度，汇总才为真实丰度
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]
geno = mat_mean$group
colnames(mat_mean_final) = geno

# 生成compare的database用于注释
mat_mean_high = mat_mean_final[rownames(filtered),]
write.table(paste("OTUID\t",sep=""), file=paste("database.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
suppressWarnings(write.table(round(mat_mean_high,5), file=paste("database.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T))

print("Selected high abundance OTUs number")
print(dim(mat_mean_high)[1])#461

print(colSums(mat_mean_high))

print(paste("你正在使用秩和检验！Now, you are using wilcoxon test!", sep=" "))

compare_DA = function(compare){
  # 筛选比较组wilcox
  group_list = as.vector(as.matrix(compare))
  SampAvsB=paste(group_list[1] ,"-", group_list[2], sep="")
  idx = design$group %in% group_list
  sub_design=design[idx, , drop = F]
  #	sub_dat=as.matrix(otutab[,rownames(sub_design)])
  #
  #	# wilcoxon秩合检验，需要先标准化
  #	# normlization to percentage
  #	if (TRUE){
  #		sub_norm = t(t(sub_dat)/colSums(sub_dat,na=T))*100
  #	}else{
  #		sub_norm = as.matrix(otutab) * 1 # 数据类型一致，计算后矩阵
  #	}
  norm = as.data.frame(t(norm))
  sub_norm = as.matrix(norm[rownames(filtered),])
  # 建立两组的矩阵
  idx = sub_design$group %in% group_list[1]
  GroupA = sub_norm[,rownames(sub_design[idx,,drop=F])]
  idx = sub_design$group %in% group_list[2]
  GroupB = sub_norm[,rownames(sub_design[idx,,drop=F])]
  
  nrDAO = data.frame(list=rownames(sub_norm), row.names =rownames(sub_norm) )
  # 对每行OUT/基因进行秩合检验
  # dim(nrDAO)[1]
  for ( i in 1:dim(nrDAO)[1]){
    FC = (mean(GroupA[i,])+0.0001)/(mean(GroupB[i,])+0.0001)
    nrDAO[i,2]=log2(FC)
    nrDAO[i,3]=log2(max(c(GroupA[i,],GroupB[i,]))*10000)
    nrDAO[i,4]= wilcox.test(as.numeric(GroupA[i,]),as.numeric(GroupB[i,]))$p.value
  }	
  nrDAO=nrDAO[,-1]
  colnames(nrDAO)=c("logFC", "logCPM", "PValue")
  nrDAO$FDR = p.adjust(nrDAO$PValue, method="fdr", dim(nrDAO)[1])    #p.adjust就是计算FDR的包，这个可要记得了
  
  
  # 整理数据格式
  nrDAO$logFC=round(nrDAO$logFC,3)
  nrDAO$logCPM=round(nrDAO$logCPM,3)
  nrDAO$level = ifelse(nrDAO$logFC>log2(1.2) & nrDAO$PValue<0.05 & nrDAO$FDR<0.2, "Enriched",ifelse(nrDAO$logFC<log2(1.2)*(-1) & nrDAO$PValue<0.05 & nrDAO$FDR<0.2, "Depleted","NotSig"))
  nrDAO$level=factor(nrDAO$level,levels = c("Enriched","Depleted","NotSig"))
  
  # Add MeanA and MeanB in percentage
  # calculate groupA mean
  A_list = subset(sub_design, group %in% group_list[1])
  A_norm = sub_norm[, rownames(A_list)]
  A_mean = as.data.frame(rowMeans(A_norm))
  colnames(A_mean)=c("MeanA")
  # calculate groupB mean
  B_list = subset(sub_design, group %in% group_list[2])
  B_norm = sub_norm[, rownames(B_list)]
  B_mean = as.data.frame(rowMeans(B_norm))
  colnames(B_mean)=c("MeanB")
  # merge and reorder
  Mean = round(cbind(A_mean, B_mean, A_norm, B_norm),3)
  Mean = Mean[rownames(nrDAO),]   
  
  # 存在物种注释，添加至Mean前
  if (file.exists("../data/taxonomy_8.txt")){
    tax = read.table("../data/taxonomy_8.txt", header=T, row.names= 1, sep="\t", comment.char = "") 
    tax = tax[rownames(nrDAO),]
    Mean=cbind(tax, Mean)
  }
  
  output=cbind(nrDAO,Mean)
  
  # write all OTU for volcano plot and manhattan plot
  write.table(paste(group_list[1],"_", group_list[2], "\t",sep=""), file=paste(SampAvsB, "_all.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
  suppressWarnings(write.table(output,file=paste(SampAvsB, "_all.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T))
  
  # 计算上、下调OTUs数量，写入统计文件
  NoE= dim(output[output$level=="Enriched",])[1]
  NoD= dim(output[output$level=="Depleted",])[1]
  NoN= dim(output[output$level=="NotSig",])[1]
  suppressWarnings(write.table(paste( SampAvsB, NoE, NoD, NoN, sep="\t"), file=paste("summary.txt",sep=""), append = T, quote = F, sep = '\t', row.names = F, col.names = F))
  
  output=output[output$level!="NotSig",]
  # 保存筛选结果于sig.txt结尾文件中，无差异不输出
  if (dim(output)[1]>1){
    write.table(paste(group_list[1],"_", group_list[2], "\t",sep=""), file=paste(SampAvsB, "_sig.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)
    suppressWarnings(write.table(output, file=paste(SampAvsB, "_sig.txt",sep=""), append = T, quote = F, sep = '\t', row.names = T))
    
    # 确保有差异才写入，否则会出不完整行
    # 保存差异列表用于维恩图展示 save each group DA OTU list for venndiagram
    # 使用edgeR中的比较-相连，无法作为变量赋值，改为_
    write.table(cbind(rownames(output),paste(group_list[1],"_", group_list[2], output$level, sep="")), file=paste("diff.list", sep=""), append = TRUE, sep="\t", quote=F, row.names=F, col.names=F)
  }
}


# 记录各组间上、下调数量
write.table("GroupAvsB\tEnriched\tDepleted\tNotSig\n", file=paste("summary.txt",sep=""), append = F, quote = F, eol = "", row.names = F, col.names = F)

# 如果没有比较文件，则自动全循环
if (!file.exists("compare_gene.txt")) {
  compare_data = as.vector(unique(design$group))
  len_compare_data = length(compare_data)
  for(i in 1:(len_compare_data-1)) {
    for(j in (i+1):len_compare_data) {
      tmp_compare = as.data.frame(cbind(sampA=compare_data[i],sampB=compare_data[j]))
      compare_DA(tmp_compare)
    }
  }
}else {
  compare_data = read.table("compare_gene.txt", sep="\t", check.names=F, quote='', comment.char="")
  colnames(compare_data) = c("sampA", "sampB")
  for(i in 1:dim(compare_data)[1]){
    compare_DA(compare_data[i,])
    print(paste("Compared", compare_data[i,1], "vs", compare_data[i,2],"finished!", sep=" "))
  }
}

system("sed -i 's/Enriched/_E/;s/Depleted/_D/;' diff.list")

