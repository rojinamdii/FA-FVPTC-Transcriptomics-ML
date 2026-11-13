library(limma)
library(reshape2)
library(plyr)
library(ggplot2)
library(pheatmap)
library(gplots)
library(affy)
library(sva)

my_cel_files = ReadAffy()









order_GSE29315 <- data.frame(
  R_order = seq_along(sampleNames(my_cel_files)),
  CEL_file_read_by_R = sampleNames(my_cel_files),
  GSM = sub("^(GSM[0-9]+).*", "\\1", basename(sampleNames(my_cel_files))),
  ManualGroup = c(rep("FA", 17), rep("FVPTC", 13)),
  stringsAsFactors = FALSE
)

dir.create("ml-output", showWarnings = FALSE, recursive = TRUE)

write.csv(
  order_GSE29315,
  file.path("ml-output", "order_GSE29315_manual_check.csv"),
  row.names = FALSE
)
print(order_GSE29315)






















norm_data_with_rma = affy::rma(my_cel_files)
exprs_data = exprs(norm_data_with_rma)
final_data = exprs_data


#Annotation
annotation_file = read.csv("GPL8300-10787.csv")
annot_df = c()
annot_df = cbind(annot_df, annotation_file$Gene.Symbol)
annot_df = cbind(annot_df, annotation_file$Gene.Title)
rownames(annot_df) = annotation_file$ID
colnames(annot_df) = c("Symbol", "Title")
annot_df = data.frame(annot_df)
prime.view = function(x){
  return(annot_df[x, 1])
}
symbol = prime.view(rownames(final_data)[1:12625])
final_data = as.data.frame(final_data)
final_data$Gene = symbol
final_data = final_data[, c("Gene", setdiff(names(final_data), "Gene"))]
final_data = final_data[final_data$Gene != "",]
final_data = aggregate(.~ Gene, data = final_data, FUN = mean)
rownames(final_data) = final_data$Gene
final_data$Gene = NULL