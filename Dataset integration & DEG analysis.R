my_cel_files1 = ReadAffy()



order_GSE54958 <- data.frame(
  R_order = seq_along(sampleNames(my_cel_files1)),
  CEL_file_read_by_R = sampleNames(my_cel_files1),
  GSM = sub("^(GSM[0-9]+).*", "\\1", basename(sampleNames(my_cel_files1))),
  ManualGroup = c(rep("FA", 7), rep("FVPTC", 6)),
  stringsAsFactors = FALSE
)

dir.create("ml-output", showWarnings = FALSE, recursive = TRUE)

write.csv(
  order_GSE54958,
  file.path("ml-output", "order_GSE54958_manual_check.csv"),
  row.names = FALSE
)
print(order_GSE54958)









norm_data_with_rma1 = affy::rma(my_cel_files1)
exprs_data1 = exprs(norm_data_with_rma1)
final_data1 = exprs_data1

annotation_file1 = read.csv("GPL6244-17930-modified.csv")
annot_df1 = c()
annot_df1 = cbind(annot_df1, annotation_file1$Gene.Symbol)
annot_df1 = cbind(annot_df1, annotation_file1$GB_LIST)
rownames(annot_df1) = annotation_file1$ID
colnames(annot_df1) = c("Symbol", "Title")
annot_df1 = data.frame(annot_df1)
prime.view1 = function(x){
  return(annot_df1[x, 1])
}
symbol1 = prime.view1(rownames(final_data1)[1:32321])
final_data1 = as.data.frame(final_data1)
final_data1$Gene = symbol1
final_data1 = final_data1[, c("Gene", setdiff(names(final_data1), "Gene"))]
final_data1 = final_data1[final_data1$Gene != "#VALUE!",]
final_data1 = aggregate(.~ Gene, data = final_data1, FUN = mean)
rownames(final_data1) = final_data1$Gene
final_data1$Gene = NULL

shared_genes = intersect(rownames(final_data), rownames(final_data1))
final_data = final_data[rownames(final_data) %in% shared_genes, ]
final_data1 = final_data1[rownames(final_data1) %in% shared_genes, ]
final_data = cbind(final_data, final_data1)
write.csv(final_data, "final_data.csv")

samples = colnames(final_data)
batch = c(rep("GSE29315", 30), rep("GSE54958", 13))
condition = c(rep("FA", 17), rep("FVPTC", 13), rep("FA", 7), rep("FVPTC", 6))
metadata = data.frame(sampleID = samples, batch = batch, condition = condition)
mod = model.matrix(~condition, data=metadata)

batch_group = c(rep("GSE29315", 30), rep("GSE54958", 13))
final_data_combat = ComBat(final_data, batch_group, mod = mod)

Groups = c(rep("FA", 17), rep("FVPTC", 13), rep("FA", 7), rep("FVPTC", 6))
Design = model.matrix(~0+factor(Groups))
colnames(Design) = c("FA", "FVPTC")
fit = lmFit(final_data_combat, Design)
cont = makeContrasts(contrasts = "FVPTC - FA", levels = Design)
fit2 = contrasts.fit(fit = fit, contrasts = cont)
fit3 = eBayes(fit2, 0.01)
Result = topTable(fit3, adjust="fdr", sort.by="B", number= Inf)
Result$Gene = rownames(Result)

FVPTC.up = subset(Result, logFC > 1 & adj.P.Val < 0.05)
FVPTC.up.genes = unique(as.character(strsplit2(FVPTC.up$Gene, "///")))
write.table(FVPTC.up.genes, file = "FVPTC-FA.up.txt", quote = F, row.names = F, col.names = F)
FVPTC.down = subset(Result, logFC < -1 & adj.P.Val < 0.05)
FVPTC.down.genes = unique(as.character(strsplit2(FVPTC.down$Gene, "///")))
write.table(FVPTC.down.genes, file = "FVPTC-FA.down.txt", quote = F, row.names = F, col.names = F)
DEGs = c(FVPTC.up.genes, FVPTC.down.genes)
write.table(DEGs, file = "DEGS.txt", quote = F, row.names = F, col.names = F)

final_data <- final_data_combat
saveRDS(final_data, "final_data_clean.rds")
write.csv(final_data, "final_data_clean.csv")