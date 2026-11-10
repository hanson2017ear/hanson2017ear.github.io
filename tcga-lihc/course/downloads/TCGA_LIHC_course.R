# TCGA-LIHC 转录组分析：按课程教程逐段运行
# 姓名：陳萬傑；学号：23360036
# 顺序与“教程-敲除细胞系/教程代码.R”一致，改用 TCGA-LIHC 的配对样本。
# 完整解压后在 RStudio 打开本文件，点击 Source 可运行全部步骤。
# 逐段运行时，先用 Session > Set Working Directory > To Source File Location 设置工作目录。
# 比较始终为 GroupX / GroupY；log2FoldChange > 0 表示 GroupX 较高。

#### 第1步：设置路径并加载 R 包 ####
if (.Platform$OS.type == "windows")
  invisible(try(suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8")), silent = TRUE))
options(stringsAsFactors = FALSE, timeout = 600, warn = 1)
invisible(try(Sys.setlocale("LC_COLLATE", "C"), silent = TRUE))

project_dir <- local({
  paths <- vapply(sys.frames(), function(x) {
    if (is.character(x$ofile) && length(x$ofile) == 1) x$ofile else ""
  }, character(1))
  paths <- paths[nzchar(paths) & file.exists(paths)]
  cli <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(paths)) dirname(normalizePath(tail(paths, 1), winslash = "/")) else
    if (length(cli)) dirname(normalizePath(sub("^--file=", "", cli[1]), winslash = "/")) else getwd()
})
input_dir <- file.path(project_dir, "data")
out_dir <- file.path(project_dir, "分析结果")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "表格"), showWarnings = FALSE)
dir.create(file.path(out_dir, "中间对象"), showWarnings = FALSE)
lib_dir <- file.path(out_dir, "R_library")
dir.create(lib_dir, showWarnings = FALSE)
r_minor <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
user_lib <- file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", r_minor)
.libPaths(unique(c(lib_dir, user_lib[dir.exists(user_lib)], .libPaths())))

# 已安装的包直接使用；缺包时需联网安装。
cran_packages <- c("ggplot2", "dplyr", "pheatmap", "ragg", "jsonlite")
missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, lib = lib_dir, repos = "https://cloud.r-project.org")
bio_packages <- c("DESeq2", "clusterProfiler")
missing_bio <- bio_packages[!vapply(bio_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bio)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = lib_dir, repos = "https://cloud.r-project.org")
  BiocManager::install(missing_bio, lib = lib_dir, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({
  library(DESeq2); library(clusterProfiler); library(ggplot2)
  library(dplyr); library(pheatmap)
})
set.seed(23360036)

# 仅用于同时保存 PNG 和 PDF；分析步骤均在下方依次执行。
save_plot <- function(p, name, width = 10, height = 7) {
  if (interactive()) print(p)
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width, height = height, dpi = 320, device = ragg::agg_png, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width, height = height, device = cairo_pdf, bg = "white")
}

#### 第2步：读入 count 矩阵和分组表 ####
expr_mat <- read.csv(file.path(input_dir, "TCGA_LIHC的8个样本count.csv"), row.names = 1, check.names = FALSE)
colData <- read.csv(file.path(input_dir, "TCGA_LIHC的8个样本分组.csv"), row.names = 1, check.names = FALSE)
stopifnot(!anyDuplicated(rownames(expr_mat)), !anyDuplicated(colnames(expr_mat)),
          !anyDuplicated(rownames(colData)), setequal(colnames(expr_mat), rownames(colData)))
colData <- colData[colnames(expr_mat), , drop = FALSE]
stopifnot(identical(colnames(expr_mat), rownames(colData)))
count_mat <- as.matrix(expr_mat)
stopifnot(is.numeric(count_mat), all(is.finite(count_mat)), all(count_mat >= 0))
dim(expr_mat)
table(colData$Group)

#### 第3步：检查非整数值、低表达基因及患者配对 ####
# 140 个基因含超出 1e-6 容差的非整数值；保存清单后从主分析排除。
integer_error <- abs(count_mat - round(count_mat))
fractional <- rowSums(integer_error > 1e-6) > 0
write.csv(expr_mat[fractional, , drop = FALSE],
          file.path(out_dir, "表格", "排除的非整数基因.csv"), na = "NA", fileEncoding = "UTF-8")
keep_low <- rowSums(round(count_mat) >= 10) >= 4
keep <- !fractional & keep_low
expr_mat <- round(count_mat[keep, , drop = FALSE])
storage.mode(expr_mat) <- "integer"
filter_status <- data.frame(gene = rownames(count_mat),
  status = ifelse(fractional, "excluded_fractional", ifelse(keep_low, "modeled", "low_expression")))
write.csv(filter_status, file.path(out_dir, "表格", "基因过滤状态.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 条码前三段只用于患者配对，组织类型留到最后核对。
patient_id <- vapply(strsplit(rownames(colData), "-", fixed = TRUE),
                     function(x) paste(x[1:3], collapse = "-"), character(1))
colData$patient <- factor(patient_id, levels = sort(unique(patient_id)))
colData$Group <- factor(colData$Group, levels = c("GroupY", "GroupX"))
stopifnot(nrow(colData) == 8, nlevels(colData$patient) == 4, all(table(colData$patient, colData$Group) == 1))
colData$Patient <- paste0("P", as.integer(colData$patient))
colData$sample_label <- paste0(colData$Patient, "_", ifelse(colData$Group == "GroupX", "X", "Y"))
table(colData$patient, colData$Group)

#### 第4步：构建 DESeq2 对象并进行差异检验 ####
# 教程的 KO/WT 使用 ~ Group；本题是四位患者的配对样本，因此增加 patient。
dds <- DESeqDataSetFromMatrix(countData = expr_mat, colData = colData, design = ~ patient + Group)
dds <- DESeq(dds, parallel = FALSE)
res <- results(dds, contrast = c("Group", "GroupX", "GroupY"),
               alpha = 0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
DEG <- data.frame(gene = rownames(res), as.data.frame(res), row.names = NULL)
DEG <- DEG[order(DEG$padj, -abs(DEG$log2FoldChange), DEG$gene, na.last = TRUE), ]
# 完整结果保留 NA，显著基因另存。
DEG$change <- ifelse(is.na(DEG$padj), "NA",
  ifelse(DEG$log2FoldChange > 1 & DEG$padj < 0.05, "Up",
    ifelse(DEG$log2FoldChange < -1 & DEG$padj < 0.05, "Down", "NoSig")))
write.csv(DEG, file.path(out_dir, "DEG表格.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
sig_DEG <- DEG[DEG$change %in% c("Up", "Down"), ]
write.csv(sig_DEG, file.path(out_dir, "表格", "显著差异基因.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
write.csv(DEG[DEG$change == "Up", ], file.path(out_dir, "表格", "GroupX较高基因.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
write.csv(DEG[DEG$change == "Down", ], file.path(out_dir, "表格", "GroupY较高基因.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
head(DEG)
table(DEG$change)

# 核对独立过滤的实际选择；启用 independentFiltering 不保证一定排除基因。
if_meta <- metadata(res)
if_audit <- data.frame(genes = nrow(res), pvalue_NA = sum(is.na(res$pvalue)),
  padj_NA = sum(is.na(res$padj)), filterTheta = if_meta$filterTheta,
  filterThreshold = unname(if_meta$filterThreshold),
  genes_below_threshold = sum(res$baseMean < unname(if_meta$filterThreshold)))
write.csv(if_audit, file.path(out_dir, "表格", "独立过滤与NA核查.csv"), row.names = FALSE)


#### 第5步：火山图 ####
volcano_data <- DEG[is.finite(DEG$log2FoldChange) & is.finite(DEG$padj), ]
p_volcano <- ggplot(volcano_data, aes(log2FoldChange, -log10(pmax(padj, 1e-300)), color = change)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(Down = "blue", NoSig = "grey70", Up = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
  labs(title = "Volcano plot: GroupX vs GroupY",
       subtitle = "Paired DESeq2 | Up = GroupX higher; Down = GroupY higher",
       x = "log2(Fold Change)", y = "-log10(BH-adjusted P-value)") +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.title = element_blank())
save_plot(p_volcano, "Volcanoplot", 10, 7)

#### 第6步：前20个差异基因的表达热图 ####
# 教程对原始表达逐行标准化；这里先做 VST，降低测序深度和均值差异的影响。
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_mat <- assay(vsd)
sel_genes <- head(sig_DEG$gene, 20)
stopifnot(length(sel_genes) >= 20)
expr_scaled <- t(scale(t(vst_mat[sel_genes, , drop = FALSE])))
colnames(expr_scaled) <- colData$sample_label
ann_col <- data.frame(Group = colData$Group, Patient = colData$Patient,
                      row.names = colData$sample_label)
ann_colors <- list(Group = c(GroupY = "#1F77B4", GroupX = "#D62728"),
  Patient = c(P1 = "#009E73", P2 = "#CC79A7", P3 = "#E69F00", P4 = "#56B4E9"))
heatmap_plot <- pheatmap(expr_scaled, annotation_col = ann_col, annotation_colors = ann_colors,
  cluster_cols = TRUE, cluster_rows = TRUE, main = "Top Differential Genes (GroupX vs GroupY)",
  fontsize = 12, fontsize_col = 10, border_color = NA, silent = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
if (interactive()) { grid::grid.newpage(); grid::grid.draw(heatmap_plot$gtable) }
ragg::agg_png(file.path(out_dir, "Heatmap_top differential genes.png"), width = 10, height = 9, units = "in", res = 320)
grid::pushViewport(grid::viewport(width = .92, height = .94))
grid::grid.draw(heatmap_plot$gtable); grid::popViewport(); dev.off()
cairo_pdf(file.path(out_dir, "Heatmap_top differential genes.pdf"), width = 10, height = 9)
grid::pushViewport(grid::viewport(width = .92, height = .94))
grid::grid.draw(heatmap_plot$gtable); grid::popViewport(); dev.off()
write.csv(data.frame(gene = rownames(expr_scaled), expr_scaled),
          file.path(out_dir, "表格", "热图前20基因Zscore.csv"), row.names = FALSE, fileEncoding = "UTF-8")

#### 第7步：读取注释，转换基因 ID 并设置背景 ####
# 使用包内固定快照，不在运行时替换 GO 或 KEGG 版本。
db <- readRDS(file.path(input_dir, "database_bundle.rds"))
capture.output(db$metadata, file = file.path(out_dir, "数据库来源与版本.txt"))
gene_map <- unique(db$map[, c("SYMBOL", "ENTREZID")])
gene_map <- gene_map[!is.na(gene_map$SYMBOL) & !is.na(gene_map$ENTREZID), ]
gene_map <- gene_map[order(gene_map$SYMBOL, as.numeric(gene_map$ENTREZID)), ]
gene_map <- gene_map[!duplicated(gene_map$SYMBOL), ]
DEG$ENTREZID <- gene_map$ENTREZID[match(DEG$gene, gene_map$SYMBOL)]
up_genes <- DEG$gene[DEG$change == "Up"]
down_genes <- DEG$gene[DEG$change == "Down"]
gene_up <- gene_map[gene_map$SYMBOL %in% up_genes, ]
gene_down <- gene_map[gene_map$SYMBOL %in% down_genes, ]
tested_ids <- unique(na.omit(DEG$ENTREZID[is.finite(DEG$pvalue)]))
go_background <- intersect(tested_ids, unique(db$go$gene))
kegg_background <- intersect(tested_ids, unique(db$kegg$gene))
go_sets <- lapply(split(db$go$gene, db$go$term), function(x) intersect(unique(x), go_background))
kegg_sets <- lapply(split(db$kegg$gene, db$kegg$term), function(x) intersect(unique(x), kegg_background))
go_sizes <- lengths(go_sets); kegg_sizes <- lengths(kegg_sets)
go_eligible <- names(go_sizes)[go_sizes >= 10 & go_sizes <= 500]
kegg_eligible <- names(kegg_sizes)[kegg_sizes >= 10 & kegg_sizes <= 500]
mapping_summary <- data.frame(set = c("All", "Up", "Down"),
  input = c(nrow(DEG), length(up_genes), length(down_genes)),
  mapped = c(sum(!is.na(DEG$ENTREZID)), nrow(gene_up), nrow(gene_down)))
mapping_summary$percent <- 100 * mapping_summary$mapped / mapping_summary$input
write.csv(mapping_summary, file.path(out_dir, "表格", "ID映射成功率.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(ENTREZID = go_background), file.path(out_dir, "表格", "GO背景基因.csv"), row.names = FALSE)
write.csv(data.frame(ENTREZID = kegg_background), file.path(out_dir, "表格", "KEGG背景基因.csv"), row.names = FALSE)
mapping_summary

#### 第8步：GO BP 超几何富集 ####
# 教程使用 enrichGO；对已保存的 GO 注释表使用同包的 enricher，检验仍为超几何右尾检验。
ego_up <- enricher(gene = unique(gene_up$ENTREZID), universe = go_background,
  TERM2GENE = db$go, TERM2NAME = db$go_names, pAdjustMethod = "BH",
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)
ego_down <- enricher(gene = unique(gene_down$ENTREZID), universe = go_background,
  TERM2GENE = db$go, TERM2NAME = db$go_names, pAdjustMethod = "BH",
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)

#### 第9步：KEGG 超几何富集 ####
# 教程的 enrichKEGG 会在线读取 KEGG；这里用固定 KEGG 快照，以便复现同一版本。
ekegg_up <- enricher(gene = unique(gene_up$ENTREZID), universe = kegg_background,
  TERM2GENE = db$kegg, TERM2NAME = db$kegg_names, pAdjustMethod = "BH",
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)
ekegg_down <- enricher(gene = unique(gene_down$ENTREZID), universe = kegg_background,
  TERM2GENE = db$kegg, TERM2NAME = db$kegg_names, pAdjustMethod = "BH",
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1)

# 将零命中条目也计入 BH 检验家族，分别输出两个方向的完整表和显著表。
ora_tables <- list()
ora_objects <- list(GO_Up = ego_up, GO_Down = ego_down, KEGG_Up = ekegg_up, KEGG_Down = ekegg_down)
for (name in names(ora_objects)) {
  is_go <- startsWith(name, "GO")
  sizes <- if (is_go) go_sizes else kegg_sizes
  eligible <- if (is_go) go_eligible else kegg_eligible
  term_names <- if (is_go) db$go_names else db$kegg_names
  result <- ora_objects[[name]]@result
  full <- data.frame(ID = eligible, Description = term_names$name[match(eligible, term_names$term)],
                     set_size = as.integer(sizes[eligible]))
  full <- merge(full, result[, c("ID", "Count", "pvalue", "geneID", "GeneRatio", "BgRatio")],
                by = "ID", all.x = TRUE, sort = FALSE)
  full$Count[is.na(full$Count)] <- 0
  full$pvalue[is.na(full$pvalue)] <- 1
  full$p.adjust <- p.adjust(full$pvalue, method = "BH")
  full$RichFactor <- full$Count / full$set_size
  full$Group <- if (endsWith(name, "_Up")) "Upregulated" else "Downregulated"
  full <- full[order(full$p.adjust, full$pvalue, full$ID), ]
  ora_tables[[name]] <- full
  write.csv(full, file.path(out_dir, "表格", paste0(name, "_全部通路.csv")), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  write.csv(full[full$p.adjust < 0.05, ], file.path(out_dir, "表格", paste0(name, "_显著通路.csv")), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}

#### 第10步：GO 和 KEGG 气泡图 ####
# 维持教程的上下调分面布局；仅画 BH 校正 P<0.05 的条目。
for (database in c("GO", "KEGG")) {
  df_up <- ora_tables[[paste0(database, "_Up")]] %>% filter(p.adjust < 0.05) %>% head(10)
  df_down <- ora_tables[[paste0(database, "_Down")]] %>% filter(p.adjust < 0.05) %>% head(10)
  df_plot <- bind_rows(df_up, df_down)
  df_plot$Group <- factor(df_plot$Group, levels = c("Upregulated", "Downregulated"))
  df_plot$label <- vapply(df_plot$Description, function(x) paste(strwrap(x, width = 48), collapse = "\n"), character(1))
  p_ora <- ggplot(df_plot, aes(RichFactor, reorder(label, RichFactor))) +
    geom_point(aes(size = Count, color = -log10(p.adjust))) +
    facet_wrap(~ Group, scales = "free_y", ncol = 1) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_bw(base_size = 12) +
    labs(title = paste(database, "Pathway Enrichment (GroupX vs GroupY)"),
      subtitle = "Upregulated = GroupX higher; Downregulated = GroupY higher | BH-adjusted P < 0.05",
      x = "Rich Factor", y = NULL, color = "-log10(BH P)", size = "Gene Count") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          strip.text = element_text(face = "bold", size = 13), plot.margin = margin(12, 18, 12, 12))
  save_plot(p_ora, paste(database, "pathway"), 13, 13)
}

#### 第11步：全部受检基因的 GSEA ####
# 教程使用 log2FC 排名；本题采用 Wald statistic，并保留全部有限统计量。
rank_df <- DEG[is.finite(DEG$stat) & !is.na(DEG$ENTREZID), c("gene", "ENTREZID", "stat")]
rank_df <- rank_df[order(-abs(rank_df$stat), rank_df$gene, method = "radix"), ]
rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]
rank_df <- rank_df[order(-rank_df$stat, rank_df$gene, method = "radix"), ]
geneList <- setNames(rank_df$stat, rank_df$ENTREZID)
stopifnot(all(is.finite(geneList)), !anyDuplicated(names(geneList)), all(diff(geneList) <= 0))
write.csv(rank_df, file.path(out_dir, "表格", "GSEA全部基因排序.csv"), row.names = FALSE)

# 并列率使用未四舍五入的 Wald 值。112 是每组首个基因之外的重复数，
# 不是全部并列基因数；保留完整计数模式，以便核查重复的原因。
tie_groups <- split(seq_len(nrow(rank_df)), match(rank_df$stat, unique(rank_df$stat)))
tie_groups <- tie_groups[lengths(tie_groups) > 1]
tie_rows <- do.call(rbind, lapply(seq_along(tie_groups), function(i) {
  idx <- tie_groups[[i]]; genes <- rank_df$gene[idx]
  patterns <- apply(counts(dds)[genes, , drop = FALSE], 1, paste, collapse = "/")
  data.frame(tie_group = i, gene = genes, stat_17_digits = sprintf("%.17g", rank_df$stat[idx]),
    count_pattern = patterns, all_counts_identical = length(unique(patterns)) == 1)
}))
write.csv(tie_rows, file.path(out_dir, "表格", "Wald并列基因与原始count.csv"), row.names = FALSE, fileEncoding = "UTF-8")
tie_summary <- data.frame(rank_genes = nrow(rank_df), duplicate_excess = sum(duplicated(rank_df$stat)),
  duplicate_excess_percent = 100 * sum(duplicated(rank_df$stat)) / nrow(rank_df),
  tied_groups = length(tie_groups), genes_in_tied_groups = sum(lengths(tie_groups)))
write.csv(tie_summary, file.path(out_dir, "表格", "Wald并列汇总.csv"), row.names = FALSE)


# GSEA 使用 ENTREZID，映射后重新排序；固定快照通过 TERM2GENE 传入。
set.seed(23360036)
gse_go <- GSEA(geneList = geneList, TERM2GENE = db$go, TERM2NAME = db$go_names,
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH",
  eps = 0, verbose = FALSE, seed = FALSE, nPermSimple = 10000,
  BPPARAM = BiocParallel::SerialParam())
set.seed(23360036)
gse_kegg <- GSEA(geneList = geneList, TERM2GENE = db$kegg, TERM2NAME = db$kegg_names,
  minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH",
  eps = 0, verbose = FALSE, seed = FALSE, nPermSimple = 10000,
  BPPARAM = BiocParallel::SerialParam())
write.csv(gse_go@result, file.path(out_dir, "表格", "GSEA_GO_全部通路.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
write.csv(gse_kegg@result, file.path(out_dir, "表格", "GSEA_KEGG_全部通路.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
write.csv(gse_go@result[gse_go@result$p.adjust < 0.05, ], file.path(out_dir, "表格", "GSEA_GO_显著通路.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(gse_kegg@result[gse_kegg@result$p.adjust < 0.05, ], file.path(out_dir, "表格", "GSEA_KEGG_显著通路.csv"), row.names = FALSE, fileEncoding = "UTF-8")

#### 第12步：GO BP 和 KEGG 的 GSEA 条形图 ####
# 显著性标准固定为 BH P<0.05；每个方向按 BH P 升序选前10个，
# BH P 相同时按 |NES| 降序、ID 升序选择。仅改变展示，不改变检验结果。
for (database in c("GO", "KEGG")) {
  all_gsea <- if (database == "GO") gse_go@result else gse_kegg@result
  df_sig <- all_gsea %>% filter(is.finite(p.adjust), p.adjust < 0.05) %>%
    arrange(p.adjust, desc(abs(NES)), ID)
  df_up <- df_sig %>% filter(NES > 0) %>% head(10) %>% mutate(Group = "Upregulated")
  df_down <- df_sig %>% filter(NES < 0) %>% head(10) %>% mutate(Group = "Downregulated")
  df_plot <- bind_rows(df_up, df_down)
  write.csv(df_plot, file.path(out_dir, "表格", paste0("GSEA_", database, "_绘图条目.csv")),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  df_plot$Group <- factor(df_plot$Group, levels = c("Upregulated", "Downregulated"))
  df_plot$label <- vapply(paste0(df_plot$Description, " (", df_plot$ID, ")"),
    function(x) paste(strwrap(x, width = 54), collapse = "\n"), character(1))
  p_gsea <- ggplot(df_plot, aes(reorder(label, NES), NES, fill = Group)) +
    geom_col(width = 0.7) + coord_flip() +
    facet_wrap(~ Group, scales = "free_y", ncol = 1) +
    scale_fill_manual(values = c(Upregulated = "red", Downregulated = "blue")) +
    theme_bw(base_size = 11) +
    labs(title = paste(database, "GSEA Pathways (GroupX vs GroupY)"),
      subtitle = "Top 10 per direction by BH-adjusted P | BH P < 0.05 | NES > 0: GroupX; NES < 0: GroupY",
      x = NULL, y = "Normalized Enrichment Score (NES)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      strip.text = element_text(face = "bold", size = 13), legend.position = "none",
      plot.margin = margin(12, 18, 12, 12))
  save_plot(p_gsea, if (database == "GO") "GSEA pathway" else "GSEA KEGG pathway", 14, 13)
}

#### 第13步：汇总结果并最后核对条码 ####
# 图表的方向与差异基因表一致；样本类型仅用于最后的独立核对。
barcode_check <- data.frame(sample = rownames(colData), Group = colData$Group, Patient = colData$Patient)
barcode_check$sample_type <- substr(vapply(strsplit(barcode_check$sample, "-", fixed = TRUE), `[`, character(1), 4), 1, 2)
barcode_check$GDC_definition <- ifelse(barcode_check$sample_type == "01", "Primary Solid Tumor",
  ifelse(barcode_check$sample_type == "11", "Solid Tissue Normal", "Other"))
write.csv(barcode_check, file.path(out_dir, "表格", "条码独立核对.csv"), row.names = FALSE, fileEncoding = "UTF-8")
summary_table <- data.frame(analysis = c("DEG_Up", "DEG_Down", names(ora_tables), "GSEA_GO", "GSEA_KEGG"),
  significant = c(sum(DEG$change == "Up"), sum(DEG$change == "Down"),
    vapply(ora_tables, function(x) sum(x$p.adjust < 0.05), integer(1)),
    sum(gse_go@result$p.adjust < 0.05), sum(gse_kegg@result$p.adjust < 0.05)))
write.csv(summary_table, file.path(out_dir, "分析汇总.csv"), row.names = FALSE)
# 正文的结果索引按方法、数据库和 ID 对应原始结果表，不把不显著条目省略。
report_specs <- list(
  c("ORA", "KEGG", "Up", "hsa04110"), 
  c("ORA", "KEGG", "Up", "hsa04115"), 
  c("ORA", "KEGG", "Up", "hsa04218"), 
  c("GSEA", "KEGG", "", "hsa04110"), 
  c("GSEA", "GO", "", "GO:0007059"), 
  c("GSEA", "GO", "", "GO:0030336"), 
  c("GSEA", "GO", "", "GO:0002250"), 
  c("ORA", "KEGG", "Down", "hsa04610"), 
  c("ORA", "KEGG", "Down", "hsa04060"), 
  c("ORA", "GO", "Down", "GO:0019724"), 
  c("GSEA", "KEGG", "", "hsa04512"), 
  c("GSEA", "KEGG", "", "hsa04510"), 
  c("GSEA", "KEGG", "", "hsa04350"), 
  c("GSEA", "KEGG", "", "hsa01524"), 
  c("GSEA", "KEGG", "", "hsa01521"), 
  c("GSEA", "KEGG", "", "hsa01523"), 
  c("GSEA", "KEGG", "", "hsa00982"), 
  c("ORA", "KEGG", "Down", "hsa00380")
)
report_index <- do.call(rbind, lapply(report_specs, function(spec) {
  if (spec[1] == "ORA") {
    r <- ora_tables[[paste(spec[2], spec[3], sep = "_")]]
    r <- r[r$ID == spec[4], , drop = FALSE]
    direction <- if (spec[3] == "Up") "GroupX" else "GroupY"
    nes <- NA_real_
  } else {
    r <- if (spec[2] == "GO") gse_go@result else gse_kegg@result
    r <- r[r$ID == spec[4], , drop = FALSE]
    direction <- if (r$NES > 0) "GroupX" else "GroupY"
    nes <- r$NES
  }
  stopifnot(nrow(r) == 1)
  data.frame(method = spec[1], database = spec[2], ID = r$ID, Description = r$Description,
    direction = direction, BH_adjusted_P = r$p.adjust, NES = nes,
    significant = r$p.adjust < 0.05, stringsAsFactors = FALSE)
}))
write.csv(report_index, file.path(out_dir, "表格", "报告通路索引.csv"), row.names = FALSE, na = "NA", fileEncoding = "UTF-8")

saveRDS(list(dds = dds, res = res, DEG = DEG, vsd = vsd, colData = colData,
             ora = ora_tables, gse_go = gse_go, gse_kegg = gse_kegg, rank_df = rank_df),
        file.path(out_dir, "中间对象", "分析对象.rds"), compress = "xz")
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
print(summary_table)
print(barcode_check)
cat("输出目录：", out_dir, "\n")
