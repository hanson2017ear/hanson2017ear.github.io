# TCGA-LIHC 配对转录组分析
# 姓名：陳萬傑；学号：23360036
# 解压本文件夹后，在 RStudio 打开本文件并点击 Source；也可执行 Rscript "TCGA_LIHC_扩展分析.R"。
# 需要安装 R；首次运行缺少的 R 包会联网安装。data 文件夹需与本文件放在一起。
# 输出目录：TCGA_LIHC_results；输入文件和参数见 CONFIG。
# 分析环境：R 4.4.1、DESeq2 1.44.0、fgsea 1.30.0、OrgDb/GO.db 3.19.1。
# 使用已保存的数据库快照可复现原富集结果；重新下载时会记录来源和版本。

# 1. 参数设置与依赖包
if (.Platform$OS.type == "windows")
  invisible(try(suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8")), silent = TRUE))

# 从脚本的位置定位输入，支持 Rscript 和 RStudio 的 Source。
script_dir <- local({
  source_files <- vapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L) path else ""
  }, character(1))
  source_files <- source_files[nzchar(source_files) & file.exists(source_files)]
  cli_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(source_files)) {
    dirname(normalizePath(tail(source_files, 1), winslash = "/", mustWork = TRUE))
  } else if (length(cli_file)) {
    dirname(normalizePath(sub("^--file=", "", cli_file[1]), winslash = "/", mustWork = TRUE))
  } else {
    getwd()
  }
})

CONFIG <- list(
  project_dir = Sys.getenv("TCGA_PROJECT_DIR", script_dir),
  input_dir = "data",
  count_file = "TCGA_LIHC的8个样本count.csv",
  group_file = "TCGA_LIHC的8个样本分组.csv",
  output_dir = Sys.getenv("TCGA_OUTPUT_DIR", "TCGA_LIHC_results"),
  database_snapshot = "data/database_bundle.rds",
  existing_library = Sys.getenv("TCGA_EXISTING_R_LIB", ""),
  # 优先使用现有 R 库，仅安装缺失的包。
  install_missing = TRUE,
  seed = 23360036L,
  integer_tolerance = 1e-6,
  min_count = 10L,
  min_samples = 4L,
  padj_cutoff = 0.05,
  lfc_cutoff = 1,
  min_gene_set_size = 10L,
  max_gene_set_size = 500L
)

resolve_path <- function(path, root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)) path else file.path(root, path)
}

ensure_packages <- function(cfg, out, need_database_packages = FALSE) {
  library_dir <- file.path(out, "R_library")
  dir.create(library_dir, recursive = TRUE, showWarnings = FALSE)
  r_minor <- paste(R.version$major, strsplit(R.version$minor, ".", fixed = TRUE)[[1]][1], sep = ".")
  # Windows 下同时查找常规用户库和 R_LIBS_USER 指定的库。
  windows_user_libraries <- if (.Platform$OS.type == "windows") c(
    file.path(Sys.getenv("LOCALAPPDATA"), "R", "win-library", r_minor),
    file.path(Sys.getenv("USERPROFILE"), "AppData", "Local", "R", "win-library", r_minor)) else character()
  candidates <- c(cfg$existing_library, file.path(cfg$project_dir, "作业成果", "R_library"),
    file.path(cfg$project_dir, "作业成果", "Rlib"), windows_user_libraries, .libPaths(), library_dir)
  .libPaths(unique(candidates[nzchar(candidates) & dir.exists(candidates)]))
  cat("实际 R 库：\n"); print(.libPaths())
  packages <- c("DESeq2", "ggplot2", "pheatmap", "ggrepel", "jsonlite", "ragg", "fgsea")
  if (need_database_packages) packages <- c(packages, "AnnotationDbi", "org.Hs.eg.db", "GO.db")
  # R 4.4 对应原分析使用的 Bioconductor 3.19；其他 R 版本使用匹配的源。
  if (getRversion() >= "4.4" && getRversion() < "4.5") {
    repos <- c(BioC = "https://bioconductor.posit.co/packages/3.19/bioc",
      Annotation = "https://bioconductor.posit.co/packages/3.19/data/annotation",
      CRAN = "https://cloud.r-project.org")
  } else {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      if (!cfg$install_missing) stop("缺少 R 包：", package, "。请安装后重新执行。")
      cat("安装缺少的包：", package, "\n")
      tryCatch(install.packages(package, lib = library_dir, repos = repos,
        dependencies = NA, type = if (.Platform$OS.type == "windows") "binary" else "source"),
        error = function(e) message("直接安装未完成：", conditionMessage(e)))
      if (!requireNamespace(package, quietly = TRUE)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager", lib = library_dir, repos = "https://cloud.r-project.org")
        if (getRversion() >= "4.4" && getRversion() < "4.5")
          BiocManager::install(version = "3.19", ask = FALSE, update = FALSE, lib = library_dir)
        BiocManager::install(package, ask = FALSE, update = FALSE, lib = library_dir)
      }
    }
    if (!requireNamespace(package, quietly = TRUE)) stop("依赖安装仍未完成：", package)
    cat(package, as.character(packageVersion(package)), "\n")
  }
  suppressPackageStartupMessages({
    library(DESeq2); library(ggplot2); library(pheatmap); library(ggrepel)
    library(jsonlite); library(ragg); library(fgsea)
  })
}

run_differential <- function(cfg, out) {
  write_csv <- function(x, name) write.csv(x, file.path(out, "结果", paste0(name, ".csv")),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  set.seed(cfg$seed)
  input <- resolve_path(cfg$input_dir, cfg$project_dir)
  countfile <- file.path(input, cfg$count_file)
  groupfile <- file.path(input, cfg$group_file)
  if (!all(file.exists(c(countfile, groupfile)))) stop("找不到输入 CSV，请检查 CONFIG$input_dir。")
  # 原样复制输入；分析只写入 output_dir。
  copies <- file.path(out, "原始输入", basename(c(countfile, groupfile)))
  for (i in seq_along(copies)) {
    original <- c(countfile, groupfile)[i]
    if (normalizePath(original, winslash = "/") != normalizePath(copies[i], winslash = "/", mustWork = FALSE))
      if (!file.copy(original, copies[i], overwrite = TRUE)) stop("无法保留原始输入：", original)
  }
  write_csv(data.frame(file = c(countfile, groupfile), md5 = as.character(tools::md5sum(c(countfile, groupfile)))), "input_md5")
  # 2. 输入检查与患者配对
  raw <- read.csv(countfile,row.names=1,check.names=FALSE)
  meta <- read.csv(groupfile,row.names=1,check.names=FALSE)
  stopifnot(!anyDuplicated(rownames(raw)),!anyDuplicated(colnames(raw)),!anyDuplicated(rownames(meta)),setequal(colnames(raw),rownames(meta)))
  original_order_identical <- identical(colnames(raw),rownames(meta))
  meta <- meta[colnames(raw),,drop=FALSE]
  stopifnot(identical(colnames(raw),rownames(meta)),all(meta$Group %in% c("GroupX","GroupY")))
  cts <- as.matrix(raw)
  stopifnot(is.numeric(cts),all(is.finite(cts)),all(cts>=0),all(nzchar(rownames(cts))))
  # 条码前三段用于患者配对；样本类型在表达分析后核对。
  patient_full <- vapply(strsplit(rownames(meta),"-",fixed=TRUE),function(x) paste(x[1:3],collapse="-"),character(1))
  meta$patient <- factor(patient_full,levels=sort(unique(patient_full)))
  meta$Group <- factor(meta$Group,levels=c("GroupY","GroupX"))
  stopifnot(all(table(meta$patient,meta$Group)==1),nlevels(meta$patient)==4,nrow(meta)==8)
  meta$patient_label <- paste0("P",as.integer(meta$patient))
  meta$sample_label <- paste0(meta$patient_label,"_",ifelse(meta$Group=="GroupX","X","Y"))
  write_csv(data.frame(sample=rownames(meta),meta),"sample_metadata")
  # 含非整数 count 的基因另存清单，并从主分析中排除。
  integer_error <- abs(cts-round(cts)); bad <- rowSums(integer_error>cfg$integer_tolerance)>0
  write_csv(data.frame(gene=rownames(cts)[bad],cts[bad,,drop=FALSE],max_integer_deviation=apply(integer_error[bad,,drop=FALSE],1,max)),"excluded_fractional_genes")
  near_cells <- which(integer_error>0 & integer_error<=cfg$integer_tolerance,arr.ind=TRUE)
  write_csv(data.frame(gene=rownames(cts)[near_cells[,1]],sample=colnames(cts)[near_cells[,2]],original=cts[near_cells],rounded=round(cts[near_cells]),absolute_error=integer_error[near_cells]),"near_integer_repair")
  # 低表达过滤：至少 4/8 个样本的 count >=10。
  keep_low <- rowSums(round(cts)>=cfg$min_count)>=cfg$min_samples
  keep <- !bad & keep_low
  intcts <- round(cts[keep,,drop=FALSE]); stopifnot(all(intcts <= .Machine$integer.max))
  storage.mode(intcts)<-"integer"
  status <- ifelse(bad,"excluded_fractional",ifelse(!keep_low,"low_expression","modeled"))
  write_csv(data.frame(gene=rownames(cts),status=status),"gene_filter_status")
  qc <- list(input_genes=nrow(cts),samples=ncol(cts),negative=sum(cts<0),missing=sum(is.na(cts)),duplicate_ids=anyDuplicated(rownames(cts)),original_order_identical=original_order_identical,
   fractional_cells=sum(integer_error>cfg$integer_tolerance),fractional_genes=sum(bad),near_integer_cells=nrow(near_cells),max_integer_deviation=max(integer_error),near_max=max(integer_error[integer_error<=cfg$integer_tolerance]),
   low_expression_excluded=sum(!bad & !keep_low),modeled_genes=sum(keep),filter="count >=10 in at least 4 of 8 samples",integer_tolerance=cfg$integer_tolerance,
   design="~ patient + Group",contrast="GroupX versus GroupY",seed=cfg$seed)
  print(qc)

  # 3. 配对差异表达分析
  dds <- DESeqDataSetFromMatrix(intcts,meta,design=~patient+Group)
  dds <- DESeq(dds,quiet=FALSE,parallel=FALSE)
  res <- results(dds,contrast=c("Group","GroupX","GroupY"),alpha=cfg$padj_cutoff,pAdjustMethod="BH",independentFiltering=TRUE)
  de <- data.frame(gene=rownames(res),as.data.frame(res),row.names=NULL)
  de$direction <- ifelse(is.na(de$padj),"NA",ifelse(de$padj<cfg$padj_cutoff & de$log2FoldChange>cfg$lfc_cutoff,"GroupX_higher",ifelse(de$padj<cfg$padj_cutoff & de$log2FoldChange< -cfg$lfc_cutoff,"GroupY_higher","Not_significant")))
  de <- de[order(de$padj,-abs(de$log2FoldChange),de$gene,na.last=TRUE),]
  sig <- de[de$direction %in% c("GroupX_higher","GroupY_higher"),]
  write_csv(de,"DESeq2_all_modeled_genes")
  write_csv(sig,"DEG_significant")
  write_csv(sig[sig$direction=="GroupX_higher",],"DEG_GroupX_higher")
  write_csv(sig[sig$direction=="GroupY_higher",],"DEG_GroupY_higher")
  allinput <- merge(data.frame(gene=rownames(cts),input_status=status),de,by="gene",all.x=TRUE,sort=FALSE)
  write_csv(allinput,"all_input_genes_with_status")
  qc$finite_pvalue <- sum(is.finite(de$pvalue)); qc$padj_NA <- sum(is.na(de$padj)); qc$pvalue_NA <- sum(is.na(de$pvalue))
  qc$up_GroupX <- sum(sig$direction=="GroupX_higher"); qc$up_GroupY <- sum(sig$direction=="GroupY_higher")
  qc$size_factors <- setNames(as.list(sizeFactors(dds)),meta$sample_label)
  qc$independent_filter_threshold <- unname(metadata(res)$filterThreshold)
  qc$cooks_flagged_genes <- sum(apply(assays(dds)[["cooks"]],1,max,na.rm=TRUE)>qf(.99,ncol(model.matrix(design(dds),colData(dds))),ncol(dds)-ncol(model.matrix(design(dds),colData(dds)))))
  write_csv(data.frame(gene=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),"normalized_counts")

  # 4. VST 转换与表达图表
  # VST仅用于可视化；差异模型始终使用整数count。
  vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
  vmat <- assay(vsd); colnames(vmat)<-meta$sample_label
  write_csv(data.frame(gene=rownames(vmat),vmat),"VST_expression")
  saveRDS(list(dds=dds,res=res,de=de,sig=sig,vsd=vsd,meta=meta,qc=qc,status=status),file.path(out,"中间对象/differential_analysis.rds"),compress="xz")
  save_plot <- function(p,name,w=9,h=6){
   ggsave(file.path(out,"图表",paste0(name,".png")),p,width=w,height=h,dpi=320,device=ragg::agg_png,bg="white")
   ggsave(file.path(out,"图表",paste0(name,".pdf")),p,width=w,height=h,device=cairo_pdf,bg="white")
  }
  theme_set(theme_bw(base_size=12)+theme(plot.title=element_text(face="bold"),plot.subtitle=element_text(size=10),plot.margin=margin(12,15,12,12)))
  cols <- c(GroupX="#D55E00",GroupY="#0072B2")
  # 火山图只绘制padj和效应均有限的点，完整结果表保留全部NA。
  vp <- de[is.finite(de$padj)&is.finite(de$log2FoldChange),]
  vp$display <- factor(ifelse(vp$direction=="GroupX_higher","GroupX higher",ifelse(vp$direction=="GroupY_higher","GroupY higher","Not significant")),levels=c("GroupX higher","GroupY higher","Not significant"))
  lab <- unique(c(head(sig$gene,8),intersect(c("GPC3","MET","CCND1","MYC","EPCAM","MKI67","ALB","CYP3A4"),sig$gene)))
  p <- ggplot(vp,aes(log2FoldChange,-log10(pmax(padj,1e-300)),color=display))+geom_point(size=.9,alpha=.55)+
   geom_vline(xintercept=c(-1,1),linetype=2,color="grey50")+geom_hline(yintercept=-log10(.05),linetype=2,color="grey50")+
   ggrepel::geom_text_repel(data=vp[vp$gene %in% lab,],aes(label=gene),size=3,max.overlaps=25,seed=cfg$seed,show.legend=FALSE)+
   scale_color_manual(values=c("GroupX higher"=unname(cols["GroupX"]),"GroupY higher"=unname(cols["GroupY"]),"Not significant"="grey75"))+
   labs(title="Differential expression: GroupX versus GroupY",subtitle=sprintf("Paired DESeq2 | BH padj < 0.05 and |log2FC| > 1 | X higher: %s; Y higher: %s",qc$up_GroupX,qc$up_GroupY),x="log2 fold change (positive = GroupX higher)",y="-log10(BH-adjusted P value)",color=NULL)+theme(legend.position="bottom")
  save_plot(p,"01_volcano",10,7)
  # PCA采用VST中变异度最高的500基因，显示同病人的配对连线。
  rv <- apply(vmat,1,var); pca <- prcomp(t(vmat[head(order(rv,decreasing=TRUE),500),]),center=TRUE,scale.=FALSE)
  pv <- pca$sdev^2/sum(pca$sdev^2)
  pc <- data.frame(sample=meta$sample_label,Group=meta$Group,patient=meta$patient_label,PC1=pca$x[,1],PC2=pca$x[,2])
  write_csv(pc,"PCA_coordinates"); qc$PCA_variance <- as.list(pv[1:2])
  p <- ggplot(pc,aes(PC1,PC2,color=Group,shape=patient))+geom_line(aes(group=patient),color="grey70",linewidth=.5)+geom_point(size=4)+
   ggrepel::geom_text_repel(aes(label=sample),size=3,seed=cfg$seed)+scale_color_manual(values=cols)+
   labs(title="PCA of VST expression",subtitle="Top 500 variable genes; lines connect matched patients",x=sprintf("PC1 (%.1f%%)",100*pv[1]),y=sprintf("PC2 (%.1f%%)",100*pv[2]))
  save_plot(p,"02_PCA",8.5,6)
  ann <- data.frame(Group=meta$Group,Patient=meta$patient_label,row.names=meta$sample_label)
  anncol <- list(Group=cols,Patient=setNames(c("#009E73","#CC79A7","#E69F00","#56B4E9"),paste0("P",1:4)))
  save_heat <- function(mat,name,title,w,h,scale="none",...){
   ph <- pheatmap(mat,annotation_col=ann,annotation_colors=anncol,main=title,scale=scale,silent=TRUE,fontsize=10,border_color=NA,...)
   ragg::agg_png(file.path(out,"图表",paste0(name,".png")),width=w,height=h,units="in",res=320); grid::pushViewport(grid::viewport(width=.94,height=.96)); grid::grid.draw(ph$gtable); dev.off()
   cairo_pdf(file.path(out,"图表",paste0(name,".pdf")),width=w,height=h); grid::pushViewport(grid::viewport(width=.94,height=.96)); grid::grid.draw(ph$gtable); dev.off()
  }
  cormat <- cor(vmat,method="pearson"); write_csv(data.frame(sample=rownames(cormat),cormat),"VST_sample_correlation")
  save_heat(cormat,"03_sample_correlation","Sample correlation (Pearson, all VST genes)",8.5,7,display_numbers=TRUE,number_format="%.3f",color=colorRampPalette(c("white","#0072B2"))(100))
  top <- head(sig$gene,40); zmat <- t(scale(t(vmat[top,,drop=FALSE])))
  stopifnot(nrow(zmat)>=20,all(is.finite(zmat)))
  write_csv(data.frame(gene=rownames(zmat),zmat),"heatmap_top40_rowZ")
  save_heat(zmat,"04_DEG_heatmap","Top 40 DE genes (VST, row Z score)",9,11,color=colorRampPalette(c("#0072B2","white","#D55E00"))(101))
  saveRDS(list(dds=dds,res=res,de=de,sig=sig,vsd=vsd,meta=meta,qc=qc,status=status),file.path(out,"中间对象/differential_analysis.rds"),compress="xz")
  jsonlite::write_json(qc,file.path(out,"结果/analysis_summary.json"),auto_unbox=TRUE,pretty=TRUE,na="null")
  # 将非整数值取整后重跑，检查其对主分析结果的影响。
  senscts <- round(cts[keep_low,,drop=FALSE]); storage.mode(senscts)<-"integer"
  sensdds <- DESeq(DESeqDataSetFromMatrix(senscts,meta,~patient+Group),quiet=TRUE,parallel=FALSE)
  sr <- as.data.frame(results(sensdds,contrast=c("Group","GroupX","GroupY"),alpha=cfg$padj_cutoff,pAdjustMethod="BH"))
  sr$gene <- rownames(sr); rownames(sr)<-NULL
  sr$significant <- !is.na(sr$padj)&sr$padj<cfg$padj_cutoff&abs(sr$log2FoldChange)>cfg$lfc_cutoff
  write_csv(sr,"sensitivity_rounded_counts_DESeq2")
  common <- merge(de,sr,by="gene",suffixes=c("_main","_rounded"))
  sig_main <- sig$gene; sig_sens <- sr$gene[sr$significant]
  sens <- data.frame(main_DEG=nrow(sig),rounded_DEG=sum(sr$significant),overlap=length(intersect(sig_main,sig_sens)),jaccard=length(intersect(sig_main,sig_sens))/length(union(sig_main,sig_sens)),log2FC_pearson=cor(common$log2FoldChange_main,common$log2FoldChange_rounded,use="complete.obs"),direction_discordant=sum(sign(common$log2FoldChange_main)!=sign(common$log2FoldChange_rounded),na.rm=TRUE))
  write_csv(sens,"sensitivity_summary"); print(sens)

  # 检查结果表的基因数、比较方向和显著基因集合。
  stopifnot(nrow(sig) == qc$up_GroupX + qc$up_GroupY,
    nrow(allinput) == nrow(cts), all(!(rownames(intcts) %in% rownames(cts)[bad])),
    !anyNA(intcts), all(intcts >= 0), all(intcts == round(intcts)))
  cat("差异分析完成：", nrow(de), "个建模基因，", nrow(sig), "个显著基因。\n")
  list(dds = dds, res = res, de = de, sig = sig, vsd = vsd, meta = meta, qc = qc,
    status = status)
}

# 5. 数据库快照与基因映射
load_or_build_database <- function(cfg, analysis_out) {
  out <- file.path(analysis_out, "数据库")
  snapshot <- resolve_path(cfg$database_snapshot, cfg$project_dir)
  if (nzchar(cfg$database_snapshot) && file.exists(snapshot)) {
    cat("使用固定数据库快照：", normalizePath(snapshot, winslash = "/"), "\n")
    bundle <- readRDS(snapshot)
    stopifnot(all(c("map", "go", "go_names", "kegg", "kegg_names", "metadata") %in% names(bundle)),
      all(c("SYMBOL", "ENTREZID") %in% names(bundle$map)),
      all(c("term", "gene") %in% names(bundle$go)),
      all(c("term", "gene") %in% names(bundle$kegg)))
    saveRDS(bundle, file.path(out, "database_bundle.rds"), compress = "xz")
    capture.output(bundle$metadata, file = file.path(out, "数据库来源与版本.txt"))
    writeLines(paste("Snapshot MD5:", tools::md5sum(snapshot)), file.path(out, "snapshot_md5.txt"))
    for (nm in c("map", "go", "go_names", "kegg", "kegg_names"))
      write.table(bundle[[nm]], file.path(out, paste0(nm, ".tsv")), sep = "\t",
        quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
    return(bundle)
  }
  cat("未找到固定快照；以现有 OrgDb/GO.db 和官方 KEGG REST 建立数据库。\n")
  suppressPackageStartupMessages({library(AnnotationDbi); library(org.Hs.eg.db); library(GO.db)})
  # 一对多映射完整保留，实际分析再按预先指定的规则去重。
  map <- unique(AnnotationDbi::select(org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "SYMBOL"),
    keytype = "SYMBOL", columns = "ENTREZID"))
  map <- map[!is.na(map$SYMBOL) & !is.na(map$ENTREZID), c("SYMBOL", "ENTREZID")]
  map <- map[order(map$SYMBOL, as.numeric(map$ENTREZID)), ]
  rownames(map) <- NULL

  # GOALL 已由 OrgDb 将直接注释向祖先术语传播；仅保留 Biological Process。
  g <- AnnotationDbi::select(org.Hs.eg.db,
    keys = AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"),
    keytype = "ENTREZID", columns = c("GOALL", "ONTOLOGYALL"))
  g <- g[!is.na(g$GOALL) & !is.na(g$ONTOLOGYALL) & g$ONTOLOGYALL == "BP", ]
  go <- unique(data.frame(term = g$GOALL, gene = g$ENTREZID))
  go <- go[order(go$term, as.numeric(go$gene)), ]
  gn <- AnnotationDbi::select(GO.db, keys = unique(go$term),
    keytype = "GOID", columns = c("TERM", "ONTOLOGY"))
  go_names <- unique(data.frame(term = gn$GOID, name = gn$TERM))
  go_names <- go_names[!is.na(go_names$name), ]
  go <- go[go$term %in% go_names$term, ]
  rownames(go) <- rownames(go_names) <- NULL

  kegg_urls <- c(link = "https://rest.kegg.jp/link/hsa/pathway",
    names = "https://rest.kegg.jp/list/pathway/hsa")
  kegg_files <- file.path(out, c("kegg_link_hsa_pathway.tsv", "kegg_list_pathway_hsa.tsv"))
  for (i in seq_along(kegg_files)) {
    if (!file.exists(kegg_files[i]) || file.info(kegg_files[i])$size == 0) {
      options(timeout = max(600, getOption("timeout")))
      download.file(kegg_urls[i], kegg_files[i], mode = "wb", method = "libcurl")
    }
  }
  kl <- read.delim(kegg_files[1], header = FALSE, sep = "\t", quote = "", colClasses = "character")
  kn <- read.delim(kegg_files[2], header = FALSE, sep = "\t", quote = "", colClasses = "character")
  stopifnot(ncol(kl) == 2L, ncol(kn) == 2L,
    all(grepl("^path:hsa", kl[[1]])), all(grepl("^hsa:", kl[[2]])))
  kegg <- unique(data.frame(term = sub("^path:", "", kl[[1]]), gene = sub("^hsa:", "", kl[[2]])))
  kegg_names <- unique(data.frame(term = sub("^path:", "", kn[[1]]),
    name = sub(" - Homo sapiens \\(human\\)$", "", kn[[2]])))
  kegg <- kegg[kegg$term %in% kegg_names$term, ]
  metadata <- list(
    created = format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE),
    org_Hs_eg_db_version = as.character(packageVersion("org.Hs.eg.db")),
    GO_db_version = as.character(packageVersion("GO.db")),
    org_Hs_eg_db_metadata = AnnotationDbi::metadata(org.Hs.eg.db),
    GO_db_metadata = AnnotationDbi::metadata(GO.db),
    go_annotation = "OrgDb GOALL BP including ancestor propagation; unique term-gene pairs",
    mapping = "Official SYMBOL only, no ALIAS; preserve all one-to-many pairs for downstream deduplication",
    kegg_source_urls = kegg_urls,
    kegg_download_date = format(as.Date(file.info(kegg_files)$mtime, tz = "Asia/Shanghai")),
    kegg_raw_md5 = tools::md5sum(kegg_files))
  bundle <- list(map = map, go = go, go_names = go_names,
    kegg = kegg, kegg_names = kegg_names, metadata = metadata)
  saveRDS(bundle, file.path(out, "database_bundle.rds"), compress = "xz")
  for (nm in c("map", "go", "go_names", "kegg", "kegg_names"))
    write.table(bundle[[nm]], file.path(out, paste0(nm, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
  dput(metadata, file = file.path(out, "database_metadata.R"))
  capture.output(metadata, file = file.path(out, "数据库来源与版本.txt"))
  cat("行数：\n"); print(vapply(bundle[1:5], nrow, integer(1)))
  bundle
}

run_enrichment <- function(s, b, cfg, out) {
  write_csv <- function(x, name) write.csv(x, file.path(out, "结果", name),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  save_plot <- function(p, name, width = 13, height = 10) {
    ggsave(file.path(out, "图表", paste0(name, ".png")), p, width = width,
      height = height, dpi = 320, bg = "white", limitsize = FALSE)
    ggsave(file.path(out, "图表", paste0(name, ".pdf")), p, width = width,
      height = height, device = cairo_pdf, bg = "white", limitsize = FALSE)
  }
  wrap_text <- function(x, width = 56) vapply(x,
    function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
  pal <- c(GroupX_higher = "#D55E00", GroupY_higher = "#0072B2")
  dir_labels <- c(GroupX_higher = "GroupX higher (log2FC > 1)",
    GroupY_higher = "GroupY higher (log2FC < -1)")
  theme_set(theme_bw(base_size = 11) + theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10), panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 12), legend.position = "bottom",
    plot.margin = margin(12, 18, 12, 12)))

  de <- as.data.frame(s$de)
  stopifnot(all(c("gene", "baseMean", "log2FoldChange", "stat", "pvalue", "padj") %in% names(de)))
  de$gene <- as.character(de$gene)
  stopifnot(!anyNA(de$gene), !anyDuplicated(de$gene))
  # SYMBOL 一对多时，取数值最小的 ENTREZID。

  map_all <- b$map
  map_all <- map_all[!is.na(map_all$SYMBOL) & !is.na(map_all$ENTREZID), ]
  map_all <- unique(map_all)
  map_all <- map_all[order(map_all$SYMBOL, as.numeric(map_all$ENTREZID)), ]
  mapping_ambiguities <- map_all[duplicated(map_all$SYMBOL) |
    duplicated(map_all$SYMBOL, fromLast = TRUE), ]
  write_csv(mapping_ambiguities, "ID_mapping_one_to_many_all_candidates.csv")
  map_one <- map_all[!duplicated(map_all$SYMBOL), ]
  de$ENTREZID <- map_one$ENTREZID[match(de$gene, map_one$SYMBOL)]
  de$GO_BP_annotated <- !is.na(de$ENTREZID) & de$ENTREZID %in% b$go$gene
  de$KEGG_annotated <- !is.na(de$ENTREZID) & de$ENTREZID %in% b$kegg$gene
  de$finite_pvalue <- is.finite(de$pvalue)
  de$finite_Wald <- is.finite(de$stat)
  de$is_DEG <- !is.na(de$padj) & de$padj < cfg$padj_cutoff & abs(de$log2FoldChange) > cfg$lfc_cutoff
  de$enrichment_direction <- ifelse(de$is_DEG & de$log2FoldChange > cfg$lfc_cutoff,
    "GroupX_higher", ifelse(de$is_DEG & de$log2FoldChange < -cfg$lfc_cutoff,
    "GroupY_higher", "Not_in_DEG_foreground"))
  write_csv(de[, c("gene", "ENTREZID", "GO_BP_annotated", "KEGG_annotated",
    "finite_pvalue", "finite_Wald", "is_DEG", "enrichment_direction")],
    "ID_mapping_for_analysis_genes.csv")
  write_csv(de[is.na(de$ENTREZID), c("gene", "baseMean", "stat", "pvalue", "padj")],
    "ID_mapping_unmapped_symbols.csv")

  # 映射率同时报告 SYMBOL 成功率与去重后的 ENTREZID 数目。
  idx_sets <- list(All_model_genes = rep(TRUE, nrow(de)),
    Tested_finite_pvalue = de$finite_pvalue,
    Ranked_finite_Wald = de$finite_Wald,
    GroupX_higher_DEGs = de$enrichment_direction == "GroupX_higher",
    GroupY_higher_DEGs = de$enrichment_direction == "GroupY_higher")
  map_summary <- do.call(rbind, lapply(names(idx_sets), function(nm) {
    z <- de[idx_sets[[nm]], ]; n <- nrow(z)
    data.frame(set = nm, input_symbols = n,
      mapped_symbols = sum(!is.na(z$ENTREZID)),
      mapping_success_percent = if(n) 100 * mean(!is.na(z$ENTREZID)) else NA_real_,
      unique_Entrez = length(unique(na.omit(z$ENTREZID))),
      GO_BP_mapped_symbols = sum(z$GO_BP_annotated),
      GO_BP_unique_Entrez = length(unique(z$ENTREZID[z$GO_BP_annotated])),
      KEGG_mapped_symbols = sum(z$KEGG_annotated),
      KEGG_unique_Entrez = length(unique(z$ENTREZID[z$KEGG_annotated])))
  }))
  write_csv(map_summary, "ID_mapping_summary.csv")
  symbol_by_entrez <- split(de$gene[!is.na(de$ENTREZID)], de$ENTREZID[!is.na(de$ENTREZID)])
  entrez_to_symbols <- function(ids) paste(unique(unlist(symbol_by_entrez[as.character(ids)])), collapse = ";")
  databases <- list(GO_BP = list(t2g = b$go, names = b$go_names),
    KEGG = list(t2g = b$kegg, names = b$kegg_names))
  min_size <- cfg$min_gene_set_size
  max_size <- cfg$max_gene_set_size
  seed <- cfg$seed
  ora_all <- list(); bg_summary <- list(); eligible_sets <- list()


  # 6. GO 和 KEGG 超几何富集
  # 超几何检验：P[X >= k]。BH 在每个数据库/方向全部10–500背景内基因集合上校正，
  # 零命中集合的 P=1，也计入 BH 校正。
  for (db in names(databases)) {
    dbi <- databases[[db]]
    universe <- unique(na.omit(de$ENTREZID[de$finite_pvalue]))
    universe <- sort(intersect(universe, unique(dbi$t2g$gene)))
    sets <- lapply(split(dbi$t2g$gene, dbi$t2g$term), function(z) intersect(unique(z), universe))
    sizes <- lengths(sets)
    sets <- sets[sizes >= min_size & sizes <= max_size]
    eligible_sets[[db]] <- sets
    bg_summary[[db]] <- data.frame(database = db, universe_Entrez = length(universe),
      eligible_terms = length(sets), minGSSize = min_size, maxGSSize = max_size)
    write_csv(data.frame(ENTREZID = universe,
      symbols = vapply(universe, function(z) entrez_to_symbols(z), character(1))),
      paste0(db, "_ORA_background.csv"))
    write_csv(data.frame(ID = names(sets), size_in_background = lengths(sets)),
      paste0(db, "_ORA_tested_term_sizes.csv"))
    for (direction in names(pal)) {
      foreground <- unique(na.omit(de$ENTREZID[de$enrichment_direction == direction]))
      foreground <- intersect(foreground, universe)
      N <- length(universe); n <- length(foreground)
      hits <- lapply(sets, intersect, y = foreground)
      k <- lengths(hits); M <- lengths(sets)
      pv <- if (n > 0 && N > 0) phyper(k - 1L, M, N - M, n, lower.tail = FALSE) else rep(1, length(M))
      result <- data.frame(ID = names(sets),
        Description = dbi$names$name[match(names(sets), dbi$names$term)],
        direction = direction, Count = k,
        GeneRatio = paste0(k, "/", n), BgRatio = paste0(M, "/", N),
        GeneRatio_numeric = if(n) k / n else NA_real_,
        RichFactor = k / M,
        FoldEnrichment = if(n) (k / n) / (M / N) else NA_real_,
        pvalue = pv, p.adjust = p.adjust(pv, method = "BH"),
        foreground_size = n, background_size = N, set_size = M,
        geneID = vapply(hits, paste, collapse = ";", character(1)),
        symbols = vapply(hits, entrez_to_symbols, character(1)))
      result <- result[order(result$p.adjust, result$pvalue, result$ID), ]
      rownames(result) <- NULL
      nm <- paste(db, direction, sep = "_")
      ora_all[[nm]] <- result
      write_csv(result, paste0("ORA_", nm, "_all.csv"))
      write_csv(result[!is.na(result$p.adjust) & result$p.adjust < 0.05, ],
        paste0("ORA_", nm, "_significant.csv"))
      cat("ORA", db, direction, "background", N, "foreground", n,
        "terms", nrow(result), "BH<0.05", sum(result$p.adjust < 0.05, na.rm = TRUE), "\n")
    }
  }
  write_csv(do.call(rbind, bg_summary), "ORA_background_summary.csv")


  # 7. 受检基因的 GSEA
  # GSEA：全部有限 Wald 统计量（包括 padj 为 NA 者）；多个 SYMBOL 对应同一 Entrez 时，
  # 取 |stat| 最大者，再按 gene 字典序打破同值。映射完成后重新按 stat 降序排列。
  rank_df <- de[de$finite_Wald & !is.na(de$ENTREZID), c("gene", "ENTREZID", "stat", "padj")]
  rank_df <- rank_df[order(-abs(rank_df$stat), rank_df$gene, method = "radix"), ]
  rank_excluded <- rank_df[duplicated(rank_df$ENTREZID), ]
  rank_df <- rank_df[!duplicated(rank_df$ENTREZID), ]
  rank_df <- rank_df[order(-rank_df$stat, rank_df$gene, method = "radix"), ]
  rownames(rank_df) <- NULL
  ranks <- setNames(rank_df$stat, rank_df$ENTREZID)
  stopifnot(!anyDuplicated(names(ranks)), all(is.finite(ranks)),
    all(diff(ranks) <= 0))
  write_csv(rank_df, "GSEA_ranked_genes_Wald.csv")
  write_csv(rank_excluded, "GSEA_duplicate_Entrez_discarded_rows.csv")
  cat("GSEA ranked unique Entrez:", length(ranks), "; tied-stat fraction:",
    sum(duplicated(ranks)) / length(ranks), "\n")
  gsea_all <- list(); gsea_sets <- list()
  for (db in names(databases)) {
    dbi <- databases[[db]]
    pathways <- lapply(split(dbi$t2g$gene, dbi$t2g$term), unique)
    gsea_sets[[db]] <- pathways
    set.seed(seed)
    cat("GSEA start", db, format(Sys.time(), tz = "Asia/Shanghai"), "\n")
    f <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks,
      minSize = min_size, maxSize = max_size, eps = 0,
      nPermSimple = 10000, nproc = 1, scoreType = "std")
    r <- as.data.frame(f)
    r$ID <- r$pathway
    r$Description <- dbi$names$name[match(r$ID, dbi$names$term)]
    r$direction <- ifelse(r$NES > 0, "GroupX_higher", "GroupY_higher")
    r$leadingEdge_symbols <- vapply(r$leadingEdge, entrez_to_symbols, character(1))
    r$leadingEdge <- vapply(r$leadingEdge, paste, collapse = ";", character(1))
    r <- r[, c("ID", "Description", "direction", "size", "ES", "NES", "pval", "padj",
      "log2err", "leadingEdge", "leadingEdge_symbols")]
    r <- r[order(r$padj, -abs(r$NES), r$ID), ]
    rownames(r) <- NULL
    gsea_all[[db]] <- r
    write_csv(r, paste0("GSEA_", db, "_all.csv"))
    write_csv(r[!is.na(r$padj) & r$padj < 0.05, ], paste0("GSEA_", db, "_significant.csv"))
    cat("GSEA complete", db, "terms", nrow(r), "BH<0.05", sum(r$padj < 0.05, na.rm = TRUE), "\n")
  }


  # 8. 富集图表与结果核对
  # 富集图显示校正后显著的条目；没有显著结果时输出说明页。
  empty_plot <- function(title, explanation) ggplot() +
    annotate("text", x = 0, y = 0, label = explanation, size = 5, lineheight = 1.2) +
    labs(title = title, subtitle = "Fixed cutoff: BH-adjusted P < 0.05") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void(base_size = 12)
  for (db in names(databases)) {
    sig <- do.call(rbind, lapply(names(pal), function(d) {
      z <- ora_all[[paste(db, d, sep = "_")]]
      head(z[!is.na(z$p.adjust) & z$p.adjust < 0.05, ], 8)
    }))
    write_csv(sig, paste0(db, "_ORA_plot_data.csv"))
    title <- paste(db, "over-representation analysis")
    if (!nrow(sig)) {
      p <- empty_plot(title, "No terms pass BH-adjusted P < 0.05\nin either higher-expression gene set.")
      save_plot(p, paste0(db, "_ORA_bubble"), height = 5)
      save_plot(p, paste0(db, "_ORA_bar"), height = 5)
    } else {
      sig$direction <- factor(sig$direction, levels = names(pal))
      sig$plot_label <- paste0(wrap_text(sig$Description), " (", sig$ID, ")")
      sig$plot_label <- factor(sig$plot_label, levels = rev(unique(sig$plot_label)))
      sig$neglogFDR <- -log10(pmax(sig$p.adjust, .Machine$double.xmin))
      base <- ggplot(sig, aes(x = neglogFDR, y = plot_label, color = direction)) +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y", labeller = labeller(direction = dir_labels)) +
        scale_color_manual(values = pal, guide = "none") +
        labs(title = title,
          subtitle = "GroupX versus GroupY; top 8 significant terms per direction; database-specific tested-gene background",
          x = "-log10(BH-adjusted P)", y = NULL)
      height <- max(6.5, 2.3 + nrow(sig) * 0.60)
      bubble <- base + geom_point(aes(size = Count), alpha = 0.9) +
        scale_size_continuous(range = c(3, 9), name = "Overlapping genes")
      bar <- base + geom_col(aes(fill = direction), width = 0.7, show.legend = FALSE) +
        scale_fill_manual(values = pal)
      save_plot(bubble, paste0(db, "_ORA_bubble"), height = height)
      save_plot(bar, paste0(db, "_ORA_bar"), height = height)
    }
  }
  curve_records <- list()
  for (db in names(gsea_all)) {
    r <- gsea_all[[db]]
    sig <- r[!is.na(r$padj) & r$padj < 0.05, ]
    tops <- do.call(rbind, lapply(names(pal), function(d) head(sig[sig$direction == d, ], 8)))
    write_csv(tops, paste0(db, "_GSEA_plot_data.csv"))
    if (!nrow(tops)) {
      p <- empty_plot(paste(db, "GSEA"), "No gene sets pass BH-adjusted P < 0.05.\nNES > 0: GroupX; NES < 0: GroupY.")
      save_plot(p, paste0(db, "_GSEA_bar"), height = 5)
    } else {
      tops$direction <- factor(tops$direction, levels = names(pal))
      tops$plot_label <- paste0(wrap_text(tops$Description), " (", tops$ID, ")")
      tops$plot_label <- factor(tops$plot_label, levels = rev(unique(tops$plot_label)))
      p <- ggplot(tops, aes(x = NES, y = plot_label, fill = direction)) +
        geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
        geom_col(width = 0.7, show.legend = FALSE) + scale_fill_manual(values = pal) +
        facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
        labs(title = paste(db, "GSEA using all finite Wald statistics"),
          subtitle = "BH-adjusted P < 0.05; positive NES: GroupX higher; negative NES: GroupY higher",
          x = "Normalized enrichment score (NES)", y = NULL)
      save_plot(p, paste0(db, "_GSEA_bar"), height = max(6.5, 2.3 + nrow(tops) * 0.6))
    }
    for (d in names(pal)) {
      selected <- head(sig[sig$direction == d, ], 1)
      # 用细胞周期曲线展示增殖相关的富集结果。
      if (db == "KEGG" && d == "GroupX_higher" && any(sig$ID == "hsa04110"))
        selected <- sig[sig$ID == "hsa04110", ][1, ]
      if (!nrow(selected)) next
      id <- selected$ID[1]
      p <- fgsea::plotEnrichment(gsea_sets[[db]][[id]], ranks) +
        labs(title = paste0(db, ": ", wrap_text(selected$Description[1], 75)),
          subtitle = sprintf("%s | NES = %.3f | BH-adjusted P = %.3g | %s",
            id, selected$NES[1], selected$padj[1], d),
          x = "Rank among all mapped genes (left: GroupX; right: GroupY)",
          y = "Running enrichment score") +
        theme(plot.title = element_text(size = 14), plot.subtitle = element_text(size = 11))
      p$layers[[1]]$aes_params$colour <- unname(pal[d])
      p$layers[[3]]$aes_params$colour <- "grey55"
      p$layers[[4]]$aes_params$colour <- "grey55"
      save_plot(p, paste0(db, "_GSEA_curve_", d), width = 11, height = 6.2)
      curve_records[[paste(db, d, sep = "_")]] <- selected
    }
  }
  if (length(curve_records)) write_csv(do.call(rbind, curve_records), "GSEA_representative_curves.csv")

  counts <- do.call(rbind, c(lapply(names(ora_all), function(nm) {
    z <- ora_all[[nm]]
    data.frame(analysis = "ORA", database_direction = nm, tested_terms = nrow(z),
      significant_BH_005 = sum(z$p.adjust < 0.05, na.rm = TRUE))
  }), lapply(names(gsea_all), function(nm) {
    z <- gsea_all[[nm]]
    data.frame(analysis = "GSEA", database_direction = nm, tested_terms = nrow(z),
      significant_BH_005 = sum(z$padj < 0.05, na.rm = TRUE))
  })))
  write_csv(counts, "enrichment_result_counts.csv")
  gsea_direction_counts <- do.call(rbind, lapply(names(gsea_all), function(db) {
    z <- gsea_all[[db]]
    data.frame(database = db, GroupX_higher = sum(z$padj < 0.05 & z$NES > 0, na.rm = TRUE),
      GroupY_higher = sum(z$padj < 0.05 & z$NES < 0, na.rm = TRUE),
      missing_padj = sum(is.na(z$padj)))
  }))
  write_csv(gsea_direction_counts, "GSEA_significant_direction_counts.csv")
  parameters <- list(contrast = "GroupX versus GroupY", ORA = "Upper-tail hypergeometric",
    ORA_universe = "finite DESeq2 pvalue, mapped and annotated in the corresponding database",
    ORA_BH_family = "all eligible terms separately within each database and each direction, including zero-overlap terms",
    cutoff = "BH-adjusted P < 0.05", minGSSize = min_size, maxGSSize = max_size,
    GSEA_ranking = "all finite Wald statistics; max abs(stat) per ENTREZID; gene lexical tie-break; re-sort descending",
    GSEA_engine = "fgseaMultilevel", GSEA_version = as.character(packageVersion("fgsea")),
    GSEA_seed = seed, GSEA_eps = 0, GSEA_nPermSimple = 10000, GSEA_nproc = 1,
    database_metadata = b$metadata)
  saveRDS(list(ora = ora_all, gsea = gsea_all, ranks = ranks, rank_df = rank_df,
    mapping_summary = map_summary, background_summary = bg_summary,
    eligible_ORA_sets = eligible_sets, parameters = parameters, counts = counts),
    file.path(out, "中间对象", "enrichment_analysis.rds"), compress = "xz")
  capture.output(parameters, file = file.path(out, "日志", "enrichment_parameters.txt"))

  cat("富集结果：\n"); print(counts)
  list(ora = ora_all, gsea = gsea_all, ranks = ranks, rank_df = rank_df,
    mapping_summary = map_summary, background_summary = bg_summary, counts = counts)
}

# 9. 表达证据汇总与条码核对
infer_then_crosscheck <- function(differential, enrichment, out) {
  de <- differential$de
  # 汇总实际检出的增殖、肿瘤相关基因和肝功能基因。
  tumor_genes <- c("EZH2", "FOXM1", "GPC3", "BIRC5", "MKI67")
  liver_genes <- c("CYP2C19", "CYP3A4", "BCHE")
  evidence <- de[de$gene %in% c(tumor_genes, liver_genes), ]
  evidence$role <- ifelse(evidence$gene %in% tumor_genes, "Tumor-related/proliferation", "Liver function")
  write.csv(evidence, file.path(out, "结果", "biological_inference_genes.csv"),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  significant_higher <- function(genes, group) sum(de$gene %in% genes & de$direction == paste0(group, "_higher"))
  cellcycle <- enrichment$gsea$KEGG[enrichment$gsea$KEGG$ID == "hsa04110", ]
  retinol <- enrichment$gsea$KEGG[enrichment$gsea$KEGG$ID == "hsa00830", ]
  pathway_support <- function(group) {
    sign <- if (group == "GroupX") 1 else -1
    nrow(cellcycle) == 1 && nrow(retinol) == 1 &&
      is.finite(cellcycle$padj) && cellcycle$padj < 0.05 && cellcycle$NES * sign > 0 &&
      is.finite(retinol$padj) && retinol$padj < 0.05 && retinol$NES * sign < 0
  }
  x_support <- significant_higher(tumor_genes, "GroupX") >= 3 &&
    significant_higher(liver_genes, "GroupY") >= 2 && pathway_support("GroupX")
  y_support <- significant_higher(tumor_genes, "GroupY") >= 3 &&
    significant_higher(liver_genes, "GroupX") >= 2 && pathway_support("GroupY")
  tumor_group <- if (x_support && !y_support) "GroupX" else if (y_support && !x_support) "GroupY" else "Unresolved"
  # 这里的组别推断用于本批样本；表达上调不能作为突变证据。
  inference <- c(paste("Expression/pathway inference before sample-type decoding:", tumor_group, "is tumor-like."),
    "Criteria: >=3 significant tumor-related genes, >=2 opposite liver-function genes, concordant Cell cycle/Retinol GSEA.",
    "Gene function sources: EZH2 PMID 17596871; FOXM1 PMID 15082532; GPC3 PMID 16024626; BIRC5 PMID 9256286.",
    "Bulk immune pathway enrichment cannot establish cell fractions, exhaustion, treatment response or causality.",
    "Gene expression changes do not demonstrate mutations; adjacent tissue is not necessarily healthy liver.")
  writeLines(inference, file.path(out, "日志", "inference_before_barcode.txt"))
  cat(paste(inference, collapse = "\n"), "\n")
  # 最后按 NCI GDC 的样本类型定义核对条码。
  # https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
  x <- data.frame(sample = rownames(differential$meta), differential$meta, check.names = FALSE)
  x$sample_type_code <- substr(vapply(strsplit(x$sample, "-", fixed = TRUE), `[`, character(1), 4), 1, 2)
  x$GDC_definition <- ifelse(x$sample_type_code == "01", "Primary Solid Tumor",
    ifelse(x$sample_type_code == "11", "Solid Tissue Normal", "Other"))
  x$expression_inference <- if (tumor_group == "Unresolved") "Unresolved" else
    ifelse(as.character(x$Group) == tumor_group, "Tumor", "Adjacent non-tumor")
  x$concordant <- if (tumor_group == "Unresolved") NA else
    (x$expression_inference == "Tumor" & x$sample_type_code == "01") |
    (x$expression_inference == "Adjacent non-tumor" & x$sample_type_code == "11")
  write.csv(x, file.path(out, "结果", "independent_barcode_crosscheck.csv"),
    row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
  print(x[, c("sample", "Group", "patient_label", "sample_type_code", "GDC_definition", "concordant")])
  if (tumor_group == "Unresolved" || !all(x$concordant))
    warning("组别证据不足或与条码不一致，请检查本次结果，不强行给出预设结论。")
  invisible(x)
}

main <- function(cfg = CONFIG) {
  options(stringsAsFactors = FALSE, timeout = 600, warn = 1)
  # 统一排序环境以复现同值统计量的字典序；在不支持该locale的平台继续使用现有UTF-8环境。
  try(suppressWarnings(Sys.setlocale("LC_COLLATE", "C")), silent = TRUE)
  try(suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8")), silent = TRUE)
  cfg$project_dir <- normalizePath(cfg$project_dir, winslash = "/", mustWork = TRUE)
  out <- resolve_path(cfg$output_dir, cfg$project_dir)
  for (directory in c("", "结果", "图表", "日志", "数据库", "中间对象", "原始输入"))
    dir.create(file.path(out, directory), recursive = TRUE, showWarnings = FALSE)
  out <- normalizePath(out, winslash = "/", mustWork = TRUE)
  log_connection <- file(file.path(out, "日志", "complete_analysis.log"), "wt", encoding = "UTF-8")
  initial_sinks <- sink.number()
  sink(log_connection, split = TRUE); sink(log_connection, type = "message")
  on.exit({
    sink(type = "message")
    while (sink.number() > initial_sinks) sink()
    close(log_connection)
  }, add = TRUE)
  cat("完整单文件分析开始：", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), "\n")
  print(cfg); cat(R.version.string, "\n")
  tryCatch({
    snapshot <- resolve_path(cfg$database_snapshot, cfg$project_dir)
    ensure_packages(cfg, out, need_database_packages = !nzchar(cfg$database_snapshot) || !file.exists(snapshot))
    # 核对报告采用的分析参数。
    stopifnot(cfg$seed == 23360036L, cfg$integer_tolerance == 1e-6,
      cfg$min_count == 10L, cfg$min_samples == 4L, cfg$padj_cutoff == 0.05,
      cfg$lfc_cutoff == 1, cfg$min_gene_set_size == 10L, cfg$max_gene_set_size == 500L)
    differential <- run_differential(cfg, out)
    bundle <- load_or_build_database(cfg, out)
    enrichment <- run_enrichment(differential, bundle, cfg, out)
    crosscheck <- infer_then_crosscheck(differential, enrichment, out)
    # 检查表格、排名和绘图文件是否完整。
    png_files <- list.files(file.path(out, "图表"), pattern = "[.]png$", full.names = TRUE)
    pdf_files <- list.files(file.path(out, "图表"), pattern = "[.]pdf$", full.names = TRUE)
    stopifnot(nrow(differential$de) == differential$qc$modeled_genes,
      !anyDuplicated(enrichment$rank_df$ENTREZID), all(is.finite(enrichment$ranks)),
      all(diff(enrichment$ranks) <= 0), length(png_files) >= 10,
      all(file.info(c(png_files, pdf_files))$size > 1000),
      setequal(tools::file_path_sans_ext(basename(png_files)),
        tools::file_path_sans_ext(basename(pdf_files))))
    summary <- list(DESeq2 = differential$qc, enrichment = enrichment$counts,
      ranked_genes = length(enrichment$ranks), barcode_concordant = sum(crosscheck$concordant, na.rm = TRUE))
    jsonlite::write_json(summary, file.path(out, "结果", "complete_analysis_summary.json"),
      auto_unbox = TRUE, pretty = TRUE, na = "null")
    capture.output(sessionInfo(), file = file.path(out, "日志", "sessionInfo.txt"))
    cat("分析结束：", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), "\n输出目录：", out, "\n")
    print(summary)
    invisible(summary)
  }, error = function(e) {
    cat("分析中止：", conditionMessage(e), "\n")
    capture.output(sessionInfo(), file = file.path(out, "日志", "sessionInfo_error.txt"))
    stop(e)
  })
}

main()
