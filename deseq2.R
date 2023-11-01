suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-s", "--samples"), type="character", default="samples.tsv", help="Samples file in TSV format [default \"%default\"]", dest="samples_file"),
  make_option(c("-e", "--experimental_conditions"), type="character", default="mut", help="Comma-separated experimental conditions [default \"%default\"]", dest="exp_cond"),
  make_option(c("-c", "--control_conditions"), type="character", default="sib", help="Comma-separated control conditions [default \"%default\"]", dest="con_cond"),
  make_option(c("-a", "--annotation"), type="character", default="annotation.tsv", help="Annotation file in TSV format [default \"%default\"]", dest="annotation_file"),
  make_option(c("-d", "--star_dir"), type="character", default="star", help="Top level STAR directory [default \"%default\"]", dest="star_dir"),
  make_option(c("-o", "--output_dir"), type="character", default="deseq2", help="Output directory [default \"%default\"]", dest="output_dir"),
  make_option(c("-p", "--pca_genes"), type="integer", default="500", help="Number of genes used in PCA [default \"%default\"]", dest="pca_genes"),
  make_option(c("--pc_threshold"), type="double", default="1", help="Threshold for percentage variance explained by PCs [default \"%default\"]", dest="pc_threshold"),
  make_option(c("--gene_threshold"), type="double", default="0.1", help="Threshold for percentage variance explained by genes [default \"%default\"]", dest="gene_threshold")
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(tidyverse))
options(readr.show_progress=FALSE)
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(DESeq2))

custom_palette <- function(palette_size = 8) {
  colour_blind_palette <- c(
    "vermillion" = rgb(0.8, 0.4, 0),
    "blue_green" = rgb(0, 0.6, 0.5),
    "blue" = rgb(0,0.45,0.7),
    "yellow" = rgb(0.95, 0.9, 0.25),
    "sky_blue" = rgb(0.35, 0.7, 0.9),
    "purple" = rgb(0.8, 0.6, 0.7),
    "black" = rgb(0, 0, 0),
    "orange" = rgb(0.9, 0.6, 0),
    "grey60" = "#999999",
    "grey20" = "#333333"
  )
  if (palette_size <= 10) {
    palette <- colour_blind_palette[seq_len(palette_size)]
    palette <- unname(palette)
  } else {
    palette <- scales::hue_pal()(palette_size)
  }

  return(palette)
}

# Get samples
num_columns <- count_fields(opt$samples_file, tokenizer_tsv(), n_max=1)
if (num_columns == 2) {
  samples <- read_tsv(
    opt$samples_file,
    col_names=c("sample", "condition"),
    col_types=cols(
      sample=col_character(),
      condition=col_factor()
    )
  )
  stop_for_problems(samples)
} else if (num_columns == 3) {
  samples <- read_tsv(
    opt$samples_file,
    col_names=c("sample", "condition", "group"),
    col_types=cols(
      sample=col_character(),
      condition=col_factor(),
      group=col_factor()
    )
  )
  stop_for_problems(samples)
} else {
  stop("Samples file must have two or three columns")
}

# Get comparison
opt$exp_cond <- str_split(opt$exp_cond, ",")[[1]]
opt$con_cond <- str_split(opt$con_cond, ",")[[1]]
#samples$condition <- fct_relevel(samples$condition, sort)
#if (num_columns == 3) {
#  samples$group <- fct_relevel(samples$group, sort)
#}
samples_deseq2 <- samples
samples_deseq2$condition <- fct_collapse(
  samples_deseq2$condition,
  exp=opt$exp_cond,
  con=opt$con_cond
)
samples_deseq2$condition <- fct_relabel(
  samples_deseq2$condition,
  ~ str_replace_all(.x, "-", "_")
)
if (num_columns == 3) {
  samples_deseq2$group <- fct_relabel(
    samples_deseq2$group,
    ~ str_replace_all(.x, "-", "_")
  )
}
samples_deseq2$condition <- fct_relevel(samples_deseq2$condition, "exp", "con")

# Get annotation
annotation <- read_tsv(
  opt$annotation_file,
  col_names=c(
    "Gene", "Chr", "Start", "End", "Strand",
    "Biotype", "Name", "Description"
  ),
  col_types=cols(
    Gene=col_character(),
    Chr=col_factor(),
    Start=col_integer(),
    End=col_integer(),
    Strand=col_integer(),
    Biotype=col_factor(),
    Name=col_character(),
    Description=col_character()
  )
)

# Get counts
counts <-tibble(
  gene=character(),
  sample=factor(samples$sample, levels=unique(samples$sample))[0],
  stranded2=integer(),
  .rows=0
)
for (sample in samples$sample) {
  tmp_counts <- read_tsv(
    str_c(opt$star_dir, "/", sample, "/ReadsPerGene.out.tab"),
    skip=4,
    col_names=c("gene", "unstranded", "stranded1", "stranded2"),
    col_types=cols(
      gene=col_character(),
      unstranded=col_skip(),
      stranded1=col_skip(),
      stranded2=col_integer()
    )
  )
  stop_for_problems(tmp_counts)
  tmp_counts <- add_column(
    tmp_counts,
    sample=factor(sample, levels=levels(counts$sample)),
    .before="stranded2"
  )
  counts <- bind_rows(counts, tmp_counts)
}
rm(tmp_counts)

# Aggregate counts
counts <- pivot_wider(counts, names_from=sample, values_from=stranded2)

# DESeq2
dir.create(opt$output_dir, showWarnings=FALSE, recursive=TRUE)
if (num_columns == 2) {
  design <- formula(~ condition)
} else if (num_columns == 3) {
  design <- formula(~ group + condition + group:condition)
}
dds <- DESeqDataSetFromMatrix(
  column_to_rownames(counts, var="gene"),
  column_to_rownames(samples_deseq2, var="sample"),
  design=design
)
rm(counts)
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition", "exp", "con"), alpha=0.05, tidy=TRUE)

# Output
write_tsv(
  enframe(dds$sizeFactor),
  str_c(opt$output_dir, "/", "size-factors.tsv"),
  col_names=FALSE
)
res <- inner_join(
  select(res, Gene=row, pval=pvalue, adjp=padj, log2fc=log2FoldChange),
  annotation,
  by="Gene"
)
res <- inner_join(
  res,
  rownames_to_column(
    rename_all(
      as.data.frame(counts(dds)),
      paste0,
      " count"
    ),
    var="Gene"
  ),
  by="Gene"
)
res <- inner_join(
  res,
  rownames_to_column(
    rename_all(
      as.data.frame(counts(dds, normalized=TRUE)),
      paste0,
      " normalised count"
    ),
    var="Gene"
  ),
  by="Gene"
)
res <- arrange(res, adjp, pval, Gene)
write_tsv(res, str_c(opt$output_dir, "/", "all.tsv"))
write_tsv(filter(res, adjp < 0.05), str_c(opt$output_dir, "/", "sig.tsv"))

# QC
rld <- rlog(dds, blind=TRUE)
colData(rld)$condition <- samples$condition
pdf(str_c(opt$output_dir, "/", "qc.pdf"))
if (num_columns == 2) {
  plotPCA(rld, intgroup="condition")
} else if (num_columns == 3) {
  plotPCA(rld, intgroup=c("condition", "group"))
}
plotMA(dds)
coef <- "condition_exp_vs_con"
if (!(coef %in% resultsNames(dds))) {
  coef <- "condition_con_vs_exp"
}
plotMA(lfcShrink(dds, coef=coef, quiet=TRUE))
plotDispEsts(dds)
dev.off()

# PCA
rv <- rowVars(assay(rld))
select <- order(rv, decreasing=TRUE)[seq_len(min(opt$pca_genes, length(rv)))]
pca <- prcomp(t(assay(rld)[select,]))
propVarPC <- pca$sdev^2 / sum( pca$sdev^2 )
aload <- abs(pca$rotation)
propVarGene <- sweep(aload, 2, colSums(aload), "/")

# Output PC coordinates above threshold
lastSigPC <- sum(propVarPC * 100 >= opt$pc_threshold)
x <- pca$x[,1:lastSigPC]
x <- rownames_to_column(as.data.frame(x), var="Sample")
write_tsv(x, str_c(opt$output_dir, "/", "PCs.tsv"))

# Output genes contributing most to each PC
for (i in seq.int(lastSigPC)) {
  propPC <- as.data.frame(propVarGene[,i] * 100)
  colnames(propPC) <- "% variance explained"
  resPC <- inner_join(
    res,
    rownames_to_column(propPC, var="Gene"),
    by="Gene"
  )
  resPC <- arrange(resPC, desc(`% variance explained`), Gene)
  resPC <- filter(resPC, `% variance explained` >= opt$gene_threshold)
  write_tsv(resPC, str_c(opt$output_dir, "/", "PC-", i, ".tsv"))
}

# Plot PCA
colour_palette <- custom_palette(nlevels(samples$condition))
names(colour_palette) <- levels(samples$condition)
pdf(str_c(opt$output_dir, "/", "pca.pdf"))
intgroup <- "condition"
if (num_columns == 3) {
  intgroup <- c("condition", "group")
}
intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop=FALSE])
group <- factor(apply( intgroup.df, 1, paste, collapse=" : "))
if (num_columns == 3) {
  colour_palette <- custom_palette(nlevels(group))
}
short_sample_names <- colnames(rld)
for (i in seq.int(lastSigPC - 1)) {
  first <- i
  second <- i + 1
  d <- data.frame(first=pca$x[,first], second=pca$x[,second], group=group, intgroup.df, name=colnames(rld))
  p <- ggplot(data=d, aes_string(x="first", y="second", fill="group")) +
    coord_cartesian(clip="off") +
    geom_point(size=3, shape=21, colour="black") +
    scale_fill_manual(values = colour_palette) +
    geom_text_repel(aes(label=short_sample_names), segment.alpha=0.3, hjust=0, vjust=0, size=2, show.legend=FALSE) +
    xlab(paste0("PC", first, ": ", round(propVarPC[first] * 100, 1), "% variance")) +
    ylab(paste0("PC", second, ": ", round(propVarPC[second] * 100, 1), "% variance")) +
    theme_minimal(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90", colour = "grey90"), legend.position="bottom")
  print(p)
}
dev.off()

# Plot counts
colour_palette <- custom_palette(nlevels(samples$condition))
names(colour_palette) <- levels(samples$condition)
pdf(str_c(opt$output_dir, "/", "counts.pdf"))
for (i in seq.int(nrow(filter(res, adjp < 0.05)))) {
  countData <- res[i,grepl("normalised", names(res))]
  names(countData) <- gsub(" normalised.*$", "", names(countData))
  counts <- pivot_longer(countData, everything(), names_to="sample", values_to="count")
  counts <- inner_join(counts, samples, by="sample")
  counts$sample <- fct_relevel(counts$sample, samples$sample)
  title <- sprintf("%s / %s\n%s:%d-%d\np = %.3g", res$Gene[i], res$Name[i], res$Chr[i], res$Start[i], res$End[i], res$adjp[i])
  p <- ggplot(data = counts, aes_(x = as.name("sample"), y = as.name("count"))) +
    geom_point(aes_(fill = as.name("condition")), size = 3, shape = 21, colour = "black") +
    scale_fill_manual(values = colour_palette) +
    labs(title = title, x = NULL, y = "Normalised Counts") +
    theme_minimal(base_size = 12) +
    theme(strip.background = element_rect(fill = "grey90", colour = "grey90"), axis.text.x = element_text(angle = 90))
  print(p)
}
dev.off()
