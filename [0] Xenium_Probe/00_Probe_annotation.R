
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# -------- helper: normalize Potri ID to Potri. --------
norm_potri <- function(x) {
  x <- as.character(x)
  x <- gsub("^Potri[-_]", "Potri.", x)  # Potri- / Potri_ -> Potri.
  x
}

# =========================
# 0) Paths
# =========================
# Path to the original 300-gene file
GENE_TSV  <- "/home/woodydrylab/FileShare/20260121_Xenium/output-XETG00360__0079031__BIO1_TISSUE_1_and_2__20260115__224443/grid05um_out/genes.tsv"
# Path to the functional annotation reference table
FUNC_PATH <- "/home/woodydrylab/FileShare/20260121_Xenium/Ptrichocarpa_v4.1_complete table_v1.1.csv"
# Output path
OUTPUT_FILE <- "/home/woodydrylab/FileShare/20260121_Xenium/Xenium_300genes_Annotation.tsv"

# =========================
# 1) Load functional annotation table
# =========================
cat("Loading functional annotation...\n")
func <- fread(FUNC_PATH, sep = ",")
stopifnot("Gene.ID" %in% colnames(func))
func[, Gene.ID := as.character(Gene.ID)]

# Keep key annotation columns
keep_func_cols <- intersect(
  c("Gene.ID","TF","TF.family_v4.1","TF.family_v3.0","Common.name","SCW",
    "Best.hit.Arabidopsis","Annotated.function",
    "Leaves_Ave","Xylem_Ave","Phloem_Ave","YoungShoot_Ave",
    "Xylem.specific.gene", "Xylem.high.expression.gene"),
  colnames(func)
)
func <- func[, ..keep_func_cols]
setkey(func, Gene.ID)
func <- unique(func, by = "Gene.ID")

# =========================
# 2) Load Xenium Genes & Merge
# =========================
cat("Processing Xenium genes.tsv...\n")
# Note: Xenium genes.tsv often has no header, and the first column is usually the gene name (ID)
# Here we assume a simple file structure and name the column "gene" after reading
genes_dt <- fread(GENE_TSV, header = TRUE, col.names = "gene")

# Normalize the Potri ID and append the .v4.1 suffix to match Gene.ID in the functional table
genes_dt[, gene_norm := norm_potri(gene)]
genes_dt[, Gene.ID := paste0(gene_norm, ".v4.1")]

# Perform the merge (left join)
annotated_300 <- merge(genes_dt, func, by = "Gene.ID", all.x = TRUE)

# Reorder columns so identifiers and core annotations appear first
front_cols <- intersect(
  c("gene", "Gene.ID", "Best.hit.Arabidopsis", "Annotated.function", "Common.name", "TF", "SCW"),
  colnames(annotated_300)
)
other_cols <- setdiff(names(annotated_300), c(front_cols, "gene_norm"))
setcolorder(annotated_300, c(front_cols, other_cols))

# Remove the temporary intermediate column used during conversion
annotated_300[, gene_norm := NULL]

# =========================
# 3) Output
# =========================
fwrite(annotated_300, OUTPUT_FILE, sep = "\t")
cat("\nDone!\n")
cat("Total genes processed:", nrow(annotated_300), "\n")
cat("Mapped with function:", sum(!is.na(annotated_300$Annotated.function)), "\n")
cat("Output saved to:", OUTPUT_FILE, "\n")

