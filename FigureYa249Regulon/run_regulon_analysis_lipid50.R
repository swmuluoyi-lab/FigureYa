source("install_dependencies.R")

library(RTN)
library(snow)
library(ComplexHeatmap)
library(ClassDiscovery)
library(RColorBrewer)
library(gplots)

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

standarize.fun <- function(indata = NULL, halfwidth = NULL, centerFlag = TRUE, scaleFlag = TRUE) {
  outdata <- t(scale(t(indata), center = centerFlag, scale = scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata > halfwidth] <- halfwidth
    outdata[outdata < (-halfwidth)] <- -halfwidth
  }
  return(outdata)
}

message("[1/8] Loading input files...")
hepg2_oa_expr <- read.table("exp.mtx.txt", sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
hepg2_group <- read.table("easyinput_exp_phenotype.txt", sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
lipid_tfs <- read.table("easyinput_regulon.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

message("[2/8] Intersecting regulons with expression matrix genes...")
regulatory_elements <- intersect(lipid_tfs$regulon, rownames(hepg2_oa_expr))
missing_regulons <- setdiff(lipid_tfs$regulon, regulatory_elements)

message(sprintf("Input regulons: %d | Matched: %d | Missing: %d", nrow(lipid_tfs), length(regulatory_elements), length(missing_regulons)))
if (length(regulatory_elements) < 20) {
  stop("Too few matched regulons. Please inspect easyinput_regulon.txt and exp.mtx.txt")
}

message("[3/8] Constructing TNI object...")
rtni_hepg2_oa <- tni.constructor(
  expData = as.matrix(log2(hepg2_oa_expr + 1)),
  regulatoryElements = regulatory_elements
)

message("[4/8] Running permutation/bootstrap (long step)...")
options(cluster = snow::makeCluster(spec = 4, "SOCK"))
rtni_hepg2_oa <- tni.permutation(rtni_hepg2_oa, pValueCutoff = 1e-5, nPermutations = 1000)
rtni_hepg2_oa <- tni.bootstrap(rtni_hepg2_oa, nBootstraps = 1000)
stopCluster(getOption("cluster"))

message("[5/8] Running DPI filter...")
rtni_hepg2_oa <- tni.dpi.filter(rtni_hepg2_oa, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15)
save(rtni_hepg2_oa, file = "rtni_hepg2_oa_lipid50.RData")

message("[6/8] Running two-sided GSEA for regulon activity...")
rtnigsea_hepg2_oa <- tni.gsea2(rtni_hepg2_oa, regulatoryElements = regulatory_elements)
hepg2_regulon_activity <- tni.get(rtnigsea_hepg2_oa, what = "regulonActivity")
save(hepg2_regulon_activity, file = "hepg2_regact_lipid50.RData")

message("[7/8] Preparing matrix for heatmap...")
group_label <- hepg2_group[colnames(hepg2_regulon_activity), 1]
regulon_activity_zscore <- hepg2_regulon_activity[, order(group_label)]
regulon_activity_zscore <- standarize.fun(regulon_activity_zscore, halfwidth = 1.5)

message("[8/8] Drawing heatmap and exporting summary...")
pdf("regulon_heatmap_lipid50.pdf", width = 9, height = 7)
heatmap.2(
  as.matrix(regulon_activity_zscore),
  scale = "none",
  Colv = FALSE,
  trace = "none",
  density.info = "none",
  col = colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(100),
  margins = c(9, 9),
  key.title = "regulon\nactivity"
)
dev.off()

summary_lines <- c(
  sprintf("Date: %s", Sys.Date()),
  sprintf("Dataset: HepG2 OA model (Model) vs OA+Sesamin treatment (Treat)"),
  sprintf("Samples in expression matrix: %d", ncol(hepg2_oa_expr)),
  sprintf("Genes in expression matrix: %d", nrow(hepg2_oa_expr)),
  sprintf("Input regulons: %d", nrow(lipid_tfs)),
  sprintf("Matched regulons: %d", length(regulatory_elements)),
  sprintf("Missing regulons: %d", length(missing_regulons)),
  ifelse(length(missing_regulons) > 0, paste("Missing list:", paste(missing_regulons, collapse = ", ")), "Missing list: None")
)
writeLines(summary_lines, con = "regulon_analysis_lipid50_summary.txt")

message("Done. Outputs:")
message(" - rtni_hepg2_oa_lipid50.RData")
message(" - hepg2_regact_lipid50.RData")
message(" - regulon_heatmap_lipid50.pdf")
message(" - regulon_analysis_lipid50_summary.txt")
