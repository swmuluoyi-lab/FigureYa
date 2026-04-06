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
tcgaBLCA <- read.table("exp.mtx.txt", sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
pheno <- read.table("easyinput_exp_phenotype.txt", sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, header = TRUE)
tfs <- read.table("easyinput_regulon.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

message("[2/8] Intersecting regulons with expression matrix genes...")
regulatoryElements <- intersect(tfs$regulon, rownames(tcgaBLCA))
missing_regulons <- setdiff(tfs$regulon, regulatoryElements)

message(sprintf("Input regulons: %d | Matched: %d | Missing: %d", nrow(tfs), length(regulatoryElements), length(missing_regulons)))
if (length(regulatoryElements) < 20) {
  stop("Too few matched regulons. Please inspect easyinput_regulon.txt and exp.mtx.txt")
}

message("[3/8] Constructing TNI object...")
rtni_tcgaBLCA <- tni.constructor(
  expData = as.matrix(log2(tcgaBLCA + 1)),
  regulatoryElements = regulatoryElements
)

message("[4/8] Running permutation/bootstrap (long step)...")
options(cluster = snow::makeCluster(spec = 4, "SOCK"))
rtni_tcgaBLCA <- tni.permutation(rtni_tcgaBLCA, pValueCutoff = 1e-5, nPermutations = 1000)
rtni_tcgaBLCA <- tni.bootstrap(rtni_tcgaBLCA, nBootstraps = 1000)
stopCluster(getOption("cluster"))

message("[5/8] Running DPI filter...")
rtni_tcgaBLCA <- tni.dpi.filter(rtni_tcgaBLCA, eps = 0, sizeThreshold = TRUE, minRegulonSize = 15)
save(rtni_tcgaBLCA, file = "rtni_tcgaBLCA_lipid50.RData")

message("[6/8] Running two-sided GSEA for regulon activity...")
rtnigsea_tcgaBLCA <- tni.gsea2(rtni_tcgaBLCA, regulatoryElements = regulatoryElements)
MIBC_regact <- tni.get(rtnigsea_tcgaBLCA, what = "regulonActivity")
save(MIBC_regact, file = "MIBC_regact_lipid50.RData")

message("[7/8] Preparing matrix for heatmap...")
grp <- pheno[colnames(MIBC_regact), 1]
regact <- MIBC_regact[, order(grp)]
regact <- standarize.fun(regact, halfwidth = 1.5)

message("[8/8] Drawing heatmap and exporting summary...")
pdf("regulon_heatmap_lipid50.pdf", width = 9, height = 7)
heatmap.2(
  as.matrix(regact),
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
  sprintf("Samples in expression matrix: %d", ncol(tcgaBLCA)),
  sprintf("Genes in expression matrix: %d", nrow(tcgaBLCA)),
  sprintf("Input regulons: %d", nrow(tfs)),
  sprintf("Matched regulons: %d", length(regulatoryElements)),
  sprintf("Missing regulons: %d", length(missing_regulons)),
  ifelse(length(missing_regulons) > 0, paste("Missing list:", paste(missing_regulons, collapse = ", ")), "Missing list: None")
)
writeLines(summary_lines, con = "regulon_analysis_lipid50_summary.txt")

message("Done. Outputs:")
message(" - rtni_tcgaBLCA_lipid50.RData")
message(" - MIBC_regact_lipid50.RData")
message(" - regulon_heatmap_lipid50.pdf")
message(" - regulon_analysis_lipid50_summary.txt")
