#!/usr/bin/env Rscript

# Create PLINK score file (3 columns) from matched SNPs and meta-GWAS summary stats.
# USAGE:
#   Rscript scripts/30_scores/get_score.R \
#     --matched data/interim/matched/matched.tsv \
#     --sumstats data/raw/meta-GWAS-Wightman-Bellenguez-raw.tbl \
#     --out data/interim/scores/global.score \
#     --exclude-apoe true

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# -------------------- Command-line interface --------------------
parser <- ArgumentParser()
parser$add_argument("--matched", required=TRUE,
                    help="TSV with RSID + ADNI_VARIANT_ID (+ optional ADNI_REF/ADNI_ALT, CHR, POS)")
parser$add_argument("--sumstats", required=TRUE,
                    help="Union meta-GWAS summary stats TSV (Bellenguez + Wightman)")
parser$add_argument("--out", required=TRUE,
                    help="Output PLINK score path (3-col; header): SNP A1 BETA")
parser$add_argument("--exclude-apoe", default="true",
                    help="Exclude rs429358/rs7412 from scoring (true/false)")
args <- parser$parse_args()
exclude_apoe <- tolower(args$`exclude-apoe`) %in% c("true","1","yes","y")

cat("Creating PLINK score file\n")
cat(" matched:  ", args$matched,  "\n")
cat(" sumstats: ", args$sumstats, "\n")
cat(" out:      ", args$out,      "\n")
cat(" exclude APOE: ", exclude_apoe, "\n\n")

# -------------------- Read inputs --------------------
matched <- fread(args$matched)
sumstats <- fread(args$sumstats)

# Expect typical columns
if (!all(c("RSID","ADNI_VARIANT_ID") %in% names(matched))) {
  stop("matched file must contain at least: RSID, ADNI_VARIANT_ID")
}
needed_ss <- c("MarkerName","Allele1","Allele2","Effect","StdErr","P-value")
if (!all(needed_ss %in% names(sumstats))) {
  stop("sumstats must contain columns: ", paste(needed_ss, collapse=", "))
}

# Uppercase/refine alleles
sumstats <- sumstats %>% mutate(
  Allele1 = toupper(Allele1),
  Allele2 = toupper(Allele2)
)

# If ADNI_REF/ADNI_ALT not present, parse from ADNI_VARIANT_ID (chr:pos:ref:alt)
if (!all(c("ADNI_REF","ADNI_ALT") %in% names(matched))) {
  parts <- str_split(matched$ADNI_VARIANT_ID, ":", simplify = TRUE)
  if (ncol(parts) < 4) stop("ADNI_VARIANT_ID must be 'chr:pos:ref:alt'")
  matched$ADNI_REF <- toupper(parts[,3])
  matched$ADNI_ALT <- toupper(parts[,4])
} else {
  matched$ADNI_REF <- toupper(matched$ADNI_REF)
  matched$ADNI_ALT <- toupper(matched$ADNI_ALT)
}

# Optional: drop APOE rsIDs for scoring
if (exclude_apoe) {
  matched <- matched %>% filter(!RSID %in% c("rs429358","rs7412"))
}

# -------------------- Merge & allele handling --------------------
merged <- matched %>%
  left_join(sumstats, by = c("RSID" = "MarkerName")) %>%
  filter(!is.na(Effect))

# base complement
comp <- function(a) dplyr::case_when(
  a == "A" ~ "T", a == "T" ~ "A",
  a == "G" ~ "C", a == "C" ~ "G",
  TRUE ~ a
)

# Choose risk allele so that BETA is positive; align to cohort alleles
pick_risk <- function(adni_ref, adni_alt, a1, a2, eff) {
  if (is.na(eff)) return(list(a=adni_alt, b=NA_real_))
  risk <- if (eff >= 0) a1 else a2
  beta <- abs(eff)

  # direct match
  if (risk == adni_ref) return(list(a=adni_ref, b=beta))
  if (risk == adni_alt) return(list(a=adni_alt, b=beta))

  # strand-flip match
  if (risk == comp(adni_ref)) return(list(a=adni_ref, b=beta))
  if (risk == comp(adni_alt)) return(list(a=adni_alt, b=beta))

  # fallback (ALT)
  list(a=adni_alt, b=beta)
}

res <- merged %>%
  rowwise() %>%
  mutate(tmp = list(pick_risk(ADNI_REF, ADNI_ALT, Allele1, Allele2, Effect)),
         RISK_ALLELE = tmp$a,
         FINAL_EFFECT = tmp$b) %>%
  ungroup()

# -------------------- Outputs --------------------
# 3-col PLINK score
score <- res %>%
  transmute(SNP = ADNI_VARIANT_ID, A1 = RISK_ALLELE, BETA = FINAL_EFFECT) %>%
  distinct(SNP, .keep_all = TRUE)

# detailed mapping (optional but useful)
detailed <- res %>%
  select(RSID, ADNI_VARIANT_ID, ADNI_REF, ADNI_ALT,
         Allele1, Allele2, Effect, RISK_ALLELE, FINAL_EFFECT, `P-value`)

# missing rsIDs (in matched but absent in sumstats)
missing <- matched %>%
  anti_join(sumstats, by = c("RSID" = "MarkerName")) %>%
  select(RSID, ADNI_VARIANT_ID, ADNI_REF, ADNI_ALT)

# write files
fwrite(score, args$out, sep = "\t", quote = FALSE)

prefix <- sub("\\.score$","", args$out)
fwrite(detailed, paste0(prefix, "_detailed.txt"), sep = "\t", quote = FALSE)
fwrite(missing,  paste0(prefix, "_missing_in_sumstats.txt"), sep = "\t", quote = FALSE)

# summary text
summary_text <- paste0(
  "=== PLINK Score File Creation Summary ===\n",
  "Input:\n",
  " Matched SNPs: ", nrow(matched), "\n",
  " Summary statistics: ", nrow(sumstats), "\n\n",
  "Output:\n",
  " SNPs with effect sizes: ", nrow(score), "\n",
  " SNPs missing from summary stats: ", nrow(missing), "\n\n",
  "Files created:\n",
  " PLINK score file: ", args$out, "\n",
  " Detailed mapping: ", prefix, "_detailed.txt\n",
  " Missing SNPs: ", prefix, "_missing_in_sumstats.txt\n"
)
writeLines(summary_text, paste0(prefix, "_summary.txt"))
cat(summary_text)
