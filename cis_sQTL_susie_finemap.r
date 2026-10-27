library(data.table)
library(vcfR)
library(susieR)
library(Rfast)

root_dir_prior <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus/reported_qtl_prepare"
root_dir_bed   <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/junction_sc/junction_file/all_ds/leafcutter/total_bed"
root_dir_summ  <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus"
root_dir_vcf   <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/merged/annotation"
output_dir     <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus/susie_results"
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

prior_dirs <- list.dirs(root_dir_prior, recursive=FALSE)
cell_types <- basename(prior_dirs)
cell_types <- sub('\\.enrichment\\.prior$', '', cell_types)

script_dir <- "/public/users/wuqingqin/wuqingqin/STAR_CELL/eqtl/iqt/torus/scripts_per_celltype"
dir.create(script_dir, recursive=TRUE, showWarnings=FALSE)

for (cell in cell_types) {
  script_file <- file.path(script_dir, paste0("susie_", cell, ".R"))
  
  script_content <- paste0(
    "library(data.table)
library(vcfR)
library(susieR)
library(Rfast)

dir_prior <- '", file.path(root_dir_prior, paste0(cell,'.enrichment.prior')), "'
file_bed <- '", file.path(root_dir_bed, paste0(cell,'.merged.bed.gz')), "'
file_summ <- '", file.path(root_dir_summ, paste0(cell,'.merged.summ.gz')), "'
file_vcf <- '", file.path(root_dir_vcf, paste0(cell,'.merged_filter_annotated_final.vcf.gz')), "'
output_dir <- '", output_dir, "'

bed_all <- fread(file_bed)
summ_all <- fread(file_summ)
vcf <- read.vcfR(file_vcf)
gt <- extract.gt(vcf, element='GT')
gt_num <- apply(gt, 2, function(x) {
  sapply(x, function(g) { alleles <- as.numeric(unlist(strsplit(g, '[/|]'))); sum(alleles) })
})

prior_files <- list.files(dir_prior, pattern='prior$', full.names=TRUE)

fitted_rss_batch <- list()
batch_size <- 500
batch_counter <- 1
out_all <- list()
out_max <- list()

for (pf in prior_files) {
  locus <- sub('.prior$', '', basename(pf))
  message('Processing locus: ', locus)

  prior_df <- fread(pf)
  colnames(prior_df) <- c('snp','prior')
  prior_vec <- prior_df$prior
  names(prior_vec) <- prior_df$snp

  summ_sub <- summ_all[gene == locus]
  if (nrow(summ_sub)==0) next
  summ_sub$se <- summ_sub$beta / summ_sub$`t-stat`
  bhat <- summ_sub$beta; names(bhat) <- summ_sub$SNP
  shat <- summ_sub$se; names(shat) <- summ_sub$SNP

  bed_row <- bed_all[gene_id == locus]
  if (nrow(bed_row)==0) next
  psi <- as.numeric(bed_row[, -(1:4), with=FALSE])
  var_y <- var(psi, na.rm=TRUE)

  common_snps <- Reduce(intersect, list(names(bhat), names(prior_vec)))
  gt_sub <- gt_num[common_snps, , drop=FALSE]
  gt_sub <- gt_sub[apply(gt_sub, 1, function(x) all(!is.na(x))), , drop=FALSE]
  if (nrow(gt_sub) < 2) next

  bhat <- bhat[rownames(gt_sub)]
  shat <- shat[rownames(gt_sub)]
  prior_weights <- prior_vec[rownames(gt_sub)]
  R <- cor(t(gt_sub))

  fitted <- tryCatch(
    susie_rss(bhat=bhat, shat=shat, prior_weights=prior_weights,
               var_y=var_y, R=R, L=10,
               estimate_prior_variance=TRUE,
               estimate_residual_variance=TRUE,
               n=nrow(gt_sub)),
    error=function(e){ message('  -> SuSiE failed for ', locus); return(NULL) }
  )
  if (is.null(fitted)) next

  fitted_rss_batch[[locus]] <- fitted
  if (length(fitted_rss_batch) >= batch_size) {
    saveRDS(fitted_rss_batch,
      file=file.path(output_dir, paste0('", cell, "_fitted_rss_batch_', batch_counter, '.rds')))
    fitted_rss_batch <- list()
    batch_counter <- batch_counter + 1
    gc()
  }

  var_names <- names(bhat)
  pip <- fitted$pip
  cs_list <- fitted$sets$cs
  cs_id <- cs_size <- cs_purity <- rep(NA, length(var_names))
  for (cid in seq_along(cs_list)) {
    idx <- cs_list[[cid]]
    cs_id[idx] <- names(cs_list)[cid]
    cs_size[idx] <- length(idx)
    cs_purity[idx] <- fitted$sets$purity[cid,1]
  }
  df <- data.frame(locus_id=locus, variant_id=var_names, pip=pip,
                   cs=cs_id, cs_size=cs_size, cs_purity=cs_purity,
                   stringsAsFactors=FALSE)
  out_all[[locus]] <- df
  out_max[[locus]] <- df[which.max(df$pip), ]

  rm(gt_sub, R, bhat, shat, prior_weights, fitted)
  gc()
}

if (length(fitted_rss_batch) > 0) {
  saveRDS(fitted_rss_batch,
    file=file.path(output_dir, paste0('", cell, "_fitted_rss_batch_', batch_counter, '.rds')))
}

fwrite(rbindlist(out_all),
       file=file.path(output_dir, paste0('", cell, "_susie_all_variants.tsv')), sep='\\t')
fwrite(rbindlist(out_max),
       file=file.path(output_dir, paste0('", cell, "_susie_max_pip.tsv')), sep='\\t')

message('Finished cell type: ", cell, "')
"
  )
  
  writeLines(script_content, script_file)
}
