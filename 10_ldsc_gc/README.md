# 10_ldsc_gc

Genome-wide genetic correlation between RA and immune-mediated diseases using
cross-trait LD score regression (LDSC-rg).

We used our European-ancestry RA meta-analysis summary statistics and a panel of
publicly available, LDSC-formatted EUR GWAS summary statistics for 15
immune-mediated diseases, restricting to HapMap3 variants and excluding the MHC
region. LD scores were computed from the 1000 Genomes Project phase 3 European
reference panel. Statistical significance was assessed after Bonferroni
correction for the number of traits tested (p < 0.05/15).

## Scripts
- `ldsc_rg.sh` : munge the RA sumstats and run `ldsc.py --rg` of RA against the
  15 immune-mediated disease sumstats.
- `collect_rg.R` : parse the LDSC `--rg` log into `stats/rg_summary.csv`.
- `plot_ldsc_rg.R` : one-column heatmap of rg with a p-value column
  (`Rscript plot_ldsc_rg.R stats/rg_summary.csv fig/rg_heatmap.png`).

## Notes
- Replace `/path_to/...` with your local paths (LDSC install and `LDSCORE` data).
- `stats/immune/traits.txt` lists the 15 immune-mediated traits (one per line);
  each `stats/immune/<trait>.sumstats.gz` is an LDSC-formatted EUR sumstat.
