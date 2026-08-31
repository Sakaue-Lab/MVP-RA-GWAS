# Genome-wide genetic correlation (LDSC-rg) between RA and immune-mediated diseases.
# Uses our EUR RA meta-analysis summary statistics vs. a panel of 15 publicly
# available LDSC-formatted EUR immune-mediated disease sumstats.
# HapMap3 SNPs, MHC excluded, 1000G Phase 3 EUR LD scores.

ldsc=/path_to/ldsc/ldsc.py
munge=/path_to/ldsc/munge_sumstats.py
LDSC=/path_to/LDSCORE

# --- 1. munge our RA (EUR) sumstats ---
RA=meta_analysis_results_EURMVP_b37_updated
python ${munge} \
  --sumstats stats/${RA}.sumstats.gz \
  --merge-alleles ${LDSC}/w_hm3.noMHC.snplist \
  --chunksize 500000 \
  --out stats/${RA}.munged

# --- 2. immune-mediated disease sumstats (already LDSC-formatted, EUR) ---
# one *.sumstats.gz per trait under stats/immune/ ; see traits.txt (15 traits)
TRAITS=$(cat stats/immune/traits.txt | tr '\n' ',' | sed 's/,$//')
RG_TRAITS=""
for t in $(cat stats/immune/traits.txt); do
  RG_TRAITS="${RG_TRAITS},stats/immune/${t}.sumstats.gz"
done

# --- 3. cross-trait rg (RA is the first trait; all pairwise rg vs RA) ---
${ldsc} \
  --rg stats/${RA}.munged.sumstats.gz${RG_TRAITS} \
  --ref-ld-chr ${LDSC}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --w-ld-chr   ${LDSC}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --out stats/RA_immune_rg

# --- 4. collect rg / se / p from the LDSC log into a summary table ---
Rscript collect_rg.R stats/RA_immune_rg.log stats/rg_summary.csv
