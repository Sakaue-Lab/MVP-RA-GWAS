# munge sumstats

for categ in meta_analysis_results_allMVP_b37_updated meta_analysis_results_EURMVP_b37_updated
do
python ${munge} \
--sumstats stats/${categ}.HESS.chi2_qc.input.gz \
--merge-alleles stats/w_hm3.noMHC.snplist \
--ignore beta --chunksize 500000 \
--out stats/${categ}.munged
done

# run LDSC
ETH=EUR

for categ in meta_analysis_results_allMVP_b37_updated meta_analysis_results_EURMVP_b37_updated
do
 ${ldsc} \
 --h2 stats/${categ}.munged.sumstats.gz \
 --ref-ld-chr ../../../external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baselineLD_v2.1_ldscores/baselineLD. \
 --frqfile-chr ../../../external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq/1000G.${ETH}.QC. \
 --w-ld-chr ../../../external-data/LDSCORE/data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
 --overlap-annot \
 --thin-annot \
 --print-coefficients \
 --print-delete-vals \
 --out stats/${categ}
done

# collect
for categ in meta_analysis_results_allMVP_b37_updated meta_analysis_results_EURMVP_b37_updated
do
FILE=stats/${categ}.log
H2=`cat ${FILE} | grep "Total Observed scale h2" | awk '{print $5}'`
H2SE=`cat ${FILE} | grep "Total Observed scale h2" | awk '{print $6}' | sed 's/(//' | sed 's/)//'`
MEANCHI=`cat ${FILE} | grep "Mean Chi" | awk '{print $3}'`
INT=`cat ${FILE} | grep "Intercept" | awk '{print $2}'`
INTSE=`cat ${FILE} | grep "Intercept" | awk '{print $3}' | sed 's/(//' | sed 's/)//'`
RATIO=`cat ${FILE} | grep "Ratio" | awk '{print $2}'`
echo -e "${categ}\t${H2}\t${H2SE}\t${MEANCHI}\t${INT}\t${INTSE}\t${RATIO}" >> stats/MVP.h2.summary_02152026.txt
done