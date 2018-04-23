#/bin/bash
krakenres="$HOME/oral_metagenomes/STEP_03_Kraken_results/kraken_classification/vs_krakendb_30-09-17"
repohome="$HOME/github_flass/microbiomes"
lefseprog="$HOME/Python_lib/lefse"
ldathresh=3

cd $krakenres/
for f in `ls ./*.summary.txt` ; do
  # from a Kraken report: select only species; filter out zero-read rows; filter out Homo spiens (taxon_id=9606)
  grep -P "\tS\t" $f | grep -vP "^ +0.00\t0\t0\t" | grep -vP "\tS\t9606\t" | sed -e 's/\t */\t/g' > ${f}_presentSonly_nohuman
done

for fil in filtered_020 unfiltered; do
  mkdir -p ./${fil}.reports/
  rm ./${fil}.reports/*.txt
  ls $PWD/*.${fil}.summary.txt_presentSonly_nohuman | while read f ; do
    bn=$(basename $f)
    case $bn in
    run*)
      sample=$(echo ${bn%%.*} | cut -d'_' -f2)
      ;;
    SRS*)
      sample=$(echo ${bn%%.*} | cut -d'_' -f1)
      ;;
    VFD*)
      sample=${bn%%.*}
      ;;
    *)
      echo "unexpected sample: ${bn}"
      ;;
    esac 
    ln -s $f ./${fil}.reports/${sample}.txt
  done
  frad="$PWD/${fil}.reports/33_oral_samples_vs_krakendb_30-09-17.${fil}.presentSonly_nohuman"
  inf=($(ls $PWD/${fil}.reports/*.txt))
  kraken-biom ${inf[@]} -o ${frad}.biom
  biom add-metadata -i ${frad}.biom -o ${frad}.biom_metadata -m ${repohome}/data/sample.mapping.tsv
  biom convert -i ${frad}.biom_metadata -o ${frad}.4lefse-0.tsv --to-tsv
  ${repohome}/scripts/kraken/sub_taxid_hierarchical_str.py ${frad}.4lefse-0.tsv ${frad}.4lefse-1.tsv
  for criterion in lifestyle locality ; do
    fradc=${frad}_${criterion}
    cat ${repohome}/data/lefse_metadata_header_${criterion} > ${fradc}.4lefse.tsv
    tail -n+3 ${frad}.4lefse-1.tsv >> ${fradc}.4lefse.tsv
    ${lefseprog}/format_input.py ${fradc}.4lefse.tsv ${fradc}.lefse -c 1 -s 2 -o 1000000
    ${lefseprog}/run_lefse.py ${fradc}.lefse ${fradc}.lefse_res -o ${fradc}.lefse_out -l ${ldathresh}
    ${lefseprog}/plot_res.py ${fradc}.lefse_res ${fradc}.lefse_lda.pdf --format pdf
    ${lefseprog}/plot_res.py ${fradc}.lefse_res ${fradc}.lefse_lda.pdf --format svg
    ${lefseprog}/plot_cladogram.py ${fradc}.lefse_res ${fradc}.lefse_clado.pdf --format pdf
    ${lefseprog}/plot_features.py ${fradc}.lefse ${fradc}.lefse_res ${fradc}.lefse_feat.pdf --format pdf -f one --feature_name "Campylobacter consisus"
  done
  rm ./${fil}.reports/*.4lefse*.tsv
done
