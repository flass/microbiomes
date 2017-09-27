#! /bin/bash
# remove the hierarchical whitespace/tab indenting from Kraken abundance summary files
# usage: parsing-Kraken_summaries.sh /path/to/kraken_summary_folder
indir=$1
for nf in `ls $indir/*summary.txt` ; do
nfr=${nf%.*}
sed 's/  \+\|^ /\t/g' $nf | sed 's/\t\t/\t/g' | sed 's/^\t//g' > $nfr.tab
done
