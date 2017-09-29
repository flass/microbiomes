# Copyright Matteo Fumagalli, 2015
#
# this is to process bam file (whole and subsample) and then perform PCA and admixture

# FULL and SUB data sets

# PATHS

NGSTOOLS=/data/Software/ngsTools
ANGSD=$NGSTOOLS/angsd
NGSADMIX=/data/data/Software/NGSadmix/NGSadmix

REF=/data/data/hg19/chr1.fa.gz

# get percentiles for cutoff on global depth
MININD=24
for TYPE in FULL SUB;
do
	echo $TYPE
	$ANGSD/angsd -b Data/$TYPE.bamlist -ref $REF -out Results/$TYPE.qc -P 4 \
		-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
		-minMapQ 20 -minInd $MINDIND \
		-doQsDist 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -r chr1: # &> /dev/null
	
	Rscript $NGSTOOLS/Scripts/plotQC.R Results/$TYPE.qc # &> /dev/null

done

# PCA

TYPE=FULL
MINDEPTH=90
MAXDEPTH=310

$ANGSD/angsd -P 4 -b Data/$TYPE.bamlist -ref $REF -out Results/$TYPE \
	-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minQ 20 -minMapQ 20 -minInd $MININD -setMinDepth $MINDEPTH -setMaxDepth $MAXDEPTH -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-4 -r chr1: \
	-doGeno 32 -doPost 1

less -S Results/$TYPE.mafs.gz
N_SITES=`zcat Results/$TYPE.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
# 573024

N_SAMPLES=`wc -l Data/$TYPE.bamlist | cut -d " " -f 1`
echo $N_SAMPLES
# 24

gunzip Results/$TYPE.geno.gz

$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/$TYPE.geno -outfile Results/$TYPE.covar -nind $N_SAMPLES -nsites $N_SITES -call 0 -norm 0 &> /dev/null

Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/$TYPE.covar -c 1-2 -a pops.clst -o Results/$TYPE.pca12.pdf
Rscript $NGSTOOLS/Scripts/plotPCA.R -i Results/$TYPE.covar -c 3-4 -a pops.clst -o Results/$TYPE.pca34.pdf

evince Results/$TYPE.pca12.pdf
evince Results/$TYPE.pca34.pdf

TYPE=SUB
MINDEPTH=40
MAXDEPTH=150

$ANGSD/angsd -P 4 -b Data/$TYPE.bamlist -ref $REF -out Results/$TYPE \
	-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minQ 20 -minMapQ 20 -minInd $MININD -setMinDepth $MINDEPTH -setMaxDepth $MAXDEPTH -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-4 -r chr1: \
	-doGeno 32 -doPost 1

less -S Results/$TYPE.mafs.gz
N_SITES=`zcat Results/$TYPE.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES
# 510111 

N_SAMPLES=`wc -l Data/$TYPE.bamlist | cut -d " " -f 1`
echo $N_SAMPLES
# 24

gunzip Results/$TYPE.geno.gz

$NGSTOOLS/ngsPopGen/ngsCovar -probfile Results/$TYPE.geno -outfile Results/$TYPE.covar -nind $N_SAMPLES -nsites $N_SITES -call 0 -norm 0 &> /dev/null

# open R

# procrustes analyses to rotate SUB coordinates

library(vegan)

full=as.matrix(read.table("Results/FULL.covar", head=F, stringsAsFact=F, sep="\t"))[1:24,1:24]
eig_full <- eigen(full, symm=TRUE)$vec;

sub=as.matrix(read.table("Results/SUB.covar", head=F, stringsAsFact=F, sep="\t"))[1:24,1:24]
eig_sub <- eigen(sub, symm=T)$vec

proc<-protest(eig_full[,1:4], eig_sub[,1:4], scale=T, scores = "sites", permutations = 99999)

protest(eig_full[,1:4], matrix(runif(24*4),ncol=4,nrow=24), scale=T, scores = "sites", permutations = 9999)

cat(proc$signif)
cat(proc$ss)

#Procrustes Sum of Squares (m12 squared):        0.4678 
#Correlation in a symmetric Procrustes rotation: 0.7295 
#Significance:  1e-04 

#Permutation: free
#Number of permutations: 9999

# close R

## Admixture

NGSADMIX=/data/data/Software/NGSadmix/NGSadmix

TYPE=FULL
MINDEPTH=30
MAXDEPTH=300

# old 30, 300, keep these

$ANGSD/angsd -P 4 -b Data/$TYPE.bamlist -ref $REF -out Results/$TYPE \
	-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minQ 20 -minMapQ 20 -minInd $MININD -setMinDepth $MINDEPTH -setMaxDepth $MAXDEPTH -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-4 -r chr1: \
	-doGlf 2 -doPost 1 

for K in 2 3 4 5 6
do
	echo $K
	$NGSADMIX -likes Results/FULL.beagle.gz -K $K -outfiles Results/FULL.admix.$K -minMaf 0.02 -P 4 &> /dev/null
done

Rscript plotAdmix.R


## GENETIC DISTANCES

TYPE=FULL
MINDEPTH=30
MAXDEPTH=300

$ANGSD/angsd -P 4 -b Data/$TYPE.bamlist -ref $REF -out Results/$TYPE \
	-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
	-minMapQ 20 -minQ 20 -minInd $MININD -setMinDepth $MINDEPTH -setMaxDepth $MAXDEPTH -doCounts 1 \
	-GL 1 -doMajorMinor 1 -doMaf 1 -skipTriallelic 1 \
	-SNP_pval 1e-4 -r chr1: \
	-doGeno 8 -doPost 1 &> /dev/null

N_SITES=`zcat Results/$TYPE.mafs.gz | tail -n+2 | wc -l`
echo $N_SITES # 604154

$NGSTOOLS/ngsDist/ngsDist -verbose 1 -geno Results/$TYPE.geno.gz -probs -n_ind $N_SAMPLES -n_sites $N_SITES -labels ../pops.label -o Results/$TYPE.dist -n_threads 4 &> /dev/null
less -S Results/$TYPE.dist








