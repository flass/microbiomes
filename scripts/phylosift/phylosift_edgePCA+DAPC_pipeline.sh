#!/bin/bash

### pipeline to perform DAPC (Jombart, T., Devillard, S. & Balloux, F. Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genet. 11, 94 (2010).)
### on edge difference matrix from pplacer/guppy (Matsen, F. A. & Evans, S. N. Edge Principal Components and Squash Clustering: Using the Special Structure of Phylogenetic Placement Data for Sample Comparison. PLoS ONE 8, e56859 (2013).)
### as used in Phylosift (Darling, A. E. et al. PhyloSift: phylogenetic analysis of genomes and metagenomes. PeerJ 2, e243 (2014).)

# fill in relevant paths to software
export repohome='/path/to/repository/microbiomes'
export reposcripts="$repohome/scripts/phylosift"
export phylosift='/path/to/phylosift/installation'
export bins='/path/to/pplacer/binaries'

# fill in relevant paths to data
export refpkg=$phylosift/phylosift_v1.0.1/markers
export resjplace='/path/to/guppy/output/jplace_files'
export resepca='/path/to/output/of/this/script'

# dataset names
export prefixrestricted='33samples'
export prefix="alledges.$prefixrestricted"
export subprefixrestricted='24samples'
export subprefix="alledges.$subprefixrestricted"

## get placement mass diferences across edge (do not prune non-significant edges)
$bins/guppy splitify -o $prefixrestricted.edgediff --epsilon 0 --kappa 1 --prefix alledges. --out-dir $resepca $resjplace/*.jplace

## prepare annotated tree in reference package
export marker='concat.updated'

# copy Phylosift reference package of marker genes 'concat.updated' (i.e. concatenate of 33 universal genes)
if [ -e $refpkg/$marker.annotated ] then
  rm -r $refpkg/$marker.annotated
fi
cp -r $refpkg/$marker $refpkg/$marker.annotated

python $reposcripts/annotate_phylosift_reference_tree.py $refpkg $resjplace

## build local PosgreSQL database using NCBI Taxonomy  dumps provided with Phylosift installation

psql < $reposcripts/taxonomy.sql

## list edges for which to correct edge mass difference signs because of re-rooting to Tree of Life root
python $reposcripts/list_edges_inverted_by_rerooting.py $resepca

## perform edge PCA and DAPC (taking into account edge mass difference sign reversal for flagged edges)
R --no-save --no-restore < $reposcripts/edgePCA.r

for pcasca in 'scaled_abundances' 'abundance-weighted' ; do
  for dataset in `ls $resepca/*.dataset.* -d` ; do
    ## tree representation
    # represent original variable (absolute) contribution to DAPC LD axis on phyloXML tree
    $bins/guppy heat -c $refpkg/$marker.annotated -o $resepca/$dataset/$prefix.dapc.$pcasca.var.contr.4.PC.xml --min-fat 0 $resepca/$dataset/$prefix.dapc.4.PC.var.contr.csv
    # represent original variable coordinates in DAPC LD axis on phyloXML tree
    $bins/guppy heat -c $refpkg/$marker.annotated -o $resepca/$dataset/$prefix.dapc.$pcasca.var.vect.4.PC.xml --min-fat 0 $resepca/$dataset/$prefix.dapc.4.PC.var.vect.csv

    # represent R-computed PCA on phyloXML tree
    $bins/guppy heat -c $refpkg/$marker.annotated -o $resepca/$dataset/$prefix.ePCA.$pcasca.4.PC.xml --min-fat 0 $resepca/$dataset/$prefix.ePCA.$pcasca.4.PC.csv

    # represent R-computed PCA on phyloXML tree using in-house colors (color wheel depicting PC-pair planes)
	python $reposcripts/colour_reftree_as_PC_plane.py $respca $pcasca/$dataset
  done
done
