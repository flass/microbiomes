#!/usr/bin/Rscript --vanilla

# graphical and metadata parameters for Philippnnies oral metagenome scripts

#~ cargs = commandArgs(trailingOnly=T)
#~ if (len(cargs)>0){
#~ 	metaghomedir = cargs[1]
#~ }else{
#~ 	metaghomedir = file.path(homedir, 'oral_metagenomes')
#~ }
# assumes 'metaghomedir' is defined in local environment calling this script

sampleref = read.table(file.path(metaghomedir, 'Philippines_Sample_List.tab'), sep='\t', header=T, stringsAsFactors=F)

pll = c('populations', 'lifestyles', 'localities', 'batches')
coulpop = c('black', 'slateblue', 'blue', 'violetred1', 'chartreuse3', 'goldenrod', 'red2')
names(coulpop) = c('Aeta', 'Agta', 'Batak', 'Tagbanua', 'Zambal', 'Casigurani', 'American')
coullif = c('royalblue', 'limegreen', 'red2')
names(coullif) = c('HG', 'TF', 'WC')
coulloc = c('orange', 'darkgreen', 'purple', 'red2')
names(coulloc) = c('Palawan_Mount', 'Luzon_Mount', 'Luzon_Coast', 'USA')
lcoul = list(coulpop, coullif, coulloc, NULL)
names(lcoul) = pll
pchlif = c(16, 17, 16, 17, 16, 17, 15)
names(pchlif) = names(coulpop)
pchpop = c(17, 17, 17, 16, 16, 16, 15)
names(pchpop) = names(coulpop)
pchloc = c(17, 17, 17, 16, 16, 16, 15)
names(pchloc) = names(coulpop)
lpch = list(pchpop, pchlif, pchloc, NULL)
names(lpch) = pll

lifeshort = c('HG', 'TF', 'WC') ; names(lifeshort) = c('Hunter-gatherer', 'Traditional farmer', 'Western control')

getIndividualFactors = function(indata=NULL, individual.labels=NULL){
	if (is.null(individual.labels)){
		individuals = colnames(indata) 
	}else{ 
		individuals = individual.labels 
	}
	populations = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Population[sampleref$Sample==x]) }else{ return("American") }}))
	lifestyles = factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(lifeshort[sampleref$Lifestyle1][sampleref$Sample==x]) }else{ return("WC") }}), levels=lifeshort)
	localities = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Locality[sampleref$Sample==x]) }else{ return("USA") }}))
	batches = as.factor(sapply(individuals, function(x){ 
		if (x %in% sampleref$Sample){ return(paste('run', sampleref$Run[sampleref$Sample==x], sep='')) 
		}else{ return(substr(x, 1, 3)) }
	}))
	names(populations) = names(lifestyles) = names(localities) = names(batches) = individuals
	criteria = list(populations, lifestyles, localities, batches)
	names(criteria) = pll
	return(criteria)
}
