#!/usr/bin/Rscript --vanilla

## graphical and metadata parameters for mibrobiome scripts 
## these specific parameters relate to the dataset of oral metagenomes from the Philippines

# assumes the script (or its parent) is executed at the root of the repository
sampleref = read.table(file.path('data', 'Philippines_Sample_List.tsv'), sep='\t', header=T, stringsAsFactors=F)

lifecat = c('Hunter-gatherer', 'Traditional farmer', 'Western control')
lifeshort = c('HG', 'TF', 'WC') ; names(lifeshort) = lifecat
popcat = c('Aeta', 'Agta', 'Batak', 'Tagbanua', 'Zambal', 'Casigurani', 'American')
popshort = c('Ae', 'Ag', 'Ba', 'Ta', 'Za', 'Ca', 'WC') ; names(popshort) = popcat
loccat = c('Palawan_Mount', 'Luzon_Mount', 'Luzon_Coast', 'USA')
locshort = c('PM', 'LM', 'LC', 'WC') ; names(locshort) = loccat

pll = c('populations', 'lifestyles', 'localities', 'batches')
coulpop = c('black', 'slateblue', 'blue', 'violetred1', 'chartreuse3', 'goldenrod', 'red2')
names(coulpop) = popcat
coullif = c('royalblue', 'limegreen', 'red2')
names(coullif) = lifeshort
coulloc = c('orange', 'darkgreen', 'purple', 'red2')
names(coulloc) = loccat
lcoul = list(coulpop, coullif, coulloc, NULL) ; names(lcoul) = pll
pchlif = c(16, 17, 16, 17, 16, 17, 15)
names(pchlif) = names(coulpop)
pchpop = c(17, 17, 17, 16, 16, 16, 15)
names(pchpop) = names(coulpop)
pchloc = c(17, 17, 17, 16, 16, 16, 15)
names(pchloc) = names(coulpop)
lpch = list(pchpop, pchlif, pchloc, NULL) ; names(lpch) = pll


getIndividualFactors = function(indata=NULL, individual.labels=NULL){
	if (is.null(individual.labels)){
		individuals = colnames(indata) 
	}else{ 
		individuals = individual.labels 
	}
#~ 	populations = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Population[sampleref$Sample==x]) }else{ return("American") }}))
#~ 	lifestyles = factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(lifeshort[sampleref$Lifestyle1][sampleref$Sample==x]) }else{ return("WC") }}), levels=lifeshort)
#~ 	localities = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Locality[sampleref$Sample==x]) }else{ return("USA") }}))
	populations = as.factor(sapply(individuals, function(x){ sampleref$Population[sampleref$Sample==x] }))
	lifestyles = factor(sapply(individuals, function(x){ lifeshort[sampleref$Lifestyle1][sampleref$Sample==x] }), levels=lifeshort)
	localities = as.factor(sapply(individuals, function(x){ sampleref$Locality[sampleref$Sample==x] }))
	batches = as.factor(sapply(individuals, function(x){ 
		r = sampleref$Run[sampleref$Sample==x]
		if (!is.na(r)){ return(paste('run', r, sep='')) 
		}else{ return(substr(x, 1, 3)) }
	}))
	strategy = ordered(lifeshort[lifestyles])
	names(populations) = names(lifestyles) = names(localities) = names(batches) = individuals
	criteria = list(populations, lifestyles, localities, batches)
	names(criteria) = pll
	return(criteria)
}
