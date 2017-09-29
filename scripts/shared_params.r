#!/usr/bin/Rscript --vanilla --silent

## graphical and metadata parameters for mibrobiome scripts 
## these specific parameters relate to the dataset of oral metagenomes from the Philippines

repohome = Sys.getenv('repohome')
# if not specified, assumes the script (or its parent) is executed at the root of the repository
if ( repohome == "" ){ repohome = getwd() }
sampleref = read.table(file.path(repohome, 'data/Philippines_Sample_List.tsv'), sep='\t', header=T, stringsAsFactors=F)

lifecat = c('Hunter-gatherer', 'Traditional farmer', 'Western control')
lifeshort = c('HG', 'TF', 'WC') ; names(lifeshort) = lifecat ; names(lifecat) = lifeshort
popcat = c('Aeta', 'Agta', 'Batak', 'Tagbanua', 'Zambal', 'Casigurani', 'American')
popshort = c('Ae', 'Ag', 'Ba', 'Ta', 'Za', 'Ca', 'WC') ; names(popshort) = popcat ; names(popcat) = popshort
loccat = c('Palawan_Mount', 'Luzon_Mount', 'Luzon_Coast', 'USA')
locshort = c('PM', 'LM', 'LC', 'WC') ; names(locshort) = loccat ; names(loccat) = locshort
batchcat = c(1:3, 'SRS', 'VFD')

pll = c('populations', 'lifestyles', 'localities', 'batches')
coulpop = c('black', 'slateblue', 'blue', 'violetred1', 'chartreuse3', 'goldenrod', 'red2')
names(coulpop) = popshort
coullif = c('royalblue', 'limegreen', 'red2')
names(coullif) = lifeshort
coulloc = c('orange', 'darkgreen', 'purple', 'red2')
names(coulloc) = locshort
lcoul = list(coulpop, coullif, coulloc, NULL) ; names(lcoul) = pll
pchlif = c(16, 17, 16, 17, 16, 17, 15)
names(pchlif) = names(coulpop)
pchpop = c(17, 17, 17, 16, 16, 16, 15)
names(pchpop) = names(coulpop)
pchloc = c(17, 17, 17, 16, 16, 16, 15)
names(pchloc) = names(coulpop)
lpch = list(pchpop, pchlif, pchloc, NULL) ; names(lpch) = pll

critcat = list(popcat, lifecat, loccat, batchcat) ; names(critcat) = pll

getIndividualFactors = function(indata=NULL, individual.labels=NULL){
	if (is.null(individual.labels)){
		individuals = colnames(indata) 
	}else{ 
		individuals = individual.labels 
	}
	populations = factor(sapply(individuals, function(x){ popshort[sampleref$Population][sampleref$Sample==x] }), levels=popshort)
	lifestyles = factor(sapply(individuals, function(x){ lifeshort[sampleref$Lifestyle1][sampleref$Sample==x] }), levels=lifeshort)
	localities = factor(sapply(individuals, function(x){ locshort[sampleref$Locality][sampleref$Sample==x] }), levels=locshort)
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
