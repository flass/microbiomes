#!/usr/bin/Rscript --vanilla --silent

## !!! This script requires a local PosgreSQL build of the NCBI Taxonomy database !!!
## Flat files available within the Phylosift package (better use the same version as used when running Phylosift), or at ftp.ncbi.nlm.nih.gov/pub/taxonomy

library('ade4')
library('MASS')
library('adegenet')
library('gplots')
library('ggplot2')
#~ library('vegan')	prefer to use `ecodist::mantel` rather than `vegan::mantel`
library('ecodist')
library('ordPens')

repohome = Sys.getenv('repohome')
# if not specified, assume script is run from the repository top folder
if ( repohome == "" ){ repohome = getwd() }

# folder where all data are stored
cargs = commandArgs(trailingOnly=T)
if (length(cargs)>0){
	metaghomedir = cargs[1]
}else{
	metaghomedir = '~/oral_metagenomes'
}
stopifnot((file.exists(metaghomedir) && file.info(metaghomedir)$isdir))

# load shared params and functions
source(file.path(repohome, 'scripts/shared_params.r'), local=TRUE)

# load GPS-related function
source(file.path(repohome, 'scripts/gpscoords.r'))

colorWheel = function(angle){
 R = (abs((angle + 180)%%360 - 180) / 180)
 G = (abs((angle + 300)%%360 - 180) / 180)
 B = (abs((angle + 60 )%%360 - 180) / 180)
 return(rgb(R,G,B))
}

#####
# define local folder structure

prefix = Sys.getenv()['prefix']
subprefix = Sys.getenv()['subprefix']
prefixrestricted = Sys.getenv()['prefixrestricted']
subprefixrestricted = Sys.getenv()['subprefixrestricted']
phylosiftmarkerdbdir = Sys.getenv()['refpkg']
phylosiftmarker = Sys.getenv()['marker']
guppyresdir = Sys.getenv()['resjplace']
epcaresdir = Sys.getenv()['resepca']

if (is.na(prefix)){ prefix = 'alledges.33samples' }
if (is.na(prefixrestricted)){ prefixrestricted = '33samples' }
if (is.na(subprefix)){ subprefix = 'alledges.24samples' }
if (is.na(subprefixrestricted)){ subprefixrestricted = '24samples' }
if (is.na(phylosiftmarkerdbdir)){ phylosiftmarkerdbdir = file.path(metaghomedir, 'phylosift_v1.0.1/markers') }
if (is.na(phylosiftmarker)){ phylosiftmarker = 'concat.updated' }
if (is.na(guppyresdir)){ guppyresdir = file.path(metaghomedir, 'STEP_05_guppy') }
if (is.na(epcaresdir)){ epcaresdir = file.path(guppyresdir, '33samples.epca_results_std_settings') }


## load the taxon names associated with the edges of the reference tree from NCBI Taxonomy database
## !!! requires a local PosgreSQL build of the NCBI Taxonomy database !!!
nfedgenames = sprintf( "%s/%s.annotated/%s.taxonmap.names", phylosiftmarkerdbdir, phylosiftmarker, phylosiftmarker)
if (!file.exists(nfedgenames)){
	library('DBI')
	library('RPostgreSQL')
	drv = dbDriver("PostgreSQL")
	con = dbConnect(drv, user="lassalle", dbname="taxonomy", host="localhost")
	edgenames = sapply(1:dim(edgediff)[2], function(e){
		edge = e - 1 # edges are numberred from 0
		enames = dbGetQuery(con, sprintf("SELECT name_txt FROM taxonomy.phylosift_concat_reference INNER JOIN taxonomy.names USING (tax_id) WHERE branch_id=%d AND name_class='scientific name';", edge))
		if (dim(enames)[1]>0) return(enames[1,]) else return("") 
	})
	write(edgenames, file=nfedgenames)
}else{
	edgenames = readLines(nfedgenames)
}

########
# question meaning of sign of edge diff values relative to definition (Evans & Matsen PLOS One 2013)
# cf. pplacer user forum thread: https://groups.google.com/forum/#!topic/pplacer-users/qn89jUL5U-A
########
# edgediff definition in Figure 3 caption in  (Evans & Matsen, PLOS One 2013): "For every edge of the tree, the difference is taken between the number of reads on the non-root side the number of reads on the root side"
# edgediff definition in pplacer/guppy splitify manual (https://matsen.github.io/pplacer/generated_rst/guppy_splitify.html#guppy-splitify, last accessed 9/06/2016): 
# "The (s,e) entry of this matrix is the difference between the distribution of mass on either side of edge e for the sample s. Specifically, it is the amount of mass on the distal (non-root) side of edge e minus the amount of mass on the proximal side." 
# hence, most edges being remote from the root (half the edges of a tree are leaves), there should be a huge majority of edgediff values being close to -1, as most distal parts carry negligible mass ~ 0, and thus the remaining proximal part carry almost all the mass ~1 
# thus we should have: E(mass(distal) - mass(proximal)) ~ (0 - 1) = -1
# indeed most masses have value 1, the opposite of the expected value given the manual
#### however, Methods section of PLOS One paper describes it the other way: "Removing a given internal edge e from the tree splits T into two components: 
# T+(e) containing the root and T-(e) without. For a probability measure P on T, define the corresponding edge mass difference (= edgediff) as:
# dP(e) = P(T+(e)) - P(T-(e))"
# this results in less positive edgediff values in samples enriched for the clade 
####
# NB: while PCA eigenvectors may have either signs (oposite vectors are always equivalent solutions),
# using ade4 implementation of PCA provides eigenvectors which sign is indicative of greater values of the variable in the samples with positive values of the eigen vector.
# thus, we choose to have less negative edgediff values in samples enriched for the clade, by taking the opposite value of guppy output as input of the PCA
edgediff = -1 * read.csv(paste(epcaresdir, paste(prefix, 'edgediff', sep='.'), sep='/'), row.names=1, header=F, sep=',')
#~ hist(data.matrix(edgediff))

# change sign of edges which are on the path from arbitrary root to ToL root, which otherwise would result in aberant inversion of signal for these edges
invertededges = read.table(sprintf('%s/alledges.33samples.old2newrootpath.invertededges', epcaresdir), sep='\t', head=F, colClasses=c('numeric', 'character'))
print(unique(edgenames[invertededges[,1]+1]))
for (i in invertededges[,1]+1){
	edgediff[,i] = edgediff[,i] * -1
}

guppy.epca.proj = read.csv(paste(epcaresdir, paste(prefixrestricted, 'proj', sep='.'), sep='/'), row.names=1, header=F, sep=',')

individuals = as.factor(sapply(rownames(edgediff), function(x){
	stsp = strsplit(strsplit(as.character(x), split='\\.')[[1]][2], split='_')[[1]]
	if (substr(stsp[1], 1, 3)=='run'){ return(stsp[2])
	}else{ return(stsp[1])}
}))

rownames(edgediff) = as.character(individuals)

# all samples' metadata (from loaded script 'shared_params.r')
criteria = getIndividualFactors(individual.labels=individuals)
attach(criteria)

## diferent datasets
filipinos = lifestyles!='WC'
nrindiv = !(individuals %in% c('SRS015055', 'SRS013942'))	# HMP samples pairs (SRS014468, SRS015055) and (SRS019120, SRS013942) are from the same individuals at different time points
samplesets = list(filipinos, nrindiv, rep(TRUE, length(individuals)))
names(samplesets) = c('philippines.dataset.24samples', 'meta.analysis.dataset.31samples', 'meta.analysis.dataset.33samples')

# geolocalization
gps.locs = read.table(file.path(repohome, "data/philippines.24samples.coordinates.tab"), head=T, sep='\t')
gps.locs$dec.long = apply(gps.locs[paste('long', c('deg', 'min', 'sec'), sep='.')], 1, minutesseconds2decimal.coords)
gps.locs$dec.lat = apply(gps.locs[paste('lat', c('deg', 'min', 'sec'), sep='.')], 1, minutesseconds2decimal.coords)

filpop = popcat[popcat!='American']
average.coords = lapply(filpop, function(pop){
	average.gps.coords(gps.locs[gps.locs$label==pop, paste('dec', c('long', 'lat'), sep='.')], rad=F)
})
names(average.coords) = filpop
average.coords[['Tagbanua']] = average.coords[['Batak']]
print(average.coords)

for (namsam in names(samplesets)){
cat(sprintf("sample set: %s\n", namsam))
sampleset = samplesets[[namsam]]
dir.create(file.path(epcaresdir, namsam), showWarnings=F)
print(individuals[sampleset])
write(as.character(individuals[sampleset]), ncol=1, file=file.path(epcaresdir, namsam, "sample.list"))

K = length(which(sampleset)) - 1
#~ npca = K
npca = 4
ntopedges = 20
# dudi.pca + lda <=> dapc
# PCA
pcaedge.scale_center = dudi.pca(edgediff[sampleset,], scale=T, center=T,  nf=npca, scannf=F)
pcaedge.noscale_center = dudi.pca(edgediff[sampleset,], scale=F, center=T,  nf=npca, scannf=F) # the one corresponding to guppy's native edgePCA
#~ pcas = list(pcaedge.scale_center, pcaedge.noscale_center) ; names(pcas) = c('scaled_abundances', 'abundance-weighted')
pcas = list(pcaedge.noscale_center) ; names(pcas) = c('abundance-weighted')

for (pcasca in names(pcas)){
	print(pcasca)
	if (pcasca=='scaled_abundances'){ scalingedgevect = 750
	}else{ if (pcasca=='abundance-weighted'){ scalingedgevect = 5 
	}else{ scalingedgevect = 1 }}
	pcali = pcas[[pcasca]]$li
	pcaco = pcas[[pcasca]]$c1
	pcaeig = pcas[[pcasca]]$eig
	pcaeig.percent = round(100*pcaeig/sum(pcaeig), digits=1)
	cat(sprintf('%d PCs, %% variance explained: %s\n', npca, paste(pcaeig.percent[1:npca], collapse=' ')))
	nfpdfsam = file.path(epcaresdir, namsam, paste('ePCA', pcasca, 'pdf', sep='.'))
	pdf(nfpdfsam, width=20, height=20)
	for (k in 1:length(criteria)){
		print(names(criteria)[k])
		crit = droplevels(criteria[[k]][sampleset])
		for (npc in 1:(npca%/%2)){
			xax=((npc-1)*2)+1 ; yax=npc*2 ; pcs = c(xax, yax)
			print(pcs)
			subtitle = sprintf("ePCA, axis %d (%g%%) and %d (%g%%)\ngrouped by %s", xax, pcaeig.percent[xax], yax, pcaeig.percent[yax], names(criteria)[k])
			## basic scatter plot
			scatter(pcas[[pcasca]], xax=xax, yax=yax, clab.col=0, sub=subtitle, possub="bottomleft")
			
			## just project the individuals as dots
			s.class(pcali, xax=xax, yax=yax, fac=crit, col=lcoul[[k]][levels(crit)],
			 cellipse=0, cstar=0, cpoint=3, clabel=0, grid=F, axesell=F, pch=lpch[[k]][as.character(populations)], sub=subtitle, possub="bottomleft")

			## build the list of edges with the most significant variation on these axis
			# take the union of top edges at both extremities of the spectrum in both axes
			topedges = union(c(head(order(pcaco[,xax]), n=ntopedges), head(order(pcaco[,xax], decreasing=T), n=ntopedges)),
			 c(head(order(pcaco[,yax]), n=ntopedges), head(order(pcaco[,yax], decreasing=T), n=ntopedges)))
			## samples elispsis and points
			s.class(pcali, xax=xax, yax=yax, fac=crit, col=lcoul[[k]][levels(crit)], sub=subtitle, possub="bottomleft",
			 csub=3, cellipse=1, cpoint=3, clabel=2.5, grid=F, axesell=F, pch=lpch[[k]][as.character(populations)]
			)
			## samples elispsis and points + top variables vectors and labels
			s.class(pcali, xax=xax, yax=yax, fac=crit, col=lcoul[[k]][levels(crit)], sub=subtitle, possub="bottomleft",
			 csub=3, cellipse=1, cpoint=3, clabel=2.5, grid=F, axesell=F, pch=lpch[[k]][as.character(populations)]
			)
			s.arrow(pcaco[topedges,]*scalingedgevect, xax=xax, yax=yax, label=edgenames[topedges], add.plot=T)
			
			## samples elispsis + selected top variables vectors (coloured) and labels
			s.class(pcali, xax=xax, yax=yax, fac=crit, col=lcoul[[k]][levels(crit)], sub=subtitle, possub="bottomleft", 
			 csub=3, cellipse=1, cpoint=0, clabel=2.5, grid=F, axesell=F, cstar=0
			)
			
			# selects only the highest diverging edge within those sharing one same name
			nrtopedges = sapply(unique(edgenames[topedges]), function(en){
				samenameedges = topedges[edgenames[topedges]==en]
				# which edge vector with max norm
				samenameedges[which.max(sqrt(pcaco[samenameedges,xax]^2 + pcaco[samenameedges,yax]^2))]
			})
			
			x = pcaco[nrtopedges,xax]*scalingedgevect
			y = pcaco[nrtopedges,yax]*scalingedgevect
			arrows(0, 0, x, y, clabel=NULL,
			 lwd=10, angle=10, col=sapply(1:length(nrtopedges), function(e){ colorWheel(360*atan2(x[e], y[e])/(2*pi)) }))
	#~ 		s.label(pcaco[nrtopedges,]*scalingedgevect, xax=xax, yax=yax, label=edgenames[nrtopedges], add.plot=T, clabel=3)
			for (e in 1:length(nrtopedges)){
				text(x[e], y[e], label=edgenames[nrtopedges][e], adj=as.integer(c(x[e]<0, y[e]<0)), cex=2)
			}
			 
		}
	}
	dev.off()
	cat(sprintf("ePCA graphics in PDF file: '%s'\n", nfpdfsam))
	
	lifestyle = droplevels(lifestyles[sampleset])
	locality = droplevels(localities[sampleset])
	for (i in 1:4){
		cat(sprintf("ePCA, axis %d:\n", i))
		pc = pcali[,i]
		cat(sprintf("ANOVA: PC%d  ~ lifestyles\n", i))
		print(summary(aov(lm(pc  ~ lifestyle))))
		cat(sprintf("ANOVA: PC%d  ~ lifestyles * localities\n", i))
		print(summary(aov(lm(pc  ~ lifestyle * locality))))
		cat(sprintf("ANOVA: PC%d  ~ localities * lifestyles\n", i))
		print(summary(aov(lm(pc  ~ locality * lifestyle))))
#~ 		print(ordPens::ordAOV(as.integer(lifestyles[sampleset]), pc))
		cat(sprintf("ordered ANOVA: PC%d  ~ lifestyles\n", i))
		print(ordAOV(as.integer(lifestyle), pc))
	}
	# write R-based PC table for projection on tree
	write.table(t(pcaco), file=file.path(epcaresdir, namsam, paste('ePCA', pcasca, npca, 'PC.csv', sep='.')), sep=',', col.names=F, row.names=T, quote=T)

	# write RGB components + width associates with each edge, based on pairs of PCs

	for (npc in 1:2){
		xax=((npc-1)*2)+1 ; yax=npc*2
		colwidtable = apply(pcaco[,c(xax, yax)], 1, function(xy){ c(colorWheel(360*atan2(xy[1], xy[2])/(2*pi)), sqrt(sum(xy^2))) })
		write.table(t(colwidtable), file=paste(epcaresdir, sprintf('ePCA-%s.PC%g%g.colourwidth.csv', pcasca, xax, yax), sep='/'), sep='\t', col.names=F, row.names=F, quote=F)
	}

	## LDA
	cat(sprintf("LDA: microbiome PC1-%d ~ lifestyles\n", npca))
	crit = 'lifestyles'
	disc.pcali = lda(as.matrix(pcali), grouping=lifestyle)
	# test significacne of DA
	cat(sprintf("MANOVA: microbiome PC1-%d ~ lifestyles\n", npca))
	manova.pcali = manova(as.matrix(pcali) ~ lifestyle)
	print(summary(manova.pcali, test='Pillai'))
	print(summary(manova.pcali, test='Roy'))

	# PCA column normed scores (principal axis ??) * linear coefficients for transformation of PC axis into LD axis
	nda = nlevels(lifestyle)-1
	print(disc.pcali$scaling)
	disc.pcali.varvect = sapply(1:(nda), function(ld){ rowSums(sapply(1:4, function(k){ pcaco[,k] * disc.pcali$scaling[k,ld] })) })
	colnames(disc.pcali.varvect) = colnames(disc.pcali$scaling)
	write.table(t(disc.pcali.varvect), file=file.path(epcaresdir, namsam, paste('dapc', pcasca, crit, npca, 'PC.var.vect.csv', sep='.')), sep=',', col.names=F, row.names=T, quote=T)

	# DAPC
	cat(sprintf("DAPC: microbiome PC1-%d ~ lifestyles\n", npca))
	dapcls = dapc(edgediff[sampleset,], lifestyle, n.pca=npca, n.da=nda)
	print(summary(dapcls))
	print(dapcls$loadings) # linear combination of PC axis into LD axis
	write.table(t(dapcls$var.contr), file=file.path(epcaresdir, namsam, paste('dapc', pcasca, crit, npca, 'PC.var.contr.csv', sep='.')), sep=',', col.names=F, row.names=T, quote=T)

	
#~ 	## calculate the part of total variance explained by DA axis
#~ 	print(RVAideMemoire::MVA.synt(disc.pcali))
	npdfdapc = file.path(epcaresdir, namsam, paste('dapc', pcasca, npca, 'PC.pdf', sep='.'))
	pdf(npdfdapc, width=15, height=15)
	scatter(dapcls, col=coullif[colnames(dapcls$posterior)], pch=pchlif[as.character(populations[sampleset])], sub=paste('DAPC, considering', npca, 'axis'))
	if (nda > 2){
		loadingplot(dapcls$var.contr, axis=2)
	}	
	bp = barplot(t(dapcls$posterior), beside=F, col=coullif[colnames(dapcls$posterior)], names.arg=rep('', length(individuals[sampleset])), main="Posterior Probabilities of individuals belonging to a group")
	legend('topright', colnames(dapcls$posterior), fill=coullif)
	mtext(text=individuals[sampleset], at=bp, side=1, line=1, col=coullif[lifestyle], las=2)
	dev.off()
	npdfdapcmini = sub('PC.pdf', 'PC.mini.pdf', npdfdapc, fixed=T)
	pdf(npdfdapcmini, width=5, height=5)
	scatter(dapcls, col=coullif[colnames(dapcls$posterior)], pch=pchlif[as.character(populations[sampleset])])	
	dev.off()
	cat(sprintf("edgeDAPC graphics in PDF files:\n'%s'\n'%s'\n", npdfdapc, npdfdapcmini))

	## correlation of microbiome distances with host genetic distances
	genettag = 'FULL'
	hostgendist = as.dist(read.table(file.path(repohome, 'data/philippines.24samples.genetic.dist'), skip=2, row.names=1))
	# using all data
	genetindiv = attr(hostgendist, 'Labels')
	nphi = length(genetindiv)
	totalmicrobdist = list(dist(edgediff[genetindiv,], method="euclidean")) ; names(totalmicrobdist) = 'microbiome.total'
	# using only 1 of the first 4 PC axis
	onePCmicrobdist = lapply(1:4, function(pc){
		as.dist(sapply(genetindiv, function(x){
			sapply(genetindiv, function(y){
				abs(pcali[x,pc] - pcali[y,pc])
			})
		}))
	})
	names(onePCmicrobdist) = paste('microbiome.PC', 1:4, sep='.')
	# using combinations of 2 of the first 4 PC axis
	twoPCmicrobdist = unlist(lapply(1:3, function(pca){
		ld = lapply((pca+1):4, function(pcb){
			as.dist(sapply(genetindiv, function(x){
				sapply(genetindiv, function(y){
					sqrt(sum((pcali[x,pca] - pcali[y,pca])^2 + (pcali[x,pcb] - pcali[y,pcb])^2))
				})
			}))
		})
		names(ld) = paste('microbiome.PC', pca, (pca+1):4, sep='.')
		return(ld)
	}), recursive=F)
	
	allmicrobdist = c(totalmicrobdist, onePCmicrobdist, twoPCmicrobdist)
	
	## taking into account geography
	geographicdist = as.dist(sapply(genetindiv, function(ind1){
		i = which(individuals[sampleset]==ind1)
		gd = sapply(genetindiv, function(ind2){
			j = which(individuals[sampleset]==ind2)
				d = geodetic.distance(average.coords[[popcat[as.character(populations[sampleset][i])]]], average.coords[[popcat[as.character(populations[sampleset][j])]]])
		})
		names(gd) = genetindiv
		return(gd)
	}))

	# test influence of host genetics
	mantelgenmicrobdists = t(sapply(allmicrobdist, function(microbdist){ 
		ecodist::mantel(microbdist ~ hostgendist, nperm=9999, mrank=F)
#~ 		vegan::mantel(xdis=microbdist, ydis=hostgendist, permutations=9999, method=tolower(mantelmethod))
	}))

	# test influence of geography
	mantelgeomicrobdists = t(sapply(allmicrobdist, function(microbdist){ 
		ecodist::mantel(microbdist ~ geographicdist, nperm=9999, mrank=F)
#~ 		vegan::mantel(xdis=microbdist, ydis=geographicdist, permutations=9999, method=tolower(mantelmethod))
	}))

	# test influence of host genetics removing variance explained by geography (partial Mantel's test)
	partialmantelgenmicrobdists = t(sapply(allmicrobdist, function(microbdist){ 
		ecodist::mantel(microbdist ~ hostgendist + geographicdist, nperm=9999, mrank=F)
#~ 		vegan::mantel.partial(xdis=microbdist, ydis=hostgendist, zdis=geographicdist, permutations=9999, method=tolower(mantelmethod))
	}))	

	# test influence of geography removing variance explained by host genetics (partial Mantel's test)
	partialmantelgeomicrobdists = t(sapply(allmicrobdist, function(microbdist){ 
		ecodist::mantel(microbdist ~ geographicdist + hostgendist, nperm=9999, mrank=F)
#~ 		vegan::mantel.partial(xdis=microbdist, ydis=hostgendist, zdis=geographicdist, permutations=9999, method=tolower(mantelmethod))
	}))	
	
	write.table(mantelgenmicrobdists, file=file.path(epcaresdir, namsam, paste('Mantel', genettag, 'GenetDistvsPCA', pcasca, 'tab', sep='.')), sep='\t')
	write.table(mantelgeomicrobdists, file=file.path(epcaresdir, namsam, paste('Mantel', genettag, 'GeographyDistvsPCA', pcasca, 'tab', sep='.')))
	write.table(partialmantelgenmicrobdists, file=file.path(epcaresdir, namsam, paste('partialMantel', genettag, 'GenetDistMinusGeographyDistvsPCA', pcasca, 'tab', sep='.')), sep='\t')
	write.table(partialmantelgeomicrobdists, file=file.path(epcaresdir, namsam, paste('partialMantel', genettag, 'GeographytDistMinusGeographyDistvsPCA', pcasca, 'tab', sep='.')), sep='\t')
	cat(sprintf("Mantel test tables written in files: '%s'\n", file.path(epcaresdir, namsam, paste('Mantel', genettag, '*', pcasca, 'tab', sep='.'))), sep='\t')
	
	hgd = as.numeric(hostgendist)
	gd = as.factor(as.numeric(geographicdist))
	plotdf = do.call(rbind, lapply(names(allmicrobdist), function(x){
		md = as.numeric(allmicrobdist[[x]])
		data.frame(microbiome.dist=md, host.genetic.dist=hgd, geographic.dist=gd, variation=as.factor(rep(x, length(md))))
	}))
	print(summary(plotdf))
	pdf(file.path(epcaresdir, namsam, paste('edgePCA', pcasca, 'correlations', genettag, 'pdf', sep='.')), width=10, height=10)
	p = ggplot(data = plotdf, aes(host.genetic.dist, microbiome.dist)) + geom_point() + stat_smooth(method="lm", se=T)
	print(p + facet_wrap( ~ variation, nrow=4))
#~ 	p = ggplot(data = plotdf, aes(geographic.dist, microbiome.dist)) + geom_violin(aes(fill = geographic.dist), trim = FALSE, adjust = .75)
	p = ggplot(data = plotdf, aes(geographic.dist, microbiome.dist)) + geom_point() + stat_smooth(method="lm", se=T)
	print(p + facet_wrap( ~ variation, nrow=4))
	dev.off()

}

# PCA made by guppy uses a diffrent projection! Also eigen vectors direction are not expected to have a particular meaning here
# The signs of vectors corresponding to the lineage of the Gamma-proteobacteria leading to Enterobacter would be wrong anyway because 
# it does not take into account the reversal of sign on lineage when re-rooting from arbitrary root in Enterobacter sp. leaf to the root of the Tree of Life. 

pdf(file.path(epcaresdir, namsam, paste('GUPPY-projected.ePCA.pdf', sep='.')), width=15, height=15)
for (k in 1:length(criteria)){
	crit = droplevels(criteria[[k]][sampleset])
	for (npc in 1:2){
		xax=((npc-1)*2)+1 ; yax=npc*2
		print(c(xax, yax))
		s.class(guppy.epca.proj[sampleset,], xax=xax, yax=yax, fac=crit,
		 col=lcoul[[k]][levels(crit)], pch=lpch[[k]][as.character(populations[sampleset])])
	}
}
dev.off()

}
