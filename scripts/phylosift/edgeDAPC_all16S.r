#!/usr/bin/R
library('ade4')
library('adegenet')
library('gplots')

epcaresdir = Sys.getenv()['resepca']

prefix = Sys.getenv()['prefix16S']
coulpop = c('black', 'blue', 'green', 'purple', 'red', 'orange', 'grey')
names(coulpop) = c('Aeta', 'Agta', 'Batak', 'Zambal', 'Casigurani', 'Tagbanua', 'American')
coullif = c('blue', 'red', 'grey')
names(coullif) = c('HG', 'F', 'C')

lifeshort = c('HG', 'F') ; names(lifeshort) = c('Hunter-gatherer', 'Agriculturalist')

sampleref = read.table(paste(metaghomedir, 'Philippines_Sample_List.tab', sep='/'), sep='\t', header=T, stringsAsFactors=F)

edgediff = read.csv(paste(epcaresdir, paste(prefix, 'edgediff', sep='.'), sep='/'), row.names=1, header=F, sep=',')

individuals = as.factor(sapply(rownames(edgediff), function(x){
	stsp = strsplit(strsplit(as.character(x), split='\\.')[[1]][1], split='_')[[1]]
	if (substr(stsp[1], 1, 3)=='run'){ return(stsp[2])
	}else{ return(stsp[1])}
}))


populations = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Population[sampleref$Sample==x]) }else{ return("American") }}))
lifestyles = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(lifeshort[sampleref$Lifestyle][sampleref$Sample==x]) }else{ return("C") }}))
localities = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Locality[sampleref$Sample==x]) }else{ return("USA") }}))
controlWGS = which(sapply(strsplit(names(individuals), split='\\.'), function(x){ substr(x[1], 1, 3)=='VFD' | (length(strsplit(x[1], split='_')[[1]])>1 & substr(x[1], 1, 3)=='SRS')  }))


rownames(edgediff) = as.character(individuals)

pcaedge.noscale_center = dudi.pca(edgediff, scale=F, center=T,  nf=8, scannf=F)
s.class(pcaedge.noscale_center$li, xax=1, yax=2, fac=populations, col=coulpop[levels(populations)], clabel=0)
s.class(pcaedge.noscale_center$li, xax=3, yax=4, fac=populations, col=coulpop[levels(populations)], clabel=0)
s.class(pcaedge.noscale_center$li, xax=1, yax=2, fac=lifestyles, col=coullif[levels(lifestyles)], clabel=0)
s.class(pcaedge.noscale_center$li, xax=3, yax=4, fac=lifestyles, col=coullif[levels(lifestyles)], clabel=0)

## PCA
eigenvalues = pcaedge.noscale_center$eig/sum(pcaedge.noscale_center$eig)
# assess bias of WGS against 16S amplicon sequencing
pdf(paste(epcaresdir, paste(prefix, 'ePCA_bias_amplicon_vs_WGS.pdf', sep='.'), sep='/'), width=10, height=10)
for (pcaplan in (1:2)){
	a = ((pcaplan-1)*2)+1
	b = a + 1
	plot(x=pcaedge.noscale_center$li[lifestyles=='C',a], y=pcaedge.noscale_center$li[lifestyles=='C',b], col='grey',
	 xlab=paste('ePCA PC ', a, ' (', round(eigenvalues[a], digits=2)*100, '%)', sep=''), ylab=paste('ePCA PC ', b, ' (', round(eigenvalues[b], digits=2)*100, '%)', sep=''))
	points(x=pcaedge.noscale_center$li[lifestyles=='HG',a], y=pcaedge.noscale_center$li[lifestyles=='HG',b], col='blue')
	points(x=pcaedge.noscale_center$li[lifestyles=='F',a], y=pcaedge.noscale_center$li[lifestyles=='F',b], col='red')
	points(x=pcaedge.noscale_center$li[controlWGS,a], y=pcaedge.noscale_center$li[controlWGS,b], col='green')
	legend('topright', pch=1, col=c('grey', 'green', 'blue', 'red'), legend=c('amplicon controls', 'WGS controls', 'WGS hunter-gatherers', 'WGS farmers'))
}
dev.off()

## DAPC
for (npca in (1:4)*2){
	pdf(paste(epcaresdir, paste(prefix, 'dapc', npca, 'PC.pdf', sep='.'), sep='/'), width=10, height=10)
	dapcls = dapc(edgediff, lifestyles, n.pca=npca, n.da=2)
	scatter(dapcls, col=coullif[levels(lifestyles)])
	dev.off()
}

# assess bias of WGS against 16S amplicon sequencing
npca = 4
pdf(paste(epcaresdir, paste(prefix, 'dapc', npca, 'PC_biasamplicon_vs_WGS.pdf', sep='.'), sep='/'), width=10, height=10)
plot(x=dapcls$ind.coord[lifestyles=='C',1], y=dapcls$ind.coord[lifestyles=='C',2], col='grey')
points(x=dapcls$ind.coord[lifestyles=='HG',1], y=dapcls$ind.coord[lifestyles=='HG',2], col='blue')
points(x=dapcls$ind.coord[lifestyles=='F',1], y=dapcls$ind.coord[lifestyles=='F',2], col='red')
points(x=dapcls$ind.coord[controlWGS,1], y=dapcls$ind.coord[controlWGS,2], col='green')
legend('topright', pch=1, col=c('grey', 'green', 'blue', 'red'), legend=c('amplicon controls', 'WGS controls', 'WGS hunter-gatherers', 'WGS farmers'))
dev.off()

