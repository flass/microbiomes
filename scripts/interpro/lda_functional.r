#!/share/apps/R/bin/R
library('ade4')
library('ggplot2')
#~ library('ordPens')
library('clinfun')
library('parallel')

nbcores = 4

options(width = 160)

cargs = commandArgs(trailingOnly=T)
nffuncmat = cargs[1]
if (is.na(nffuncmat)){
	
funcanaldir = '/path/to/functional/metagenomic/dataset'
nffuncmat = sprintf('%s/all_interpro_absolute_counts.tab', funcanaldir)
inmat = read.table(nffuncmat, header=T, row.names=1, sep='\t', quote='"')
funcmat = data.matrix(t(inmat[,-1]))
relfuncmat = t(apply(funcmat, 1, function(x){ x / sum(x) }))
description = as.character(inmat[,1]) ; names(description) = rownames(inmat)

nfnewfuncmat = sprintf('%s/all_ip_term_read_counts.csv', funcanaldir)
newfuncmat = data.matrix(read.csv(nfnewfuncmat, header=T, row.names=1))
relnewfuncmat = t(apply(newfuncmat, 2, function(x){ x / sum(x) }))

rfm = relnewfuncmat ; nffm = nfnewfuncmat

}else{
	funcmat = data.matrix(read.csv(nffuncmat, header=T, row.names=1))
	relfuncmat = t(apply(funcmat, 2, function(x){ x / sum(x) }))
	rfm = relfuncmat ; nffm = nffuncmat
}

sampleref = read.table(sprintf('%s/Philippines_Sample_List.tab', funcanaldir), header=T, sep='\t')
stratshort = c('HG', 'AG', 'WC') ; names(stratshort) = c('Hunter-gatherer', 'Agriculturalist', 'Western control')
strategy = ordered(sapply(rownames(rfm), function(s){
  if ((u<-strsplit(s, split='_')[[1]][2]) %in% sampleref[['Sample']]) stratshort[as.character(sampleref[sampleref[['Sample']]==u, 'Lifestyle'])] else 'WC'
 }), levels=stratshort)

# do a PCA to make dudi object
func.pca = dudi.pca(rfm, scannf=F, nf=ncol(rfm), scale=F)

S = combn(nlevels(strategy),2)
layout(matrix(1:(ncol(S)*2), ncol(S), 2, byrow=T))

topdiscvars = list()
topsignifdiscvars = list()
topvars = list()
for (i in 1:ncol(S)){
	s = S[,i]
	
	stratpair = levels(strategy)[s]
	pastepair = paste(stratpair, collapse='vs')
	print(stratpair)
	stratpair.i = which(strategy %in% stratpair)
	strat = droplevels(strategy[stratpair.i])
#~ 	print(strat) ; print(length(strat))
	pair.rfm = data.matrix(rfm[stratpair.i,])
#~ 	print(dim(pair.rfm))
	# do a PCA to make dudi object
	pair.func.pca = dudi.pca(pair.rfm, scannf=F, nf=nrow(pair.rfm))
#~ 	print(func.pca)
	# check how it is
	s.class(pair.func.pca$li, fac=strat, xax=1, yax=2)	# PC1,2
	s.class(pair.func.pca$li, fac=strat, xax=3, yax=4)	# PC3,4
	disc.func = discrimin(pair.func.pca, strat, scannf=F, nf=1)
	topdiscvar = disc.func$va[order(abs(disc.func$va)[,1], decreasing=T),1]
#~ 	print(topdiscvar[1:24])
	
	# perform alll independent t-tests on species abundances
	ttestpvals = apply(pair.rfm, 2, function(func){
		tt = t.test(func ~ strat)
		return(tt$p.value)
	})
	ranked.ttestpvals = ttestpvals[order(ttestpvals)]
	# Benjamini–Hochberg procedure
	signifthresh = 0.05
	m = length(ranked.ttestpvals)
	BHcor.ttestpvals = sapply(1:m, function(k){ ranked.ttestpvals[k] * m / k })
	topsignif = which(BHcor.ttestpvals <= signifthresh)
	topsignifdisc = order(ttestpvals)[topsignif]
	topsignifdiscvar = ttestpvals[topsignifdisc]
	
	outnames = names(topdiscvar)
	vardf = data.frame(functional.term=outnames, LDA.loading=topdiscvar[outnames], t.test.p=ttestpvals[outnames], t.test.FDRcor.p=BHcor.ttestpvals[outnames], description=description[outnames])
	print(vardf[1:30,])
	write.table(vardf, sprintf("%s_LDA_%s", nffm, pastepair), sep='\t', row.names=F)
	
	topvars[[pastepair]] = 	vardf	
	topdiscvars[[pastepair]] = topdiscvar
	topsignifdiscvars[[pastepair]] = topsignifdiscvar

}

# sum the loadings (sum th signed loading so opposed varaitions will cancel each other)
summary.pw.lda = as.data.frame(t(simplify2array(mclapply(1:dim(rfm)[2], function(i){
	ipterm = colnames(rfm)[i]
	cat(sprintf("\r%d", i))
	loadings = sapply(topdiscvars, `[`, ipterm)	
	ranks.lda = sapply(topdiscvars, function(x, y){which(names(x)==y)}, ipterm)	
	sum.load = sum(loadings)
	sum.rank = sum(ranks.lda)
	return(c(loadings, sum.load, ranks.lda, sum.rank))
}, mc.cores=nbcores))))
rownames(summary.pw.lda) = colnames(rfm)
colnames(summary.pw.lda) = c(paste('loading.lda', names(topdiscvars), sep='.'), 'summed.loadings', paste('rank.lda', names(topdiscvars), sep='.'), 'summed.ranks')
summary.pw.lda$description = description[rownames(summary.pw.lda)]

for (functerm in c('Ribosomal protein', 'CRISPR-associated', '([mM]ultidrug|antimicrobial) (resistance|extrusion|efflux)', 'Ketopantoate reductase')){
	t = summary.pw.lda[grep(functerm, summary.pw.lda$description), c(grep('loading', colnames(summary.pw.lda), value=T), 'description')]
	print(functerm)
	print(t[order(abs(t$summed.loadings), decreasing=T),])
}

if (length(cargs)>1){
	for (carg in cargs[-1]){
		candidatepathway = read.table(carg, sep='\t', header=F, row.names=1)
		ipterms = intersect(rownames(candidatepathway), colnames(rfm))
		pdf(paste(carg, 'counts_per_strategy.pdf', sep='.'), height=12, width=12)
		ecs = unique(as.character(candidatepathway[ipterms,1]))
		wrelabun = rowMeans(sapply(ecs, function(ec){
			# summed read abundance for all terms relating to one reaction 
			annots = ipterms[ecs==ec]
			if (length(annots)>1) rs = rowSums(rfm[, annots])
			else rs = rfm[, annots]
			# normalize per reaction
			return(rs/mean(rs))
		}))
		boxplot(wrelabun ~ strategy, ylab='weighted relative read abundance for combined pathway reactions', main=paste(ipterms, collapse=' '))
		layout(matrix(1:9,3,3,byrow=T))
		for (ipterm in ipterms){
			boxplot(rfm[, ipterm] ~ strategy, ylab='relative read abundance for functional annotation',
			 main=paste(paste(ipterm, as.character(candidatepathway[ipterm,1]), sep=' - '), description[ipterm], sep='\n'))
		}
		dev.off()
	}	
}
