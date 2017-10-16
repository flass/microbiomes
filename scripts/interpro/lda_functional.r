#!/usr/bin/Rscript --vanilla
library('ade4')
library('ggplot2')
#~ library('ordPens')
#~ library('clinfun')
library('parallel')
library('ancom.R')


repohome = Sys.getenv('repohome')
# if not specified, assume script is run from the repository top folder
if ( repohome == "" ){ repohome = getwd() }
source(file.path(repohome, 'scripts/shared_params.r'))

nbcores = 4

options(width = 160)

interestingterms = c('Ribosomal protein', 'CRISPR-associated', '([mM]ultidrug|antimicrobial) (resistance|extrusion|efflux)', 'Ketopantoate reductase')

# folder where all data are stored
cargs = commandArgs(trailingOnly=T)
if (length(cargs)>0){
	nffuncmat = cargs[1]
}else{
	nffuncmat = '~/oral_metagenomes/functional_analysis/ERP016024_IPR_abundances_v3.0.tsv'
}
if (file.exists(nffuncmat)){
	ffuncmat = file(nffuncmat)
}else{
	ebiurl = 'https://www.ebi.ac.uk/metagenomics/projects/ERP016024/download/3.0/export?contentType=text&exportValue=IPR_abundances'
	cat(sprintf("Local file '%s' does not exist.\nInstead, download dataset from '%s'.\n", nffuncmat, ebiurl))
	ffuncmat = file(ebiurl)
	nffuncmat = file.path(getwd(), 'ERP016024_IPR_abundances_v3.0.tsv') # to define output file names
}
inmat = read.table(ffuncmat, header=T, row.names=1, sep='\t', quote='"')
funcmat = data.matrix(t(inmat[,-1]))
relfuncmat = t(apply(funcmat, 1, function(x){ x / sum(x) }))
description = as.character(inmat[,1]) ; names(description) = rownames(inmat)

individuals = sapply(rownames(relfuncmat), function(err){ as.character(sampleref[sampleref[['EBI_Metagenomics_Run_ID']]==err, 'Sample']) })
	
# all samples' metadata (from loaded script 'shared_params.r')
criteria = getIndividualFactors(individual.labels=individuals)
attach(criteria)
strategy = ordered(lifeshort[lifestyles])
names(strategy) = individuals

pdf(sub('.tsv', '_counts_per_strategy.pdf', nffuncmat, fixed=T), height=6, width=12)
layout(matrix(1:2, 1, 2))

# do a PCA to make dudi object
func.pca = dudi.pca(relfuncmat, scannf=F, nf=ncol(relfuncmat), scale=F)
func.pca.percenteig = round(100*func.pca$eig/sum(func.pca$eig), digit=1)
for (pcs in list(1:2, 3:4, 5:6)){
  s.class(func.pca$li, fac=strategy, xax=pcs[1], yax=pcs[2], col=coullif, axesell=F, grid=F, sub=paste(paste('PC', pcs, ':', func.pca.percenteig[pcs]), collapse='; '))
  s.class(func.pca$li, fac=populations, xax=pcs[1], yax=pcs[2], col=coulpop, axesell=F, grid=F, sub=paste(paste('PC', pcs, ':', func.pca.percenteig[pcs]), collapse='; '))
}


S = combn(nlevels(strategy),2)
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
	stpop = droplevels(populations[stratpair.i])
	pair.relfuncmat = data.matrix(relfuncmat[stratpair.i,])
	# find non constant variables
	nocons = which(simplify2array(mclapply(1:ncol(pair.relfuncmat), function(j){ var(pair.relfuncmat[,j]) > 0 }, mc.cores=nbcores, mc.preschedule=T)))
#~ 	print(head(t(pair.relfuncmat)))
	# do a PCA to make dudi object
	pair.func.pca = dudi.pca(pair.relfuncmat[,nocons], scannf=F, nf=nrow(pair.relfuncmat))
	pair.func.pca.percenteig = round(100*pair.func.pca$eig/sum(pair.func.pca$eig), digit=1)
	for (pcs in list(1:2, 3:4, 5:6)){
	  print(pcs)
	  s.class(pair.func.pca$li, fac=strat, xax=pcs[1], yax=pcs[2], col=coullif[levels(strat)], axesell=F, grid=F, sub=paste(paste('PC', pcs, ':', func.pca.percenteig[pcs]), collapse='; '))
	  s.class(pair.func.pca$li, fac=stpop, xax=pcs[1], yax=pcs[2], col=coulpop[levels(stpop)], axesell=F, grid=F, sub=paste(paste('PC', pcs, ':', func.pca.percenteig[pcs]), collapse='; '))
	}
	disc.func = discrimin(pair.func.pca, strat, scannf=F, nf=1)
	topdiscvar = disc.func$va[order(abs(disc.func$va)[,1], decreasing=T),1]
#~ 	print(topdiscvar[1:24])
	
	# perform alll independent t-tests on species abundances
	ttestpvals = simplify2array(mclapply(nocons, function(j){
		tt = t.test(pair.relfuncmat[,j] ~ strat)
		return(tt$p.value)
	}, mc.cores=nbcores, mc.preschedule=T))
	names(ttestpvals) = colnames(relfuncmat)[nocons]
	ranked.ttestpvals = ttestpvals[order(ttestpvals)]
	cat('ranked.ttestpvals:', head(ranked.ttestpvals), '...\n')
	
#~ 	# ANCOM test (Mandal S et al. (2015). Analysis of composition of microbiomes: a novel method for studying microbial composition. Microbial Ecology in Health and Disease, 26, 1-7.)
#~ 	pair.relfuncmat.strat = as.data.frame(cbind(pair.relfuncmat[,nocons], as.character(strat)))
#~ 	print(dim(pair.relfuncmat.strat))
#~ 	print(head(pair.relfuncmat.strat[,(ncol(pair.relfuncmat.strat)-10):ncol(pair.relfuncmat.strat)]))
#~ 	ancom.result = ANCOM(pair.relfuncmat.strat, multcorr=2, ncore=nbcores)
#~ 	plot_ancom(ancom.result)
#~ 	write(ancom.result$detected, file=sprintf("%s_ANCOM_signif_diff_relabun_%s", nffuncmat, pastepair))
#~ 	cat('ancom.result$detected:', head(ancom.result$detected), '...\n')
	# computation does not return within a day; also does not use the specified multiple cores
	
	# Benjamini–Hochberg procedure
	signifthresh = 0.05
	m = length(ranked.ttestpvals)
	BHcor.ttestpvals = simplify2array(mclapply(1:m, function(k){ ranked.ttestpvals[k] * m / k }, mc.cores=nbcores, mc.preschedule=T))
	topsignif = which(BHcor.ttestpvals <= signifthresh)
	topsignifdisc = order(ttestpvals)[topsignif]
	topsignifdiscvar = ttestpvals[topsignifdisc]
	cat('topsignifdiscvar:', head(topsignifdiscvar), '...\n')
	
#~ 	outnames = names(topdiscvar[1:30])
	outnames = names(topdiscvar)
	vardf = data.frame(functional.term=outnames, LDA.loading=topdiscvar[outnames], t.test.p=ttestpvals[outnames], t.test.FDRcor.p=BHcor.ttestpvals[outnames], description=description[outnames])
	print(vardf[1:30,])
	write.table(vardf, sprintf("%s_LDA_%s", nffuncmat, pastepair), sep='\t', row.names=F)
	
	topvars[[pastepair]] = 	vardf	
	topdiscvars[[pastepair]] = topdiscvar
	topsignifdiscvars[[pastepair]] = topsignifdiscvar
}

dev.off()

# sum the loadings (sum th signed loading so opposed varaitions will cancel each other)
summary.pw.lda = as.data.frame(t(simplify2array(mclapply(1:dim(relfuncmat)[2], function(i){
	ipterm = colnames(relfuncmat)[i]
	cat(sprintf("\r%d", i))
	loadings = sapply(topdiscvars, `[`, ipterm)	
	ranks.lda = sapply(topdiscvars, function(x, y){which(names(x)==y)}, ipterm)	
	sum.load = sum(loadings)
	sum.rank = sum(ranks.lda)
#~ 	oaov = ordAOV(as.numeric(strategy), relfuncmat[, ipterm], nsim=1e+06)
	jtt = jonckheere.test(relfuncmat[,ipterm], as.integer(strategy), nperm=1e+06, alternative="two.sided")
#~ 	return(c(loadings, sum.load, ranks.lda, sum.rank, oaov$p.value))
	return(c(loadings, sum.load, ranks.lda, sum.rank, jtt$p.value))
}, mc.cores=nbcores))))
rownames(summary.pw.lda) = colnames(relfuncmat)
colnames(summary.pw.lda) = c(paste('loading.lda', names(topdiscvars), sep='.'), 'summed.loadings', paste('rank.lda', names(topdiscvars), sep='.'), 'summed.ranks', 'jonckheere.test.p.value')
summary.pw.lda$description = description[rownames(summary.pw.lda)]
save(summary.pw.lda, file=sprintf("%s_summary_LDAandJonckheere.RData", nffuncmat))
#~ print(head(summary.pw.lda[order(abs(summary.pw.lda$summed.loadings), decreasing=T), c(grep('loading', colnames(summary.pw.lda), value=T), 'description', 'jonckheere.test.p.value')], n=30))
#~ print(head(summary.pw.lda[order(summary.pw.lda$summed.ranks,         decreasing=F), c(grep('rank',    colnames(summary.pw.lda), value=T), 'description')], n=30))

for (functerm in interestingterms){
	t = summary.pw.lda[grep(functerm, summary.pw.lda$description), c(grep('loading', colnames(summary.pw.lda), value=T), 'jonckheere.test.p.value', 'description')]
	print(functerm)
	print(t[order(abs(t$summed.loadings), decreasing=T),])
}

if (length(cargs)>1){
	for (carg in cargs[-1]){
		candidatepathway = read.table(carg, sep='\t', header=F, row.names=1)
		ipterms = intersect(rownames(candidatepathway), colnames(relfuncmat))
		pdf(paste(carg, 'counts_per_strategy.pdf', sep='.'), height=12, width=12)
		ecs = unique(as.character(candidatepathway[ipterms,1]))
		wrelabun = rowMeans(sapply(ecs, function(ec){
			# summed read abundance for all terms relating to one reaction 
			annots = ipterms[ecs==ec]
			if (length(annots)>1) rs = rowSums(relfuncmat[, annots])
			else rs = relfuncmat[, annots]
			# normalize per reaction
			return(rs/mean(rs))
		}))
		boxplot(wrelabun ~ strategy, ylab='weighted relative read abundance for combined pathway reactions', main=paste(ipterms, collapse=' '))
# 		ggplot(data.frame(wrelabun=wrelabun, strategy=strategy), aes(strategy, wrelabun))+geom_boxplot()+ylab('weighted relative read abundance for combined pathway reactions')+ggtitle(paste(ipterms, collapse=' '))+geom_jitter()
		layout(matrix(1:9,3,3,byrow=T))
		for (ipterm in ipterms){
			boxplot(relfuncmat[, ipterm] ~ strategy, ylab='relative read abundance for functional annotation',
			 main=paste(paste(ipterm, as.character(candidatepathway[ipterm,1]), sep=' - '), description[ipterm], sep='\n'))
		}
# 		ggplot(as.data.frame(cbind(relfuncmat[, ipterms], strategy)), aes(strategy))+geom_boxplot()+ylab('relative read abundance for functional annotation')+ggtitle(paste(paste(ipterm, as.character(candidatepathway[ipterm,1]), sep=' - '), description[ipterm], sep='\n'))+facet_grid(.~ipterms)
		dev.off()
	}	
}
