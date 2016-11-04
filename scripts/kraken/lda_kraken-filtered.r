#!/usr/bin/Rscript
# removes quotes
library('ade4')
library('parallel')
library('gplots')
library('RColorBrewer')
library('MASS')


plotPCAabun = function(tabun.pca, criterion, couls, subsample=NULL, np=4, topdisc=NULL, scalingvarvect=1, cex=1){
	pcali = tabun.pca$li
	pcaco = tabun.pca$c1
	if (is.null(subsample)){ obj=rep(T, length(criterion)) }else{ obj=subsample }
	for (i in 0:(np-1)){
		xax = (i*2)+1
		yax = (i+1)*2
		s.class(pcali, fac=criterion[obj], xax=xax, yax=yax, col=couls)
		s.label(pcali, xax=xax, yax=yax, add.plot=T)
		if (!is.null(topdisc)){ 
			if (scalingvarvect=='auto'){
				# get the longest norm of the sample and loading vectors
				maxsamnorm = max(apply(pcali[,c(xax,yax)], 1, function(x){ sqrt(x[1]^2 + x[2]^2) }))
				maxloadnorm = max(apply(pcaco[topdisc,c(xax,yax)], 1, function(x){ sqrt(x[1]^2 + x[2]^2) }))
				svv = maxsamnorm / maxloadnorm
			}else{ svv = scalingvarvect }
			x = pcaco[topdisc,xax]*svv
			y = pcaco[topdisc,yax]*svv
			arrows(x0=0, y0=0, x1=x, y1=y)
			for (e in 1:length(topdisc)){
				text(x[e], y[e], label=topdisc[e], adj=as.integer(c(x[e]<0, y[e]<0)), cex=cex)
			}
		}
	}
}

nbcores = 4
user = Sys.info()['user']
homedir = paste('/home', user, sep='/')
metaghome = paste(homedir, 'oral_metagenomes', sep='/')

#~ lifeshort = c('C', 'F', 'HG') ; names(lifeshort) = c('Western', 'Agriculturalist', 'Hunter-gatherer')
#~ lifeshort = c('WC', 'AG', 'HG') ; names(lifeshort) = c('Western', 'Agriculturalist', 'Hunter-gatherer')
lifeshort = c('HG', 'AG', 'WC') ; names(lifeshort) = c('Hunter-gatherer', 'Agriculturalist', 'Western Controls')

pll = c('populations', 'lifestyles', 'localities')
coulpop = c('black', 'slateblue', 'blue', 'violetred1', 'chartreuse3', 'goldenrod', 'red2')
names(coulpop) = c('Aeta', 'Agta', 'Batak', 'Tagbanua', 'Zambal', 'Casigurani', 'American')
#~ coullif = c('red2', 'limegreen', 'royalblue')
coullif = c('royalblue', 'limegreen', 'red2')
names(coullif) = lifeshort
#~ coullif = c('royalblue', 'red2', 'limegreen', 'slategrey')
#~ names(coullif) = c('A', 'C', 'F', 'HG')
coulloc = c('orange', 'darkgreen', 'red2', 'purple')
names(coulloc) = c('Palawan_Mount', 'Luzon_Mount', 'USA', 'Luzon_Coast')
lcoul = list(coulpop, coullif, coulloc)
names(lcoul) = pll
pchlif = c(16, 17, 16, 17, 16, 17, 15)
names(pchlif) = names(coulpop)
pchpop = c(17, 17, 17, 16, 16, 16, 15)
names(pchpop) = names(coulpop)
pchloc = c(17, 17, 17, 16, 16, 16, 15)
names(pchloc) = names(coulpop)
lpch = list(pchpop, pchlif, pchloc)
names(lpch) = pll

orederedpops.PC3 = c('Aeta', 'Batak', 'Agta', 'Tagbanua', 'Zambal', 'Casigurani', 'American')
#~ orederedlifs.PC4 = c('HG', 'F', 'C')
orederedlifs.PC4 = c('HG', 'AG', 'WC')
orederedlocs.PC3 = c('Palawan_Mount', 'Luzon_Mount', 'Luzon_Coast', 'USA')
ordered.crits = list(orederedpops.PC3, orederedlifs.PC4, orederedlocs.PC3) ; names(ordered.crits) = pll

#~ scalingvarvect = 'auto'
scalingvarvect = 60

#~ contamsample = 'Cae45'

taxlevels = c('G', 'S')
taxonomiclevels = c('Genus', 'Species') ; names(taxonomiclevels) = taxlevels
refranks = c('genera', 'species') ; names(refranks) = taxlevels
krakenhome = paste(metaghome, 'STEP_03_Kraken_results', sep='/')
coabundir = paste(metaghome, 'STEP_08_co-abundance', sep='/')
#~ krakenhome = '/cluster/scratch3/ucbtmsp/oral_microbiome_project/STEP_03_Kraken_results'
krakenrestopdir = paste(krakenhome, 'kraken_classification/microbial_db', sep='/')
krakenresdirs = paste(krakenrestopdir, c('reads_human_removed_duskmask_4', 'HMP_reads_duskmsk_4'), sep='/')

interestingtaxa = c('Veillonella', 'Prevotella', 'Neisseria', 'Streptococcus', 'Haemophilus', 'Actinomyces', 'Treponema', 'Megasphaera')
pathogenspecies = readLines(sprintf('%s/pathogens_Warinner2014+Chen2014.txt', coabundir))

filtertags = c('filtered.020', 'unfiltered')
outdirs = paste(krakenhome, c('LDA_confidence-filtered', 'LDA_sensitive'), sep='/') ; names(outdirs) = filtertags

for (filtertag in filtertags){
	print(filtertag)
	outdir = outdirs[[filtertag]]
	if (!file.exists(outdir)){
		dir.create(outdir)
	}

	nfcleankrakenabunmat = paste(outdir, 'kraken_cumulative_readcounts_non-normalized_cleaned-taxa.RData', sep='/')
	if (file.exists(nfcleankrakenabunmat)){
		print(sprintf('load matrix \'cleankrakenabun\' and \'ranks\' vector from file %s', nfcleankrakenabunmat))
		load(nfcleankrakenabunmat)
	}else{
		nfkrakenabuntab = paste(outdir, 'kraken_cumulative_readcounts_non-normalized_all-taxa.tsv', sep='/')
	#~ 	nfkrakenabuntab = paste(outdir, 'kraken_cumulative_readcounts_non-normalized_all-taxa_names-as-column.tsv', sep='/')
		if (file.exists(nfkrakenabuntab)){
			print(sprintf('load table \'krakenabun\' from file %s', nfkrakenabuntab))
			krakenabun = read.table(nfkrakenabuntab, sep='\t', head=T, row.names=1)
	#~ 		krakenabun = read.table(nfkrakenabuntab, sep='\t', head=T)
		}else{
			# kraken results from the 29/08/2014 :
			lnfkrakenpath = do.call(c, lapply(krakenresdirs, function(p){ 
				q = list.files(p)
				if (p=='HMP_reads_duskmsk_4'){ q = q[substr(q, 1, 3) %in% c('SRS', 'VFD')] }
				return(paste(p, q, sep='/')) 
			}))
			lnfkrakenres = basename(lnfkrakenpath)
			#~ lnfkrakenrad = t(sapply(lnfkrakenres, function(x){ y = strsplit(x, split='\\.')[[1]] ; return(c(y[1], paste(y[1:(length(y)-1)], collapse='.'), y[length(y)])) }))
			lnfkrakenrad = t(sapply(lnfkrakenpath, function(p){ x = basename(p) ; y = strsplit(x, split='\\.')[[1]] ; return(c(y[1], paste(y[2:(length(y)-1)], collapse='.'), y[length(y)])) }))
			colnames(lnfkrakenrad) = c('sample', 'file.rad', 'file.ext')
			tabfs = which(lnfkrakenrad[,'file.rad']==paste('kraken', filtertag, 'summary', sep='.') & lnfkrakenrad[,'file.ext']=='tab')
			lkrakenrestab = mclapply(tabfs, function(tabf){
				nfkrakentab = lnfkrakenres[tabf]
				print(nfkrakentab)
			#~ 	ktab = read.table(paste(krakenresdir, nfkrakentab, sep='/'), sep='\t', quote='')
				ktab = read.table(rownames(lnfkrakenrad)[tabf], sep='\t', quote='')
				colnames(ktab) = c('abundance', 'cummulative.nb.reads', 'node.specific.nb.reads', 'rank', 'taxid', 'name')
				return(ktab[order(ktab$taxid),])
			}, mc.cores=nbcores, mc.preschedule=T)

			individuals = as.factor(sapply(lnfkrakenres[tabfs], function(x){
				stsp = strsplit(strsplit(x, split='\\.')[[1]][1], split='_')[[1]]
				subs = substr(stsp[1], 1, 3)
				if (subs=='run'){ return(stsp[2])
				}else{ if (subs=='VFD'){ return(stsp[1]) 
				}else{ return(stsp[1])}}
			}))

			names(lkrakenrestab) = individuals
			# check the taxids are the same in the same order, has to be FALSE
			ordertaxidok = any(as.logical(apply(simplify2array(mclapply(names(lkrakenrestab), function(samplename){ lkrakenrestab[[samplename]]$taxid }, mc.cores=nbcores, mc.preschedule=T)), 1, function(x){length(unique(x))-1})))
			if (!ordertaxidok){ print('same tax_id order in all abundance tables')
			}else{ stop('not the same tax_id order in all abundance tables') }
		
			#~ krakenabun = simplify2array(mclapply(names(lkrakenrestab), function(samplename){ lkrakenrestab[[samplename]]$abundance }, mc.cores=nbcores, mc.preschedule=T))
			#~ rownames(krakenabun) = lkrakenrestab[[1]]$name
			#~ colnames(krakenabun) = names(lkrakenrestab)
			#~ write.table(krakenabun, file=paste(outdir, 'kraken_cumulative_abundances_non-normalized_all-taxa.tsv', sep='/'), sep='\t')
			krakenabun = simplify2array(mclapply(names(lkrakenrestab), function(samplename){
				lkrakenrestab[[samplename]][['cummulative.nb.reads']]
			}, mc.cores=nbcores, mc.preschedule=T))
			rownames(krakenabun) = lkrakenrestab[[1]]$taxid
			colnames(krakenabun) = names(lkrakenrestab)
			write.table(krakenabun, file=nfkrakenabuntab, sep='\t')	# with non-redundant row names
			write(lkrakenrestab[[1]]$name, file=paste(nfkrakenabuntab, "rownames", sep='.'))
			rownames(krakenabun) = lkrakenrestab[[1]]$name
		}
		# remove taxa with no match in any sample
		presenttax = apply(krakenabun, 1, sum)>0
		# remove abundance counts from human contamination
		humanlineage = lkrakenrestab[[1]]$taxid %in% read.table(paste(metaghome, 'human_lineage_taxid.tab', sep='/'), sep='\t', header=T)$taxid
		
		samplespecifictaxa = lapply(colnames(krakenabun), function(i){ krakenabun[,i]>0 & rowSums(krakenabun[,colnames(krakenabun)!=i])==0 }) #, mc.cores=nbcores, mc.preschedule=T)
		names(samplespecifictaxa) = colnames(krakenabun)
		print(sapply(samplespecifictaxa, function(x){length(which(x))}))
		
		# deal with (apparently) soil-contaminated sample
	#~ 	contaminant = krakenabun[,contamsample]>0 & rowSums(krakenabun[,colnames(krakenabun)!=contamsample])==0
#~ 		contaminant = samplespecifictaxa[[contamsample]]
		# has many RHizobiales and pathogenic Enterobacteriaceae (Eschecrichia, Salmonella, Yersinia)
		# try keeping it for the moment
		#~ 
		cleankrakenabun = krakenabun[presenttax & !humanlineage, ]
		ranks = lkrakenrestab[[1]]$rank[presenttax & !humanlineage]
	#~ 	cleankrakenabun = krakenabun[presenttax & !contaminant & !humanlineage, ]
	#~ 	ranks = lkrakenrestab[[1]]$rank[presenttax & !contaminant & !humanlineage]
		save(cleankrakenabun, ranks, file=nfcleankrakenabunmat)
		# clean up memory space
		rm(lkrakenrestab, krakenabun)
		gc()
	}

	sampleref = read.table(paste(metaghome, 'Philippines_Sample_List.tab', sep='/'), sep='\t', header=T, stringsAsFactors=F)
	individuals = colnames(cleankrakenabun)
	populations = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Population[sampleref$Sample==x]) }else{ return("American") }})) ; names(populations) = individuals
	lifestyles = factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(lifeshort[sampleref$Lifestyle][sampleref$Sample==x]) }else{ return("WC") }}), levels=lifeshort) ; names(lifestyles) = individuals
	#~ lifestyles = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(lifeshort[sampleref$Lifestyle2][sampleref$Sample==x]) }else{ return("C") }}))
	localities = as.factor(sapply(individuals, function(x){ if (x %in% sampleref$Sample){ return(sampleref$Locality[sampleref$Sample==x]) }else{ return("USA") }})) ; names(localities) = individuals

	criteria = list(populations, lifestyles, localities)
	names(criteria) = pll

	readsets = list(c('root', 'unclassified'), c('root'), c('Bacteria', 'Archaea'))
	names(readsets) = c('all', 'classified', 'Prokaryotes')

	depths = sapply(readsets, function(rs){
		sapply(colnames(cleankrakenabun), function(s){
			log10(sum(cleankrakenabun[rs, s]))
		})
	})
	depthset = 'Prokaryotes'
	couldepth = rainbow(20, start=0.66, end=1)
	depthrange = seq(from=min(depths[, depthset]), to=max(depths[, depthset]), length.out=21)
	couldepvec = sapply(depths[, depthset], function(x){ which(sapply(2:21, function(n){ x >= depthrange[n-1] & x <= depthrange[n] })) })


	for (refrank in names(refranks)){
		print(refranks[refrank])
		# normalize read counts at this taxonomic rank so sum(all counts) = 1
		normcleankrakenabun = apply(cleankrakenabun[ranks==refrank,], 2, function(x){ x / sum(x) })	
		#~ normcleankrakenabun = apply(cleankrakenabun[ranks=='S',], 2, function(x){ x / sum(x) })	
		#~ normcleankrakenabun = apply(cleankrakenabun[ranks=='G',], 2, function(x){ x / sum(x) })	
				
		################################
		## check sequencing depth effect

		logabun = lapply(individuals, function(samp){ a = log(normcleankrakenabun[,samp]) ; a[!is.infinite(a)] })
		names(logabun) = individuals
		df.logabun = do.call(rbind, lapply(individuals, function(samp){ data.frame(Depth=depths[samp, depthset], Sample=samp, Locality=as.character(criteria$localities[samp]), log.abun=logabun[[samp]]) }) )
		#~ print(anova(lm(log.abun ~ Depth * Sample * Locality, data=df.logabun)))
		print(anova(lm(log.abun ~ Depth, data=df.logabun)))
		for (thinthresh in -10:-5){
			print(sprintf("thinning threshold: 1e%d", thinthresh))
			print(sprintf("remain %d %s", length(which(df.logabun$log.abun > thinthresh)), refranks[refrank]))
			aovthin = anova(lm(log.abun ~ Depth, data=df.logabun[df.logabun$log.abun > thinthresh,]))
			if (aovthin["Depth", "Pr(>F)"] > 0.05){
				print(aovthin)
				break
			}
		}
			
		print(sprintf("-> sets thinning threshold: 1e%d", thinthresh))
	#~ 	thinthresh = -6

		logabunthin = lapply(logabun, function(x){ x[x > thinthresh] })

		pdf(file.path(outdir, sprintf('distrib_kraken_%s_abundances.pdf', refranks[refrank])), width=15, height=15)

		boxplot(logabun, names=individuals, col=coulloc[localities], las=2, ylab=sprintf('log(%s abundances)', refranks[refrank]))
		legend('topright', legend=names(coulloc), fill=coulloc, cex=2)

		n=1
		for (samp in individuals){
			if (n==1) plot(density(logabun[[samp]]), xlim=c(-20,0), ylim=c(0, 0.3), col=coulloc[localities[n]], xlab=sprintf('log(%s abundances)', refranks[refrank]))
			else lines(density(logabun[[samp]]), col=coulloc[localities[n]])
			n = n + 1
		}
		legend('topright', legend=names(coulloc), fill=coulloc, cex=2)
		text(-9, 0.15, labels=paste(c("log.abun ~ Locality", "", capture.output(print(anova(lm(log.abun ~ Locality, data=df.logabun))))), collapse='\n'), adj=0, cex=1.5)

		boxplot(logabun, names=individuals, col=couldepth[couldepvec], las=2, ylab=sprintf('log(%s abundances)', refranks[refrank]))
		legend('topright', legend=c(sprintf("log10(reads) = %.3g", min(depths[, depthset])), sprintf("log10(reads) = %.3g", max(depths[, depthset]))), fill=couldepth[c(1,20)], cex=2, bg='white')

		boxplot(logabunthin, names=individuals, col=couldepth[couldepvec], las=2, ylab=sprintf('log(%s abundances)', refranks[refrank]))
		legend('topright', legend=c(sprintf("log10(reads) = %.3g", min(depths[, depthset])), sprintf("log10(reads) = %.3g", max(depths[, depthset]))), fill=couldepth[c(1,20)], cex=2, bg='white')

		n=1
		for (samp in individuals){
			if (n==1) plot(density(logabun[[samp]]), xlim=c(-20,0), ylim=c(0, 0.3), col=couldepth[couldepvec[n]], lwd=2, xlab=sprintf('log(%s abundances)', refranks[refrank]))
			else lines(density(logabun[[samp]]), col=couldepth[couldepvec[n]], lwd=2)
			n = n + 1
		}
		legend('topright', legend=c(sprintf("log10(reads) = %.3g", min(depths[, depthset])), sprintf("log10(reads) = %.3g", max(depths[, depthset]))), fill=couldepth[c(1,20)], cex=2, bg='white')
		text(-9, 0.15, labels=paste(c("log.abun ~ Depth", "", capture.output(print(anova(lm(log.abun ~ Depth, data=df.logabun))))), collapse='\n'), adj=0, cex=1.5)

		n=1
		for (samp in individuals){
			if (n==1) plot(density(logabunthin[[samp]]), xlim=c(-9,0), ylim=c(0, 0.3), col=couldepth[couldepvec[n]], lwd=2, xlab=sprintf('log(%s abundances)', refranks[refrank]))
			else lines(density(logabunthin[[samp]]), col=couldepth[couldepvec[n]], lwd=2)
			n = n + 1
		}
		legend('topright', legend=c(sprintf("log10(reads) = %.3g", min(depths[, depthset])), sprintf("log10(reads) = %.3g", max(depths[, depthset]))), fill=couldepth[c(1,20)], cex=2, bg='white')
		text(-3.5, 0.15, labels=paste(c(sprintf("log.abun ~ Depth (data truncated at 10e%d)", thinthresh), "", capture.output(print(anova(lm(log.abun ~ Depth, data=df.logabun[df.logabun$log.abun > thinthresh,]))))), collapse='\n'), adj=0, cex=1.5)

		dev.off()

		##############################################
		## Differential abundances of all species


		# PCA on Kraken aundances
		#~ shitty = "run1_Cae45"
		#~ cleans = which(colnames(normcleankrakenabun)!=shitty)
		#~ normcleankrakenabunnoshitty = normcleankrakenabun[,cleans]
		#~ noshit.taxa = rownames(normcleankrakenabun)[apply(normcleankrakenabunnoshitty, 1, sum)>1] 
		#~ normcleankrakenabunnoshit = normcleankrakenabun[noshit.taxa,]
		#~ pca.abun = dudi.pca(t(normcleankrakenabunnoshit), scannf=F, nf=8)
		thinednormcleankrakenabun = normcleankrakenabun[apply(normcleankrakenabun, 1, function(x){ all(x > 10^thinthresh) }),]
		
		
		## scaling of variables is made implicitly in the LDA, so use scaling in the PCA on which it is plotted
		#~ pca.abun = dudi.pca(t(normcleankrakenabun), scannf=F, nf=8)
		pca.abun = dudi.pca(t(thinednormcleankrakenabun), scannf=F, nf=8)
	#~ 	pca.abun = dudi.pca(t(thinednormcleankrakenabun), scannf=F, nf=8, scale=FALSE)
		print(c("PCA Eigen values:", pca.abun$eig/sum(pca.abun$eig)))

		####################################
		# using Linear Discriminant Analysis

		# apply LDA to the dataset with HG vs. (Farmers or Controls) only
		critpairs = list(list(c('Palawan_Mount', 'Luzon_Mount'), c('Palawan_Mount', 'Luzon_Coast')), list(c('HG', 'AG'), c('HG', 'WC'))) ; names(critpairs) = c('localities', 'lifestyles')
		for (crit in names(critpairs)){
			for (critpair in critpairs[[crit]]){
				print(critpair)
				criterion = criteria[[crit]]
				crpr = criterion %in% critpair
				faccrpr = as.factor(as.character(criterion[crpr]))
				# subselects only samples included in compared sets
				crprthinednormcleankrakenabun = thinednormcleankrakenabun[,crpr]
				# find non constant variables
				nocons = sapply(1:dim(crprthinednormcleankrakenabun)[1], function(i){ var(crprthinednormcleankrakenabun[i,]) > 0 })

				## NOT a DAPC: discrimin uses a PCA just to recycle the underlying dudi object 
				## hence scaling or not loads in PCA is irrelevant... to the LDA, but not to the PCA on which it is plotted
				#~ disc.abun = lda(t(crprthinednormcleankrakenabun[nocons,]), grouping=faccrpr, tol=1.0e-8)
				disc.abun = discrimin(dudi.pca(t(crprthinednormcleankrakenabun[nocons,]), scannf=F, nf=length(which(nocons))), fac=faccrpr, scannf=F, nf=1)
				topdiscvar = disc.abun$va[order(abs(disc.abun$va), decreasing=T),]
	#~ 			disc.abun = discrimin(pca.abun, fac=criterion, scannf=F, nf=2)
	#~ 			topdiscvar = disc.abun$va[order(abs(disc.abun$va[,1]), decreasing=T),]
				topdiscspe20 = names(topdiscvar)[1:20]
				topdiscspe50 = names(topdiscvar)[1:50]
				
				# perform alll independent t-tests on species abundances
				ttestpvals = apply(thinednormcleankrakenabun[,crpr], 1, function(speabun){
	#~ 				tt = t.test(speabun ~ criterion[crpr])
					tt = wilcox.test(speabun ~ criterion[crpr])
					return(tt$p.value)
				})
				ranked.ttestpvals = ttestpvals[order(ttestpvals)]
				# Benjamini–Hochberg procedure
				signifthresh = 0.05
				m = length(ranked.ttestpvals)
				BHcor.ttestpvvals = sapply(1:m, function(k){ ranked.ttestpvals[k] * m / k })
				topsignif = which(BHcor.ttestpvvals <= signifthresh)
				topsignifdiscspe = names(topsignif)
				topsignifdisc = order(ttestpvals)[topsignif]
				
	#~ 			if (crit == 'lifestyles'){
					pdf(sprintf('%s/top_dirscrim_%s_diffabun_%svs%s_truncdata.pdf', outdir, refranks[refrank], levels(faccrpr)[1], levels(faccrpr)[2]), height=10, width=15)
					nr = 5
					layout(matrix(1:(5*nr), nr, 5, byrow=T))
					critorder = sapply(levels(criteria[[crit]]), function(x){ which(ordered.crits[[crit]]==x) })
					print(crit)
		#~ 			topsignifdiscspe = names(signifdiscvar80)[order(abs(signifdiscvar80), decreasing=T)][1:20]
					for (spe in topdiscspe50){
	#~ 				for (spe in topsignifdiscspe){
						Spe = gsub('([^s][^p])\\.', '\\1 ', spe)
						if (Spe %in% rownames(thinednormcleankrakenabun)){
	#~ 						print(Spe)
							signifmark = ifelse(Spe %in% topsignifdiscspe, '*', '')	
							boxplot(thinednormcleankrakenabun[Spe, ] ~ criteria[[crit]], horizontal=T, col=lcoul[[crit]][levels(criteria[[crit]])], las=2,
							 main=sprintf("%s\nLD loading: %.3g %s", Spe, disc.abun$va[spe,], as.character(signifmark)), at=critorder)
						}
					}
					dev.off()
					
					if (length(topsignifdiscspe) > 1){
						DAorderedtncnkrakenabun = thinednormcleankrakenabun[order(disc.abun$va, decreasing=T), ]
						crit.ave.log.abuns = apply(DAorderedtncnkrakenabun[rownames(DAorderedtncnkrakenabun) %in% topsignifdiscspe, ], 1, function(speabun){
							sapply(levels(criterion), function(c){ log10(mean(speabun[criterion==c])) })
						})
					
						pdf(sprintf('%s/top_signif_%s_diffabun_%svs%s_truncdata.pdf', outdir, refranks[refrank], levels(faccrpr)[1], levels(faccrpr)[2]), height=10, width=8)
						par(mar=c(4, 10, 4, 4))
						barplot(6+crit.ave.log.abuns[levels(criterion),], beside=T, horiz=T, las=2, col=lcoul[[crit]], axes=F)
						axis(side=1, at=0:5, labels=(0:5)-6)
						axis(side=3, at=0:5, labels=(0:5)-6)
						mtext(side=1, sprintf("log10(%s abundances)", refranks[refrank]), line=3)
						legend('topright', leg=levels(criterion), fill=lcoul[[crit]])
						dev.off()
					}
	#~ 			}

				
				##################################
				# PCA with max discriminant taxa:
				
				#~ pdf(sprintf('%s/PCA_kraken_with_top_dirscrim_%s_diffabun_%svs%s.pdf', outdir, refranks[refrank], levels(faccrpr)[1], levels(faccrpr)[2]), height=20, width=20)
				pdf(sprintf('%s/PCA_kraken_with_top_dirscrim_%s_diffabun_%svs%s_truncdata.pdf', outdir, refranks[refrank], levels(faccrpr)[1], levels(faccrpr)[2]), height=20, width=20)
				#~ layout(matrix(1:4, 2, 2, byrow=T))
				s.class(pca.abun$li, fac=criteria[[crit]], col=lcoul[[crit]])
				plotPCAabun(pca.abun, criteria[[crit]], lcoul[[crit]], subsample=NULL, np=2, topdisc=topdiscspe50, scalingvarvect=scalingvarvect)
				dev.off()
				#~ write.table(pca.abun$c1[topdiscspe50,], file=sprintf('%s/top_dirscrim_%s_diffabun_%svs%s.txt', outdir, refranks[refrank], levels(faccrpr)[1], levels(faccrpr)[2]))
				write.table(pca.abun$c1[topdiscspe50,], file=sprintf('%s/top_dirscrim_%s_diffabun_%svs%s_truncdata.txt', outdir, refranks[refrank], levels(faccrpr)[1], levels(faccrpr)[2]))
			}
		}
	}

}
