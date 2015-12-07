#!/usr/bin/python
# -*- coding: utf8 -*-

import os
import gzip
import copy
import getopt, sys
import cPickle
import time
import glob

from parseKrona import *

def parsereadheader(line, matepairtagseps=["#", " ", "\\"], stripchars="@", inpoolpair=False, SRAtag=False):
	sc = stripchars + '\\\n\\n'
	if "'" in line: li = line.split("'")[-1]	
	else: li = line
	if inpoolpair: lmptagsep = []	# provide no separator for splitting the mate read tag, which will be kept
	else: lmptagsep = matepairtagseps
	if li==")":
		return None
	elif SRAtag:
		# header starts with @SRACCESSION.READNUM PROPER:HEADER:TAG, rely on the SRACCESSION.READNUM part
		return li.strip(sc).split(" ")[0]
	else:
		#~ print li
		for mptagsep in lmptagsep:
			# tries all potential split character in the set (by default '#' and ' ') as separator for removing the tag indicating the mate pair reads 
			if mptagsep in li:
				#~ print 'here', ':'.join(li.strip(sc).split(mptagsep)[0].split(':')[2:])
				return ':'.join(li.strip(sc).split(mptagsep)[0].split(':')[2:])
		else:
			# no separator matched, consider there was no mate pair read tag
			#~ print 'there', ':'.join(li.strip(sc).split(':')[2:])
			return ':'.join(li.strip(sc).split(':')[2:])
		
def openFastQCon(lnffastq, fastqdir, gzipped='guess'):
	dnffastqcon = {}
	for nffastq in lnffastq:
		pathffastq = "%s/%s"%(fastqdir, nffastq)
		# determines if file is gzipped or not and use the corresponding open()
		if gzipped==True or (gzipped=='guess' and nffastq.split('.')[-1]=="gz"):
			ffastq = gzip.open(pathffastq, 'r')
		else:
			ffastq = open(pathffastq, 'r')
		dnffastqcon[nffastq] = ffastq
	return dnffastqcon

def parseFastQtoOffsets(lnffastq, fastqdir, gzipped='guess', verbose=False, informat='FASTQ', inpoolpair=False, nreadtell=10000, SRAtag=False):
	"""parsing FASTQ files to index them with read headers"""
	# select sequence format to be read
	if informat=='FASTQ':
		fastqfmt = True
		nlineread = 4
		headerstartchar = '@'
	elif informat=='FASTA':
		fastqfmt = False
		nlineread = 2
		headerstartchar = '>'
	else:
		raise ValueError, "inappropriate format was specified, must be either 'FASTQ' or 'FASTA'"
	dnfindex = {}
	dnffastqcon = {}
	for nffastq in lnffastq:
		if verbose: print "parsing", nffastq
		dhead2offset = {} # dict of fastq header to offset ((start, length) position in file)
		dnffastq1con = openFastQCon([nffastq], fastqdir, gzipped=gzipped)
		ffastq = dnffastq1con[nffastq]
		nline = 0
		nreads = 0
		start = ffastq.tell()
		prevend = None
		prevstart = 0
		readstart = None
		readtag = None
		line = ffastq.readline()
		while line:
			if verbose and (nreads%nreadtell)==0 and (nline%nlineread): print "parsed %d reads"%nreads
			if (fastqfmt and nline%nlineread==0) or ((not fastqfmt) and line.startswith(headerstartchar)):
				# FASTQ header found every nlinereadth line ; FASTA header frequency is not previsible from generic fasat specifications
				nreads += 1
				if prevend!=None and nline>0:
					# store previous read information
					if readtag not in dhead2offset: dhead2offset[readtag] = (readstart, prevend-readstart)
					else: raise IndexError, "recorded several reads with header %s"%(str(readtag))
					if verbose and nreads%(nreadtell)==0: print readtag, (readstart, prevend-readstart)
				# start parsing current read
				readtag = parsereadheader(line, SRAtag=SRAtag, stripchars=headerstartchar, inpoolpair=inpoolpair)
				if str(readtag)=="": raise ValueError, "unproper header format in:\n%s"%line
				readstart = prevstart
			prevend = start
			prevstart = start
			line = ffastq.readline()
			nline += 1
			start = ffastq.tell()
		else:
			if readtag not in dhead2offset: dhead2offset[readtag] = (readstart, prevend-readstart)
			else: raise IndexError, "recorded several reads with header %s"%(str(readtag))
			if verbose: 
				print "parsed %d reads; finished"%nreads
				print readtag, (readstart, prevend-readstart)
		dnffastqcon[nffastq] = ffastq
		dnfindex[nffastq] = dhead2offset
		
	return dnfindex, dnffastqcon

def combineOffsets(loffsets):
	if len(loffsets)<=1: return loffsets
	loffsets.sort()
	lout = []
	laststart = None
	lastend = None
	for offset in loffsets:
		if not lastend:
			laststart = offset[0]
			lastend = sum(offset)
		else:
			if lastend==offset[0]:
				lastend += offset[1]
			else:
				lout.append((laststart, lastend-laststart))
				laststart = offset[0]
				lastend = sum(offset)
	else:
		lout.append((laststart, lastend-laststart))
	return lout

def fetch_reads(taxonomy, nameortaxid, dirkronamemberfiles, dnfindex, dnffastqcon, dirout, excludenameortaxids=[], informat='FASTQ', inpoolpair=False, outtag=None, cumulative=True, concat=False, outmode='a', gzipped='guess', fastIO=False, verbose=False, SRAtag=False):
	# select sequence format to be read
	if informat=='FASTQ':
		headerstartchar = '@'
	elif informat=='FASTA':
		headerstartchar = '>'
	else:
		raise ValueError, "inappropriate format was specified, must be either 'FASTQ' or 'FASTA'"
	# prepare output in separate files
	if outtag: tag = str(outtag)
	else: tag = str(nameortaxid)
	dnfout = {}
	for nffastq in dnfindex:
		if gzipped==True or (gzipped=='guess' and nffastq.split('.')[-1]=="gz"):
			fout = gzip.open("%s/%s.%s"%(dirout, tag, nffastq), outmode+'b')
		else:
			fout = open("%s/%s.%s"%(dirout, tag, nffastq), outmode)
		dnfout[nffastq] = fout
	taxnode = taxonomy.get(nameortaxid, taxonomy.taxidgetnode(nameortaxid)) # first consider the input as the name, then as a numerical taxid
	excltaxnodes = [taxonomy.get(exclnameortaxid, taxonomy.taxidgetnode(exclnameortaxid)) for exclnameortaxid in excludenameortaxids]
	if not cumulative:
		lmemb = [taxnode.members]
	else:
		lmemb = taxnode.get_cumulative_members(exclude=excltaxnodes)
	if not fastIO:
		for memb in lmemb:
			fmemb = open("%s/%s"%(dirkronamemberfiles, memb), 'r')
			if verbose: print "parsing", memb
			n = 0
			for line in fmemb:
				readtag = parsereadheader(line, SRAtag=SRAtag, stripchars=headerstartchar, inpoolpair=inpoolpair)
				for nffastq in dnfindex:
					readoffset = dnfindex[nffastq].get(readtag)
					if readoffset:
						n += 1
						ffastq = dnffastqcon[nffastq]
						ffastq.seek(readoffset[0])
						fastqread = ffastq.read(readoffset[1])
						fout = dnfout[nffastq]
						fout.write(fastqread)
			if verbose: print "fetched %d reads from %d FastQ files"%(n, len(dnfindex))
			fmemb.close()
	else:
		if verbose: print "try combining I/O operations on read files"
		# try to accelerate processing by combining read and write operations of FASTQ reads that were contiguous in the input file
		lreadtags = []
		for memb in lmemb:
			if verbose: print "parsing", memb
			fmemb = open("%s/%s"%(dirkronamemberfiles, memb), 'r')
			lreadtags += [parsereadheader(line, SRAtag=SRAtag, stripchars=headerstartchar, inpoolpair=inpoolpair) for line in fmemb]
			if verbose: print "reads to fetch: [%s, ...]"%(str(lreadtags[0:3]).strip('[]'))
			fmemb.close()
		for nffastq in dnfindex:
			dtagoffsets = dnfindex[nffastq]
			if verbose: 
				print "searching matching reads from %s (%d recorded reads)"%(nffastq, len(dtagoffsets))
				print "{%s, ...}"%(str(dict(zip(dtagoffsets.keys()[0:3], dtagoffsets.values()[0:3]))).strip('{}'))
			#~ loffsets = [dtagoffsets[readtag] for readtag in lreadtags if (readtag in dtagoffsets)]
			#~ lfoundreadtags = list(set(lreadtags) & set(dtagoffsets.keys()))
			lfoundreadtags = lreadtags
			loffsets = [dtagoffsets[readtag] for readtag in lfoundreadtags if readtag] # why readtag should be null ?
			# equivalent to:
			# loffsets = []
			# for readtag in lfoundreadtags:
			# 	loffsets.append(dtagoffsets[readtag])
			if verbose: print "%d offsets to combine"%(len(loffsets))
			combloffsets = combineOffsets(loffsets)
			if verbose: print "combined into %d offsets"%(len(combloffsets))
			ffastq = dnffastqcon[nffastq]
			fout = dnfout[nffastq]
			for readoffset in combloffsets:
				ffastq.seek(readoffset[0])
				fastqread = ffastq.read(readoffset[1])
				fout.write(fastqread)	
			if verbose: print "fetched %d reads from FastQ file '%s'"%(len(loffsets), nffastq)
	for nffastq in dnfout:
		dnfout[nffastq].close()
		
def telltime(starttime, task='total'):
	endtime = time.time()
	tottime = endtime-starttime
	print "Spent %s (%fs) in %s"%(time.strftime('%Hh %Mm %Ss', time.gmtime(tottime)), tottime, task)
	return endtime

def main(nfhtml, dirhtmlfiles, fastqdir, outputdir, lftaxa, letaxa=[], informat='FASTQ', inpoolpair=False, SRAtag=False, cumulative=True, concat=False, fastIO=False, pickled=None, outmode='a', verbose=False):
	
	starttime = time.time()
	currtime = starttime
	
	# create output directory if does not exist
	if not os.path.exists(outputdir):
		os.mkdir(outputdir)
	
	# parsing HTML Krona files
	print "Parsing taxonomy in Krona HTML result file: '%s'"%(nfhtml)
	fhtml = open(nfhtml, 'r')
	kronaparser = KronaHTMLParser() #verbose=verbose
	kronaparser.feed(fhtml.read())
	fhtml.close()
	taxonomy = kronaparser.taxonomy[0]
	if not taxonomy: raise IndexError
	#~ else: print taxonomy
	# check that indicated taxa are present in the taxonomic tree
	for nameortaxid in lftaxa+letaxa:
		taxnode = taxonomy.get(nameortaxid, taxonomy.taxidgetnode(nameortaxid))
		if not taxnode: raise IndexError, "'%s' is not present as a taxon name or taxid in the reference taxonomy"%str(nameortaxid)
		print taxnode.members
		if verbose:
			lfoundmemb = list(set(os.listdir(dirhtmlfiles)) & set(taxnode.get_cumulative_members()))
			print "found %d member files matching taxon '%s': [%s, ...]"%(len(lfoundmemb), nameortaxid, str(lfoundmemb[0:5]).strip('[]'))
	else:
		if verbose: print "all taxa listed for read fetching %s or excluding %s are present in the reference taxonomy"%(str(lftaxa), str(letaxa))
	# parse reads
	if os.path.isdir(fastqdir):
		lnffastq = os.listdir(fastqdir)
		dirfastq = fastqdir
	else:
		dirfastq, prefixfastq = fastqdir.rsplit('/', 1)
		if '*' in fastqdir:
			lnffastq = [os.path.basename(nf) for nf in glob.glob(fastqdir)]
		else:
			lnffastq = [nf for nf in os.listdir(dirfastq) if nf.startswith(prefixfastq)]
		
	print "Parsing reads in %s files: %s"%(informat, str(lnffastq))
	if pickled:
		if os.path.exists(pickled):
			dnffastqcon = openFastQCon(lnffastq, dirfastq, gzipped='guess')
			if verbose: print "load pickled data from '%s'"%pickled
			fpickle = open(pickled, 'r')
			dnfindex = cPickle.load(fpickle)
			fpickle.close()
		else:
			dnfindex, dnffastqcon = parseFastQtoOffsets(lnffastq, dirfastq, verbose=verbose, SRAtag=SRAtag, informat=informat, inpoolpair=inpoolpair)
			currtime = telltime(starttime, task='parsing read offsets')	
			if verbose: print "dump pickled data to '%s'"%pickled
			fpickle = open(pickled, 'w')
			cPickle.dump(dnfindex, fpickle, protocol=2)
			fpickle.close()
	else:
		dnfindex, dnffastqcon = parseFastQtoOffsets(lnffastq, dirfastq, verbose=verbose, SRAtag=SRAtag, informat=informat, inpoolpair=inpoolpair)	
	# fetch reads for desired taxon
	print "Will extract reads form taxa ('%s'), excluding those from taxa ('%s')"%("', '".join([str(taxon) for taxon in lftaxa]), "', '".join([str(taxon) for taxon in letaxa]))
	for i in range(len(lftaxa)):
		taxon = lftaxa[i]
		if concat and len(lftaxa)>1:
			outtag = 'all'	# generic output file prefix
			if i == 0: mode = outmode
			else: mode = 'a'	# append to file after first iteration
		else:
			mode = outmode
			outtag = taxon
		fetch_reads(taxonomy, taxon, dirhtmlfiles, dnfindex, dnffastqcon, outputdir, excludenameortaxids=letaxa, informat=informat, inpoolpair=inpoolpair, outtag=outtag, cumulative=cumulative, concat=concat, fastIO=fastIO, outmode=mode, verbose=verbose, SRAtag=SRAtag)
		currtime = telltime(currtime, task='fetching and writing reads for %s'%taxon)	

	telltime(starttime, task='total')	

def usage():
	s =  "python parseKronaGetReadsByTaxid.py [options]\n"
	s += "\textracts subset of metagenomic reads by taxonomic classification using results of Kraken classifier\n"
	s += "Mandatory Options:\n"
	s += "\t--krona.html=file.path\n\t\tpath to Krona vizualisation HTML file\n"
	s += "\t--dir.fastq=directory.path\n\t\tpath to directory containing FASTQ files for metagenomic reads.\nA restrictive list of file can be obtained by appending the name, a degennerate name or a prefix of files to read to the directory path. A degenerate name matching is signalled by a '*' wildcard as in UNIX, but in this case option value must be given as --dir.fastq=value (as opposed to --dir.fastq value) to override shell expansion.\n"
	s += "\t--dir.output=directory.path\n\t\tpath to directory where output FASTQ files will be written\n"
	s += "\t--fetch.taxa=taxon1[,taxon2 [, ...]]\n\t\tcomma-separated list of names (enclosed in quotes) or NCBI taxids of taxa for which classified reads must be extracted\n"
	s += "\t--exclude.taxa=taxon1[,taxon2 [, ...]]\n\t\tcomma-separated list of names (enclosed in quotes) or NCBI taxids of taxa for which classified reads must NOT be extracted\n"
	s += "\tA value for at least one option among --fetch.taxa and --exclude.taxa must be specified. If no taxa to be fetched have been specified by --fetch.taxa, all reads but the ones specified by --exclude.taxa will be extracted.\n"
	s += "Facultative Options:\n"
	s += "\t--dir.krona.html.files=directory.path\n\t\tpath to directory containing supporting files for Krona vizualisation HTML file ('members' files). Defaults to the the path to Krona HTML file followed by '.files' suffix.\n"
	s += "\t--input.fasta\n\t\if input sequence reads are provided in simple FASTA format (rather than the default, FASTQ)\n"
	s += "\t--input.pooled.paired\n\t\if input sequence reads are provided in a single file (rather than two separate, the default), does not trim the paired-end #1 or #2 read tag\n"
	s += "\t-c --cumulative\n\t\twhether reads assigned to lower taxa must be cumulatively extracted along with those exactly assigned to the querried taxon (defaults to False)\n"
	s += "\t-k --concatenate\n\t\twhether to pool reads extracted for all querry taxa\n"
	s += "\t-a --append\n\t\twhether ouptut must be appended to existing files (defaults to False)\n"
	s += "\t-f --fastIO\n\t\twether to use memory caching of read offsets (use more memory but might accelerate I/O operations on FASTQ files when a majority of reads are to be collected)\n"
	s += "\t-p --pickled=file.path\n\t\twether to dump pickled the offsets of reads in FASTQ files for faster running next time or loading the data if they already exist\n"
	s += "\t-s --sra\n\t\tconsider the input read to have the SRA (Sequence Read Archive) FASTQ header format like '@SRRACCESSION.READNUM PROPER:HEADER:TAG'\n"
	s += "\t-v --verbose\n\t\tprint information on processing of reads and taxonomic data (defaults to False)\n"
	s += "\t-h --help\n\t\tprint this help message and quit\n"
	# think of options for rank-specific extraction of reads
	return s
	
	
if __name__=="__main__":

	#~ homedir = os.environ['HOME']
	#~ dirresult = '%s/oral_metagenomes'%(homedir)
	#~ outputdir = '%s/taxonomy-filtered_reads'%(dirresult)
	#~ nfhtml = '%s/run1_AE08.krona.filtered.050.html'%(dirresult)
	#~ dirhtmlf = '%s/run1_AE08.krona.filtered.050.html.files'%(dirresult)
	#~ lnfnodes = os.listdir(dirhtmlf)
	#~ fastqdir = '%s/fastq'%(dirresult)
	
	incall = " ".join(sys.argv)

	try:
		options, args = getopt.getopt(sys.argv[1:], 'ckvhafsp:', \
		["krona.html=", "dir.krona.html.files=", "dir.fastq=", "dir.output=", "fetch.taxa=", "exclude.taxa=", "input.fasta", "input.pooled.paired", "cumulative", "concatenate", "pickled=", "fastIO", "sra", "verbose", "append", "help"])
	except getopt.GetoptError, err:
		print err
		print usage()
		sys.exit(2)	
	dopt = dict(options)
	#~ print dopt
	
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(1)
	verbose = (('-v' in dopt) or ('--verbose' in dopt))

	nfhtml = dopt['--krona.html']
	dirhtmlfiles = dopt.get('--dir.krona.html.files', "%s.files"%(nfhtml))
	fastqdir = dopt['--dir.fastq']
	outputdir = dopt['--dir.output']
	slftaxa = dopt.get('--fetch.taxa')
	sletaxa = dopt.get('--exclude.taxa')
	SRAtag = (('-s' in dopt) or ('--sra' in dopt))
	pickled = dopt.get('-p', dopt.get('--pickled'))
	cumulative = (('-c' in dopt) or ('--cumulative' in dopt))
	concat = (('-k' in dopt) or ('--concatenate' in dopt))
	fastIO = (('-f' in dopt) or ('--fastIO' in dopt))
	outmode = 'a' if (('-a' in dopt) or ('--append' in dopt)) else 'w'
	informat = 'FASTA' if ('--input.fasta' in dopt) else 'FASTQ'
	inpoolpair = ('--input.pooled.paired' in dopt)
	
	if verbose: print "Command call was:\n%s"%incall	
	
	if slftaxa: lftaxa = [taxon.strip('"\'')for taxon in slftaxa.split(',')]
	elif sletaxa: lftaxa = ["Root"]
	else: print "Must specify a value for at least one option among --fetch.taxa and --exclude.taxa." ; sys.exit(2)
	if sletaxa: letaxa = [taxon.strip('"\'')for taxon in sletaxa.split(',')]
	else: letaxa = []
	
	main(nfhtml, dirhtmlfiles, fastqdir, outputdir, lftaxa=lftaxa, letaxa=letaxa, informat=informat, inpoolpair=inpoolpair, SRAtag=SRAtag, cumulative=cumulative, concat=concat, fastIO=fastIO, outmode=outmode, pickled=pickled, verbose=verbose)
