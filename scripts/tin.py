#!/usr/bin/env python3
'''-------------------------------------------------------------------------------------------------
This program calculates transcript integrity number (TIN) for each transcript (or gene) in
BED file. TIN is conceptually similar to RIN (RNA integrity number) but provides transcript
level measurement of RNA quality and is more sensitive to measure low quality RNA samples:

1) TIN score of a transcript is used to measure the RNA integrity of the transcript.
2) Median TIN score across all transcripts can be used to measure RNA integrity  of that
   "RNA sample".
3) TIN ranges from 0 (the worst) to 100 (the best). TIN = 60 means: 60% of the transcript
   has been covered if the reads coverage were uniform.
4) TIN will be assigned to 0 if the transcript has no coverage or covered reads is fewer than
   cutoff.
-------------------------------------------------------------------------------------------------'''
import sys,os
import math,random
from bisect import bisect_right
from optparse import OptionParser
from qcmodule import getBamFiles
from qcmodule import BED
from numpy import mean,median,std
from time import strftime
import pysam
from bx.intervals import *
from mpire import WorkerPool

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="5.0.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"


def printlog (mesg):
	'''
	print mesg into stderr with time string appending to it.
	'''
	mesg="@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	print(mesg, file=sys.stderr)

def uniqify(seq):
	'''
	duplicated members only keep one copy. [1,2,2,3,3,4] => [1,2,3,4].
	'''
	seen = set()
	return [x for x in seq if x not in seen and not seen.add(x)]

def shannon_entropy(arg):
	'''
	calculate shannon's H = -sum(P*log(P)). arg is a list of float numbers. Note we used
	natural log here.
	'''
	lst_sum = sum(arg)
	entropy = 0.0
	for i in arg:
		entropy += (i/lst_sum) * math.log(i/lst_sum)
	if entropy == 0:
		return 0
	else:
		return -entropy

def build_bitsets(list):
	'''
	build intevalTree from list
	'''
	ranges={}
	for l in list:
		chrom =l[0]
		st = l[1]
		end = l[2]
		if chrom not in ranges:
			ranges[chrom] = Intersecter()
		ranges[chrom].add_interval( Interval( st, end ) )
	return ranges

def union_exons(refbed):
	'''
	take the union of all exons defined in refbed file and build bitset
	'''
	tmp = BED.ParseBED(refbed)
	all_exons = tmp.getExon()
	unioned_exons = BED.unionBed3(all_exons)
	exon_ranges = build_bitsets(unioned_exons)
	return exon_ranges

def union_exons_by_chrom(refbed):
	'''
	take union of exons and return plain intervals grouped by chromosome
	'''
	tmp = BED.ParseBED(refbed)
	all_exons = tmp.getExon()
	unioned_exons = BED.unionBed3(all_exons)
	ret = {}
	for chrom, st, end in unioned_exons:
		ret.setdefault(chrom, []).append((chrom, st, end))
	return ret

def estimate_bg_noise(chrom, tx_st, tx_end, samfile, e_ranges):
	'''
	estimate background noise level for a particular transcript
	'''
	intron_sig = 0.0	# reads_num * reads_len
	alignedReads = samfile.fetch(chrom,tx_st,tx_end)
	for aligned_read in alignedReads:
		if aligned_read.is_qcfail:continue 
		if aligned_read.is_unmapped:continue
		if aligned_read.is_secondary:continue
		read_start = aligned_read.pos
		if read_start < tx_st: continue
		if read_start >= tx_end: continue
		read_len = aligned_read.qlen
		if len(e_ranges[chrom].find(read_start, read_start + read_len)) > 0:
			continue
		intron_sig += read_len
	return intron_sig

def genomic_positions(refbed, sample_size):
	''' 
	return genomic positions of each nucleotide in mRNA. sample_size: Number of nucleotide
	positions sampled from mRNA.
	'''
	if refbed is None:
		print("You must specify a bed file representing gene model\n", file=sys.stderr)
		exit(0)
	
	for line in open(refbed):
		try:
			if line.startswith(('#','track','browser')):continue  
			# Parse fields from gene tabls
			fields = line.split()
			chrom     = fields[0]
			tx_start  = int( fields[1] )
			tx_end    = int( fields[2] )
			geneName      = fields[3]
			strand    = fields[5]
			cdsStart = int(fields[6]) + 1	#convert to 1-based
			cdsEnd = int(fields[7])
			exon_count = int(fields[9])
			mRNA_size = sum([int(i) for i in fields[10].strip(',').split(',')])
			geneID = '_'.join([str(j) for j in (chrom, tx_start, tx_end, geneName, strand)])
				
			exon_starts = [int(i) for i in fields[11].rstrip(',\n').split(',')]
			exon_starts = [x + tx_start for x in exon_starts]
			exon_sizes = [int(i) for i in fields[10].rstrip(',\n').split(',')]
			exon_ends = [x + y for x, y in zip(exon_starts, exon_sizes)]
			intron_size = tx_end - tx_start - mRNA_size
			if intron_size < 0:
				intron_size = 0
		except:
			print("[NOTE:input bed must be 12-column] skipped this line: " + line, end='', file=sys.stderr)
			continue
		
		chose_bases=[tx_start+1, tx_end]
		exon_bounds = []
		if mRNA_size <= sample_size:	# return all bases of mRNA
			for st,end in zip(exon_starts,exon_ends):
				chose_bases.extend(range(st+1,end+1))		#1-based coordinates on genome, include exon boundaries
			yield (geneName, chrom, tx_start, tx_end, intron_size, chose_bases)
		elif mRNA_size > sample_size:
			step_size = int(mRNA_size/sample_size)
			exon_cumsizes = []
			running = 0
			for st,end in zip(exon_starts,exon_ends):
				exon_bounds.append(st+1)
				exon_bounds.append(end)
				running += (end - st)
				exon_cumsizes.append(running)
			indx = range(0, mRNA_size, step_size)
			chose_bases = []
			for i in indx:
				exon_idx = bisect_right(exon_cumsizes, i)
				prev_size = 0
				if exon_idx > 0:
					prev_size = exon_cumsizes[exon_idx - 1]
				chose_bases.append(exon_starts[exon_idx] + 1 + (i - prev_size))
			yield (geneName, chrom, tx_start, tx_end, intron_size, uniqify(exon_bounds + chose_bases))
		
def transcript_read_metrics(samfile, chrom, tx_st, tx_end, cutoff, subtract_bg=False, e_ranges=None):
	'''
	Single-pass read scan for a transcript region.
	Returns (has_min_cov, intron_signal).
	'''
	has_min_cov = False
	intron_sig = 0.0
	read_count = set()
	try:
		alignedReads = samfile.fetch(chrom, tx_st, tx_end)
		for aligned_read in alignedReads:
			if aligned_read.is_qcfail:continue
			if aligned_read.is_unmapped:continue
			if aligned_read.is_secondary:continue
			read_start = aligned_read.pos
			if read_start < tx_st: continue
			if read_start >= tx_end: continue
			read_count.add(read_start)
			if len(read_count) > cutoff:
				has_min_cov = True
				if not subtract_bg:
					break
			if subtract_bg:
				read_len = aligned_read.qlen
				if len(e_ranges[chrom].find(read_start, read_start + read_len)) > 0:
					continue
				intron_sig += read_len
		return has_min_cov, intron_sig
	except:
		return False, 0.0

	
def genebody_coverage(samfile, chrom, positions, bg_level = 0):
	'''
	calculate coverage for each nucleotide in *positions*. some times len(cvg) < len(positions)
	because positions where there is no mapped reads were ignored.
	'''
	cvg = []
	pos_set = set(positions)
	start = positions[0] - 1
	end = positions[-1]
	
	try:
		for pileupcolumn in samfile.pileup(chrom, start, end, truncate=True):
			ref_pos = pileupcolumn.pos+1
			if ref_pos not in pos_set: continue
			if pileupcolumn.n == 0:
				cvg.append(0.0)
				continue
			cover_read = 0.0
			for pileupread in pileupcolumn.pileups:
				if pileupread.is_del: continue
				if pileupread.alignment.is_qcfail:continue 
				if pileupread.alignment.is_secondary:continue 
				if pileupread.alignment.is_unmapped:continue
				#if pileupread.alignment.is_duplicate:continue
				cover_read +=1.0
			cvg.append(cover_read)
	except:
		cvg = []
	
	if bg_level <= 0:
		return cvg
	else:
		tmp = []
		for i in cvg:
			subtracted_sig = int(i - bg_level)
			if subtracted_sig > 0:
				tmp.append(subtracted_sig)
			else:
				tmp.append(0)
		return tmp


def tin_score(cvg, l):
	'''calcualte TIN score'''
	
	if len(cvg) == 0:
		tin = 0
		return tin
	
	cvg_eff = [float(i) for i in cvg if float(i) > 0]	#remove positions with 0 read coverage
	entropy = shannon_entropy(cvg_eff)
	
	tin = 100*(math.exp(entropy)) / l
	return tin


def process_chromosome_job(job):
	'''
	Worker job for one chromosome in one BAM file.
	'''
	bam_file, tx_records, minimum_coverage, subtract_bg, exon_intervals = job
	samfile = pysam.Samfile(bam_file, "rb")
	rows = []
	tins = []

	exon_ranges = None
	if subtract_bg:
		exon_ranges = build_bitsets(exon_intervals)
	for idx, gname, i_chr, i_tx_start, i_tx_end, intron_size, pick_positions in tx_records:
		noise_level = 0.0
		has_min_cov, intron_signals = transcript_read_metrics(
			samfile,
			i_chr,
			i_tx_start,
			i_tx_end,
			minimum_coverage,
			subtract_bg,
			exon_ranges
		)
		if has_min_cov is not True:
			rows.append((idx, gname, i_chr, i_tx_start, i_tx_end, 0.0))
			continue
		if subtract_bg:
			if intron_size > 0:
				noise_level = intron_signals/intron_size
		coverage = genebody_coverage(samfile, i_chr, pick_positions, noise_level)
		tin1 = tin_score(cvg = coverage, l = len(pick_positions))
		tins.append(tin1)
		rows.append((idx, gname, i_chr, i_tx_start, i_tx_end, tin1))
	samfile.close()
	return rows, tins


def main():
	usage="%prog [options]" + '\n' + __doc__ + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_files",
				   help='Input BAM file(s). "-i" takes these input: ' \
				   '1) a single BAM file. ' \
				   '2) "," separated BAM files (no spaces allowed). ' \
				   '3) directory containing one or more bam files. ' \
				   '4) plain text file containing the path of one or more bam files '
				   '(Each row is a BAM file path). All BAM files should be sorted and' \
				   'indexed using samtools. [required]')
	parser.add_option("-r","--refgene",action="store",type="string",dest="ref_gene_model",help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]")
	parser.add_option("-c","--minCov",action="store",type="int",dest="minimum_coverage",default=10,help="Minimum number of read mapped to a transcript. default=%default")
	parser.add_option("-n","--sample-size",action="store",type="int",dest="sample_size",default=100,help="Number of equal-spaced nucleotide positions picked from mRNA. Note: if this number is larger than the length of mRNA (L), it will be halved until it's smaller than L. default=%default")
	parser.add_option("-p","--processes",action="store",type="int",dest="processes",default=1,help="Number of worker processes for chromosome-level parallel mode (MPIRE). default=%default")
	parser.add_option("-s","--subtract-background",action="store_true",dest="subtract_bg",help="Subtract background noise (estimated from intronic reads). Only use this option if there are substantial intronic reads.")
	(options,args)=parser.parse_args()

	if options.processes < 1:
		print("'-p/--processes' must be >= 1", file=sys.stderr)
		sys.exit(0)
	
	exon_intervals_by_chrom = {}
	# if '-s' was set
	if options.subtract_bg:
		exon_ranges = union_exons(options.ref_gene_model)
		exon_intervals_by_chrom = union_exons_by_chrom(options.ref_gene_model)
		
	if options.sample_size < 0:
		print("Number of nucleotide can't be negative", file=sys.stderr)
		sys.exit(0)
	elif options.sample_size >1000:
		print("Warning: '-n' is too large! Please try smaller '-n' valeu if program is running slow.", file=sys.stderr)
		
	if not (options.input_files and options.ref_gene_model):
		parser.print_help()
		sys.exit(0)

	if not os.path.exists(options.ref_gene_model):
		print('\n\n' + options.ref_gene_model + " does NOT exists" + '\n', file=sys.stderr)
		parser.print_help()
		sys.exit(0)

	transcript_records = []
	for rec in genomic_positions(refbed = options.ref_gene_model, sample_size = options.sample_size):
		gname, i_chr, i_tx_start, i_tx_end, intron_size, pick_positions = rec
		transcript_records.append((gname, i_chr, i_tx_start, i_tx_end, intron_size, sorted(pick_positions)))
	indexed_records = [(idx, rec[0], rec[1], rec[2], rec[3], rec[4], rec[5]) for idx, rec in enumerate(transcript_records)]
	chrom_groups = {}
	for rec in indexed_records:
		chrom_groups.setdefault(rec[2], []).append(rec)
	chrom_keys = sorted(chrom_groups.keys())
		
	printlog("Get BAM file(s) ...")
	bamfiles = sorted(getBamFiles.get_bam_files(options.input_files))
	
	if len(bamfiles) <= 0:
		print("No BAM file found, exit.", file=sys.stderr)
		sys.exit(0)
	else:
		print("Total %d BAM file(s):" % len(bamfiles), file=sys.stderr)
		for f in bamfiles:
			print("\t" + f, file=sys.stderr)	
	
	
	for f in bamfiles:
		printlog("Processing " + f)
		
		SUM = open(os.path.basename(f).replace('bam','') + 'summary.txt','w')
		print("\t".join(['Bam_file','TIN(mean)', 'TIN(median)','TIN(stdev)']), file=SUM)
		
		OUT = open(os.path.basename(f).replace('bam','') + 'tin.xls','w')
		print("\t".join(["geneID","chrom", "tx_start", "tx_end","TIN"]), file=OUT)

		sample_TINs = []	#sample level TIN, values are from different genes
		finish = 0

		use_parallel = (options.processes > 1 and len(chrom_keys) > 1)

		if use_parallel:
			printlog("Chromosome-parallel mode enabled with %d worker(s)" % options.processes)
			jobs = [(f, chrom_groups[k], options.minimum_coverage, options.subtract_bg, exon_intervals_by_chrom.get(k, [])) for k in chrom_keys]
			with WorkerPool(n_jobs=options.processes) as pool:
				results = pool.map(process_chromosome_job, jobs)
			all_rows = []
			for rows, tins in results:
				all_rows.extend(rows)
				sample_TINs.extend(tins)
			all_rows.sort(key=lambda x: x[0])
			for _, gname, i_chr, i_tx_start, i_tx_end, tin1 in all_rows:
				finish += 1
				print('\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, tin1)]), file=OUT)
				print(" %d transcripts finished\r" % (finish), end='', file=sys.stderr)
		else:
			samfile = pysam.Samfile(f, "rb")
			for gname, i_chr, i_tx_start, i_tx_end, intron_size, pick_positions in transcript_records:	
				finish += 1
				noise_level = 0.0
				
				has_min_cov, intron_signals = transcript_read_metrics(
					samfile,
					i_chr,
					i_tx_start,
					i_tx_end,
					options.minimum_coverage,
					options.subtract_bg,
					exon_ranges
				)
				if has_min_cov is not True:
					print('\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, 0.0)]), file=OUT)
					continue
					
				if options.subtract_bg:
					if intron_size > 0:
						noise_level = intron_signals/intron_size					

				coverage = genebody_coverage(samfile, i_chr, pick_positions, noise_level)
				
				#for a,b in zip(sorted(pick_positions),coverage):
				#	print str(a) + '\t' + str(b)
				
				tin1 = tin_score(cvg = coverage, l = len(pick_positions))
				sample_TINs.append(tin1)
				print('\t'.join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, tin1)]), file=OUT)
				print(" %d transcripts finished\r" % (finish), end='', file=sys.stderr)
			samfile.close()
		
		print("\t".join( [str(i) for i in (os.path.basename(f), mean(sample_TINs), median(sample_TINs), std(sample_TINs))]), file=SUM)
		OUT.close()
		SUM.close()

if __name__ == '__main__':
	
	main()

	
