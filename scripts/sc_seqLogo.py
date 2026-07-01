#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program generates a DNA sequence logo from fasta or fastq format file.
It is useful to visualize the nucleotide compositions of sample barcodes,
cell barcodes and molecular barcodes.
"""

import sys
from qcmodule import fastq
from optparse import OptionParser
import logging

__author__ = "Liguo Wang"
__contributor__="Liguo Wang"
__copyright__ = "Copyright 2020, Mayo Clinic"
__credits__ = []
__license__ = "GPL"
__version__="5.0.4"
__maintainer__ = "Liguo Wang"
__email__ = "wangliguo78@gmail.com"
__status__ = "Development"


def main():
	usage= __doc__
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile",action="store",type="string", dest="in_file",
help="Input DNA sequence file in FASTQ (https://en.wikipedia.org/wiki/FASTQ_format#),\
FASTA (https://en.wikipedia.org/wiki/FASTA_format) or pure sequence format. All sequences must be\
 the same length. This file can be plain text or compressed format (\".gz\", \".Z\",\
\".z\",\".bz\", \".bz2\", \".bzip2\").")
	parser.add_option("-o","--outfile",action="store",type="string", dest="out_file",
help="The prefix of output files.")
	parser.add_option("--iformat",action="store",type="string", dest="in_format",default='fq',
help="The format of input file. Must be 'fq' or 'fa'. default=\'%default\'")
	parser.add_option("--oformat",action="store",type="string", dest="out_format",default='pdf',
help="The format of output logo file. Must be 'pdf', 'png' or 'svg'. default=\'%default\'")
	parser.add_option("-n","--nseq-limit",action="store",type="int", dest="max_seq",default=None,
help="Only process this many sequences and stop. default=%default (generate logo from\
ALL sequences).")
	parser.add_option("--font-name",action="store",type="string", dest="font_name",default='sans',
help="The font of logo characters. For a list of valid font names, run logomaker.list_font_names().\
default=\'%default\'")
	parser.add_option("--stack-order",action="store",type="string", dest="stack_order",default='big_on_top',
help="Must be 'big_on_top', 'small_on_top', or 'fixed'. 'big_on_top' : nucleotide with the highest\
frequency will be on the top; 'small_on_top' : nucleotide with the lowest frequency will be on the\
top; 'fixed' : nucleotides from top to bottom are in the same order as characters appear in the\
data frame. default=\'%default\'")
	parser.add_option("--flip-below",action="store_true", dest="flip_below",default=False,
help="If set, characters below the X-axis (which correspond to negative values in the matrix)\
will be flipped upside down. default=%default")
	parser.add_option("--shade-below",action="store",type="float", dest="shade_below", default=0.0,
help="The amount of shading to use for characters drawn below the X-axis. 0 <= shade_below <= 1.\
Larger values correspond to more shading. default=%default")
	parser.add_option("--fade-below",action="store",type="float", dest="fade_below", default=0.0,
 help="The amount of fading to use for characters drawn below the X-axis. 0 <= shade_below <= 1.\
 Larger values correspond to more fading. default=%default")
	parser.add_option("--excludeN",action="store_true",dest="exclude_N", default=False,
help="If set, exclude all DNA sequences containing \"N\".")
	parser.add_option("--highlight-start",action="store",type='int', dest="hi_start",default=None,
help="Highlight logo from this position. Must be within [0, sequence_length-1].\
default=%default (no highlight)")
	parser.add_option("--highlight-end",action="store",type='int', dest="hi_end",default=None,
help="Highlight logo to this position. Must be within [0, len(logo)-1].\
default=%default (no highlight)")
	parser.add_option("--verbose",action="store_true",dest="debug", default=False,
help="If set, print detailed information for debugging.")

	(options,args)=parser.parse_args()

	#DEGUB->INFO->WARNING->ERROR->CRITICAL
	if options.debug:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
	else:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)

	for file in ([options.in_file, options.out_file]):
		if not (file):
			parser.print_help()
			sys.exit(0)
	if options.in_format not in ['fa', 'fq']:
		logging.error("--format takes 'fq' or 'fa' as argument.")
		parser.print_help()
		sys.exit(0)
	if options.shade_below < 0 or options.shade_below > 1:
		logging.error("The shade_below value must be between 0 and 1")
		sys.exit(0)
	if options.fade_below < 0 or options.fade_below > 1:
		logging.error("The fade_below value must be between 0 and 1")
		sys.exit(0)

	if options.in_format.lower() == 'fq':
		file_iter = fastq.fastq_iter(options.in_file, mode = 'seq')
	else:
		file_iter = fastq.fasta_iter(options.in_file)
	mat = fastq.seq2countMat(file_iter, step_size=10000, exclude_N = options.exclude_N, limit = options.max_seq)
	mat.to_csv(options.out_file + '.count_matrix.csv', index=True, index_label="Index")
	fastq.make_logo(mat, outfile = options.out_file, exclude_N = options.exclude_N, font_name = options.font_name, stack_order = options.stack_order, flip_below = options.flip_below, shade_below = options.shade_below, fade_below = options.fade_below, highlight_start = options.hi_start, highlight_end = options.hi_end, oformat = options.out_format)


if __name__=='__main__':
	main()
