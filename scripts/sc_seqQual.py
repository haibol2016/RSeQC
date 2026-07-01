#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program generates heatmap from a FASTQ file to visualize the sequencing quality.
"""


import sys
import logging
from qcmodule import fastq
from qcmodule import heatmap
from optparse import OptionParser

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
	parser.add_option("-i","--infile", action="store", type="string", dest="in_file",
		help="Input file in FASTQ (https://en.wikipedia.org/wiki/FASTQ_format#) format.")
	parser.add_option("-o","--outfile",action="store",type="string", dest="out_file",
		help="The prefix of output files.")
	parser.add_option("-n","--nseq-limit",action="store",type="int", dest="max_seq",default=None,
		help="Only process this many sequences and stop. default=%default (generate logo from ALL sequences).")
	parser.add_option("--cell-width",action="store",type="int", dest="cell_width",default = 12,
		help="Cell width (in points) of the heatmap. default=%default")
	parser.add_option("--cell-height",action="store",type="int", dest="cell_height",default = 10,
		help="Cell height (in points) of the heatmap. default=%default")
	parser.add_option("--font-size",action="store",type="int", dest="font_size", default = 6,
		help="Font size in points. If --display-num was set, fontsize_number = 0.8 * font_size. default=%default")
	parser.add_option("--angle",action="store",type="int", dest="col_angle",default = 45,
		help="The angle (must be 0, 45, 90, 270, 315) of column text lables under the heatmap. default=%default")
	parser.add_option("--text-color",action="store",type="string", dest="text_color", default = 'black',
		help="The color of numbers in each cell. default=%default")
	parser.add_option("--file-type",action="store",type="string", dest="file_type", default = 'pdf',
		help="The file type of heatmap. Choose one of 'pdf', 'png', 'tiff', 'bmp', 'jpeg'. default=%default")
	parser.add_option("--no-num",action="store_true", dest="no_num", default=False,
		help="if set, will not print numerical values to cells. default=%default")
	parser.add_option("--verbose",action="store_true",dest="debug",default=False,
		help="If set, will produce detailed information for debugging.")

	(options,args)=parser.parse_args()
	if options.debug:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
	else:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)

	for file in ([options.in_file, options.out_file]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	file_iter = fastq.fastq_iter(options.in_file, mode = 'qual')
	qual_mat = fastq.qual2countMat(file_iter, limit=options.max_seq)
	qual_mat = qual_mat.T
	qual_mat.sort_index(inplace=True, ascending=False)

	logging.debug("Sequence quality score matrix (raw reads count)")
	if options.debug:
		print (qual_mat)
	qual_mat_per = (qual_mat/qual_mat.sum())

	#for c in qual_mat_per.columns:
	#	qual_mat_per[c] = qual_mat_per[c].apply(lambda x: ("%.2f" % x).lstrip('0'))

	logging.debug("Sequence quality score matrix (percent of reads)")
	if options.debug:
		print (qual_mat_per)

	qual_mat.to_csv(options.out_file + '.qual_count.csv', index=True, index_label="Index")
	qual_mat_per.to_csv(options.out_file + '.qual_percent.csv', index=True, index_label="Index")


	heatmap.make_heatmap(infile = options.out_file + '.qual_percent.csv', outfile = options.out_file + '.qual_heatmap',filetype= options.file_type, cell_width = options.cell_width, cell_height = options.cell_height, col_angle = options.col_angle, font_size = options.font_size, text_color = options.text_color, no_numbers = options.no_num)


if __name__=='__main__':
	main()
