#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program generates heatmaps to visualize the **positions** (X-axis),
type of edits (Y-axis, such as "C" to "T") and  **frequencies** (color)
of error-corrected nucleotides in cell barcodes and UMIs.

"""

import sys
import logging
from qcmodule import scbam
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
__status__ = "Production"


def main():
	usage = __doc__
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--infile",action="store",type="string", dest="in_file",help="Input file in BAM foramt.")
	parser.add_option("-o","--outfile",action="store",type="string", dest="out_file",help="The prefix of output files.")
	parser.add_option("--limit",action="store",type="int", dest="reads_num",default = None, help="Number of alignments to process. default=%default")
	parser.add_option("--cr-tag",action="store",type="string", dest="CR_tag", default='CR', help="Tag of cellular barcode reported by the sequencer in BAM file. default=\'%default\'")
	parser.add_option("--cb-tag",action="store",type="string", dest="CB_tag", default='CB', help="Tag of error-corrected cellular barcode in BAM file. default=\'%default\'")
	parser.add_option("--ur-tag",action="store",type="string", dest="UR_tag", default='UR', help="Tag of UMI reported by the sequencer in BAM file. default=\'%default\'")
	parser.add_option("--ub-tag",action="store",type="string", dest="UB_tag", default='UB', help="Tag of error-corrected UMI in BAM file. default=\'%default\'")
	parser.add_option("--cell-width",action="store",type="int", dest="cell_width",default = 15, help="Points of cell width in the heatmap. default=%default")
	parser.add_option("--cell-height",action="store",type="int", dest="cell_height",default = 10, help="Points of cell height in the heatmap. default=%default")
	parser.add_option("--font-size",action="store",type="int", dest="font_size", default = 8, help="Font size. If --display-num was set, fontsize_number = 0.8 * font_size. default=%default")
	parser.add_option("--angle",action="store",type="int", dest="col_angle",default = 45, help="The angle (must be 0, 45, 90, 270, 315) of column text lables under the heatmap. default=%default")
	parser.add_option("--text-color",action="store",type="string", dest="text_color", default = 'black', help="The color of numbers in each cell. default=%default")
	parser.add_option("--file-type",action="store",type="string", dest="file_type", default = 'pdf', help="The file type of heatmap. Choose one of 'pdf', 'png', 'tiff', 'bmp', 'jpeg'. default=%default")
	parser.add_option("--verbose",action="store_true",dest="debug",default=False,help="If set, detailed running information is printed to screen.")
	parser.add_option("--no-num",action="store_true", dest="no_num", default=False, help="If set, will not print numerical values to cells. default=%default")

	(options,args)=parser.parse_args()
	if options.debug:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
	else:
		logging.basicConfig(format = "%(asctime)s [%(levelname)s]  %(message)s",datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)


	for file in ([options.in_file, options.out_file]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	scbam.barcode_edits(infile = options.in_file, outfile = options.out_file, limit = options.reads_num, CR_tag = options.CR_tag, CB_tag = options.CB_tag, UR_tag = options.UR_tag, UB_tag = options.UB_tag)

	CB_mat_file = options.out_file + '.CB_edits_count.csv'
	CB_heatmap_file = options.out_file + '.CB_edits_heatmap'
	heatmap.make_heatmap(infile = CB_mat_file, outfile = CB_heatmap_file, filetype=options.file_type, cell_width=options.cell_width, cell_height=options.cell_height, col_angle=options.col_angle, font_size=options.font_size, text_color=options.text_color, no_numbers = options.no_num, log2_scale=True)

	UMI_mat_file = options.out_file + '.UMI_edits_count.csv'
	UMI_heatmap_file = options.out_file + '.UMI_edits_heatmap'
	heatmap.make_heatmap(infile = UMI_mat_file, outfile = UMI_heatmap_file, filetype=options.file_type, cell_width=options.cell_width, cell_height=options.cell_height, col_angle=options.col_angle, font_size=options.font_size, text_color=options.text_color, no_numbers = options.no_num, log2_scale=True)


if __name__=='__main__':
	main()
