#! /usr/bin/env python3

import sys
import argparse
import pprint
import dotplot
import copy
import re
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio

pp = pprint.PrettyPrinter(indent=2)

def parse_cigar(cigar):
	currentPos = 0
	

def split_paf(paf):
	altpaf = []



def main():
	parser = argparse.ArgumentParser(description='Split input paf file at large indel')
	parser.add_argument("PAFfilename", metavar='PAF', type=str, help='input PAF file name')
	parser.add_argument("-t", metavar='t', type=int, default = 10, help='Threshold of indel (10 as default). An alignmet will be split at indels greater than this threshold')
	parser.add_argument("-o", metavar='o', type=str, help='output file name (xxx_split.paf as default for xxx.paf)')

	args = parser.parse_args()
	paf_file_name = args.PAFfilename
	threshold = args.t

	raw_paf = dotplot.minimap2_paf_parser(paf_file_name)
	main_line_scatter = []
	x_points = []
	y_points = []
	for alignment in raw_paf:
		largeindel = [x for x in dotplot.parse_cigar(alignment.cigar) if x[-1] in ["D", "I"] and int(x[:-1]) > 50]
		print(largeindel) if len(largeindel) > 0 else None
		points = dotplot.cut_alignment_at_large_indel(alignment)
		for i in range(len(points)):
			x_points.append(points[i][0])
			y_points.append(points[i][1])
	#print(x_points)
	#print(y_points)
	fig = go.Figure()
	fig.add_trace(go.Scattergl(x = x_points, y = y_points, mode='lines+markers'))
	fig.update_yaxes(autorange="reversed")

	#fig.show()






if __name__ == '__main__':
	main()