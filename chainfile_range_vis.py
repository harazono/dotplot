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
import json

pp = pprint.PrettyPrinter(indent=2)



def main():
	parser = argparse.ArgumentParser(description='convert chain file to paf format and visualize them with gene annotation')
	parser.add_argument("CHAINfilename", metavar='chain', type=str, help='Chain file name')
	parser.add_argument("grch37gene", metavar='grch37_gene', type=str, help='GRCh37 gene annotation(gff3)')
	parser.add_argument("grch38gene", metavar='grch38_gene', type=str, help='GRCh38 gene annotation(gff3)')

	args = parser.parse_args()
	chain_file_name = args.CHAINfilename
	grch37_gene_filename = args.grch37gene
	grch38_gene_filename = args.grch38gene

	chain = dotplot.chain_parser(chain_file_name)
	"""
	for each_chain in chain:
		print(each_chain)
	return None
	"""	
	grch37_gene_data = dotplot.gff3_parser(grch37_gene_filename, None)
	grch38_gene_data = dotplot.gff3_parser(grch38_gene_filename, None)

	const = {}
	with open("const.json", "r") as f:
		const = json.load(f)


	chr_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
	#chr_list = ["1"]

	current_coordinate_37 = 0
	current_coordinate_38 = 0
	start_coordinate_37 = {}
	start_coordinate_38 = {}
	length_dict_37 = const["GRCh37_chromosome_length"]
	length_dict_38 = const["GRCh38_chromosome_length"]
	chromosome_line_scatter = []
	for key in chr_list:
		cnt = 0
		if key == "X":
			cnt = 23
		elif key == "Y":
			cnt == 24
		else:
			cnt = int(key)
		length_37 = length_dict_37[str(key)]
		length_38 = length_dict_38[str(key)]
		chromosome_line_scatter.append(go.Scattergl(x=[current_coordinate_37, current_coordinate_37 + length_37], y=[-0.1, -0.1], line = dict(color = px.colors.qualitative.Pastel2[cnt % 8]), name = f"GRCh37 chr{key}"))
		chromosome_line_scatter.append(go.Scattergl(x=[current_coordinate_38, current_coordinate_38 + length_38], y=[1.1, 1.1], line = dict(color = px.colors.qualitative.Dark2[cnt % 8]), name = f"GRCh38 chr{key}"))

		start_coordinate_37[key] = current_coordinate_37
		start_coordinate_38[key] = current_coordinate_38

		current_coordinate_37 += length_37
		current_coordinate_38 += length_38


	chain_line_scatter = []

	"""
	chain_dict = {}
	for each_chain in chain:
		if each_chain.ref_chrom in chain_dict.keys():
			chain_dict[each_chain.ref_chrom].append(each_chain)
		else:
			chain_dict[each_chain.ref_chrom] = [each_chain]

	for each_ref_chromosome in chr_list:
		chain_in_this_chromosome = chain_dict[f"chr{each_ref_chromosome}"]
		for each_chain in chain_in_this_chromosome:
			if each_chain.ref_chrom != each_chain.query_chrom:
				continue
			q_key = each_chain.query_chrom.replace("chr", "")
			r_key = each_chain.ref_chrom.replace("chr", "")
			x = []
			y = []
			x.append(start_coordinate_38[r_key] + each_chain.ref_start)
			x.append(start_coordinate_37[q_key] + each_chain.query_start)
			x.append(start_coordinate_37[q_key] + each_chain.query_end)
			x.append(start_coordinate_38[r_key] + each_chain.ref_end)
			#x.append(start_coordinate_38[r_key] + each_chain.ref_start)
			x.append(None)
			y.append(1)
			y.append(0)
			y.append(0)
			y.append(1)
			#y.append(1)
			y.append(None)
			chain_line_scatter.append(go.Scattergl(x=x, y=y, mode='lines', line = dict(width=1,color="aqua"), fill = "toself"))

	"""


	
	chr_list_for_chain = [f"chr{x}" for x in chr_list]
	x = []
	y = []
	for each_chain in chain:
		if each_chain.ref_chrom not in chr_list_for_chain or each_chain.query_chrom not in chr_list_for_chain:
			#print(f"{each_chain.ref_chrom} or {each_chain.query_chrom} was ignored", file = sys.stderr)
			continue
		if each_chain.ref_chrom != each_chain.query_chrom:
			continue
		"""
		       38 start      1          4
		GRCh38 |-------------*----------*-----------------| y = 1


		GRCh37 |--------------*----------*----------------| y = 0
		       37 start       2          3
		"""
		q_key = each_chain.query_chrom.replace("chr", "")
		r_key = each_chain.ref_chrom.replace("chr", "")
		x.append(start_coordinate_38[r_key] + each_chain.ref_start)
		x.append(start_coordinate_37[q_key] + each_chain.query_start)
		x.append(start_coordinate_37[q_key] + each_chain.query_end)
		x.append(start_coordinate_38[r_key] + each_chain.ref_end)
		#x.append(start_coordinate_38[r_key] + each_chain.ref_start)
		x.append(None)
		y.append(1)
		y.append(0)
		y.append(0)
		y.append(1)
		#y.append(1)
		y.append(None)
	chain_line_scatter.append(go.Scattergl(x=x, y=y, mode='lines', line = dict(width=1,color="aqua")))
	
	chromosome_line_scatter.extend(chain_line_scatter)
	fig = go.Figure(data = chromosome_line_scatter)


	fig.update_layout(height=600, width=2400, title_text="上が38")
	fig.update_yaxes(fixedrange=True, visible=False)
	fig.update_xaxes(visible=False)


	fig.show()

if __name__ == '__main__':
	main()