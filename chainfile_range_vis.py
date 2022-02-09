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

y_axis = {
	"chain": 1.1,
	"gene": 1.5,
	"repeat": 3,
	"chrom": 4
}

class repeat_data():
	def __init__(self, chrom, start, end, desc):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.desc = desc
	def __str__(self):
		return f"chromosome: {self.chrom}, start: {self.start}, end: {self.end}, description:{self.desc}"

def bed_parser(bed_file_name):
	retArray = []
	with open(bed_file_name, "r") as f:
		for each_line in f:
			line = each_line.strip().split("\t")
			tmpbed = repeat_data(line[0], int(line[1]), int(line[2]), line[3])
			retArray.append(tmpbed)
	return retArray



def main():
	parser = argparse.ArgumentParser(description='convert chain file to paf format and visualize them with gene annotation')
	#parser.add_argument("CHAINfilename", metavar='chain', type=str, help='Chain file name')
	#parser.add_argument("grch37gene", metavar='grch37_gene', type=str, help='GRCh37 gene annotation(gff3)')
	#parser.add_argument("grch38gene", metavar='grch38_gene', type=str, help='GRCh38 gene annotation(gff3)')
	parser.add_argument("config", metavar='config', type=str, help='config file')

	args = parser.parse_args()
	config = None
	with open(args.config, "r") as f:
		config = json.load(f)
	assert config is not None, "failed to load config.json"
	chain_from_19_file_name = config.get("hg19 to hg38")
	chain_from_38_file_name = config.get("hg38 to hg19")
	grch37_gene_filename = config.get("hg19 gene")
	grch38_gene_filename = config.get("hg38 gene")
	grch37_repeat_filename = config.get("hg19 repeat")
	grch38_repeat_filename = config.get("hg38 repeat")

	chain_from_19 = dotplot.chain_parser(chain_from_19_file_name) if chain_from_19_file_name is not None else None
	chain_from_38 = dotplot.chain_parser(chain_from_38_file_name) if chain_from_38_file_name is not None else None
	grch37_gene_data = dotplot.gff3_parser(grch37_gene_filename, None) if grch37_gene_filename is not None else None
	grch38_gene_data = dotplot.gff3_parser(grch38_gene_filename, None) if grch38_gene_filename is not None else None
	grch37_repeat_data = bed_parser(grch37_repeat_filename) if grch37_repeat_filename is not None else None
	grch38_repeat_data = bed_parser(grch38_repeat_filename) if grch38_repeat_filename is not None else None

	chr_list = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
	if config.get("chromosome list") is not None:
		chr_list = config.get("chromosome list")

	const = {}
	with open("const.json", "r") as f:
		const = json.load(f)



###Add chromosome line
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
		chromosome_line_scatter.append(go.Scattergl(x=[current_coordinate_37, current_coordinate_37 + length_37], y=[-1*y_axis["chrom"], -1*y_axis["chrom"]], line = dict(color = px.colors.qualitative.Pastel2[cnt % 8]), name = f"GRCh37 chr{key}"))
		chromosome_line_scatter.append(go.Scattergl(x=[current_coordinate_38, current_coordinate_38 + length_38], y=[y_axis["chrom"], y_axis["chrom"]], line = dict(color = px.colors.qualitative.Dark2[cnt % 8]), name = f"GRCh38 chr{key}"))

		start_coordinate_37[key] = current_coordinate_37
		start_coordinate_38[key] = current_coordinate_38

		current_coordinate_37 += length_37 + length_dict_37["1"] * 0.1
		current_coordinate_38 += length_38 + length_dict_38["1"] * 0.1


###Add chain(from 19 to 38) line
	if chain_from_19 is not None:
		chain_line_scatter = []
		chr_list_for_chain = [f"chr{x}" for x in chr_list]
		x = []
		y = []
		covered_by_chain_37_x = []
		covered_by_chain_38_x = []
		covered_by_chain_37_y = []
		covered_by_chain_38_y = []
		for each_chain in chain_from_19:
			if each_chain.ref_chrom not in chr_list_for_chain or each_chain.query_chrom not in chr_list_for_chain:
				continue
			if each_chain.ref_chrom != each_chain.query_chrom:
				continue
			"""
			       38 start      1          4
			GRCh38 |-------------*----------*-----------------| y = 1 ref
			                                                           ^
			                                                           |
			GRCh37 |--------------*----------*----------------| y = 0 query
			       37 start       2          3
			"""
			q_key = each_chain.query_chrom.replace("chr", "")
			r_key = each_chain.ref_chrom.replace("chr", "")
			x.append(start_coordinate_38[r_key] + each_chain.ref_start)
			x.append(start_coordinate_37[q_key] + each_chain.query_start)
			x.append(start_coordinate_37[q_key] + each_chain.query_end)
			x.append(start_coordinate_38[r_key] + each_chain.ref_end)
			x.append(start_coordinate_38[r_key] + each_chain.ref_start)
			x.append(None)
			y.append(1)
			y.append(-1)
			y.append(-1)
			y.append(1)
			y.append(1)
			y.append(None)
			covered_by_chain_37_x.append(start_coordinate_37[q_key] + each_chain.query_start)
			covered_by_chain_37_x.append(start_coordinate_37[q_key] + each_chain.query_end)
			covered_by_chain_37_x.append(None)
			covered_by_chain_38_x.append(start_coordinate_38[r_key] + each_chain.ref_start)
			covered_by_chain_38_x.append(start_coordinate_38[r_key] + each_chain.ref_end)
			covered_by_chain_38_x.append(None)
			covered_by_chain_37_y.append(-1*y_axis["chain"])
			covered_by_chain_37_y.append(-1*y_axis["chain"])
			covered_by_chain_37_y.append(None)
			covered_by_chain_38_y.append(y_axis["chain"])
			covered_by_chain_38_y.append(y_axis["chain"])
			covered_by_chain_38_y.append(None)

		chain_line_scatter.append(go.Scattergl(x=x, y=y, mode='lines', line = dict(width=1,color="#66FFCC"), opacity = 0.7))
		chain_line_scatter.append(go.Scattergl(x=covered_by_chain_37_x, y=covered_by_chain_37_y, line = dict(width=2,color="#66FFCC"), opacity = 0.7, name = "chain from 19 to 38", marker_line_width = 2, marker_size = 4, marker_symbol = "line-ns", mode = "markers+lines"))
		chain_line_scatter.append(go.Scattergl(x=covered_by_chain_38_x, y=covered_by_chain_38_y, line = dict(width=2,color="#66FFCC"), opacity = 0.7, name = "chain from 19 to 38", marker_line_width = 2, marker_size = 4, marker_symbol = "line-ns", mode = "markers+lines"))

	chromosome_line_scatter.extend(chain_line_scatter)


###Add chain(from 38 to 19) line
	if chain_from_38 is not None:
		chain_line_scatter = []
		chr_list_for_chain = [f"chr{x}" for x in chr_list]
		x = []
		y = []
		covered_by_chain_37_x = []
		covered_by_chain_38_x = []
		covered_by_chain_37_y = []
		covered_by_chain_38_y = []
		for each_chain in chain_from_19:
			if each_chain.ref_chrom not in chr_list_for_chain or each_chain.query_chrom not in chr_list_for_chain:
				continue
			if each_chain.ref_chrom != each_chain.query_chrom:
				continue
			"""
			       38 start      2          3
			GRCh38 |-------------*----------*-----------------| y = 1 query
			                                                          |
			                                                          v
			GRCh37 |--------------*----------*----------------| y = 0 ref
			       37 start       1          4
			"""
			q_key = each_chain.query_chrom.replace("chr", "")
			r_key = each_chain.ref_chrom.replace("chr", "")
			x.append(start_coordinate_37[r_key] + each_chain.query_start)
			x.append(start_coordinate_38[q_key] + each_chain.ref_start)
			x.append(start_coordinate_38[q_key] + each_chain.ref_end)
			x.append(start_coordinate_37[r_key] + each_chain.query_end)
			x.append(start_coordinate_37[r_key] + each_chain.query_start)
			x.append(None)
			y.append(-1)
			y.append(1)
			y.append(1)
			y.append(-1)
			y.append(-1)
			y.append(None)
			covered_by_chain_37_x.append(start_coordinate_37[r_key] + each_chain.query_start)
			covered_by_chain_37_x.append(start_coordinate_37[r_key] + each_chain.query_end)
			covered_by_chain_37_x.append(None)
			covered_by_chain_38_x.append(start_coordinate_38[q_key] + each_chain.ref_start)
			covered_by_chain_38_x.append(start_coordinate_38[q_key] + each_chain.ref_end)
			covered_by_chain_38_x.append(None)
			covered_by_chain_37_y.append(-1*y_axis["chain"])
			covered_by_chain_37_y.append(-1*y_axis["chain"])
			covered_by_chain_37_y.append(None)
			covered_by_chain_38_y.append(y_axis["chain"])
			covered_by_chain_38_y.append(y_axis["chain"])
			covered_by_chain_38_y.append(None)

		chain_line_scatter.append(go.Scattergl(x=x, y=y, mode='lines', line = dict(width=1,color="#FF6699"), opacity = 0.7))
		chain_line_scatter.append(go.Scattergl(x=covered_by_chain_37_x, y=covered_by_chain_37_y, line = dict(width=2,color="#FF6699"), opacity = 0.7, name = "chain from 38 to 19", marker_line_width = 2, marker_size = 4, marker_symbol = "line-ns", mode = "markers+lines"))
		chain_line_scatter.append(go.Scattergl(x=covered_by_chain_38_x, y=covered_by_chain_38_y, line = dict(width=2,color="#FF6699"), opacity = 0.7, name = "chain from 38 to 19", marker_line_width = 2, marker_size = 4, marker_symbol = "line-ns", mode = "markers+lines"))

	chromosome_line_scatter.extend(chain_line_scatter)


### Add gene(19)
	if grch37_gene_data is not None:
		grch37_annotation_scatter = []
		grch37_x = []
		grch37_y = []
		name_list = []
		count = 0
		for each_anno in grch37_gene_data:
			chrm = each_anno.chrom
			if not chrm in chr_list:
				continue
			offset = start_coordinate_37[chrm]
			length = each_anno.end - each_anno.start
			grch37_x.append(offset + each_anno.start)
			grch37_x.append(offset + each_anno.end)
			grch37_x.append(None)
			grch37_y.append(-1*y_axis["gene"] - (count % 10) / 10.0)
			grch37_y.append(-1*y_axis["gene"] - (count % 10) / 10.0)
			#grch37_y.append(0)
			#grch37_y.append(0)
			grch37_y.append(None)
			name_list.append(each_anno.name)
			name_list.append(each_anno.name)
			name_list.append("")
			count += 1
		tmp = go.Scattergl(
				x = grch37_x,
				y = grch37_y,
				line = dict(width = 2, color = "midnightblue"),
				marker_line_color = "midnightblue", marker_color = "midnightblue", marker_line_width = 2, marker_size = 4,
				marker_symbol = "line-ns",
				name = "grch37 annotation",
				text = name_list,
				opacity = 0.7,
				mode = "markers+lines",
				textposition = "bottom center",
				showlegend = False
				)
		chromosome_line_scatter.extend([tmp])


### Add gene(38)
	if grch38_gene_data is not None:
		grch38_annotation_scatter = []
		grch38_x = []
		grch38_y = []
		name_list = []
		count = 0
		for each_anno in grch38_gene_data:
			chrm = each_anno.chrom
			if not chrm in chr_list:
				continue
			offset = start_coordinate_38[chrm]
			length = each_anno.end - each_anno.start
			grch38_x.append(offset + each_anno.start)
			grch38_x.append(offset + each_anno.end)
			grch38_x.append(None)
			grch38_y.append(y_axis["gene"] + (count % 10) / 10.0)
			grch38_y.append(y_axis["gene"] + (count % 10) / 10.0)
			grch38_y.append(None)
			name_list.append(each_anno.name)
			name_list.append(each_anno.name)
			name_list.append("")
			count += 1
		tmp = go.Scattergl(
				x = grch38_x,
				y = grch38_y,
				line = dict(width = 2, color = "midnightblue"),
				marker_line_color = "midnightblue", marker_color = "midnightblue", marker_line_width = 2, marker_size = 4,
				marker_symbol = "line-ns",
				name = "grch38 annotation",
				text = name_list,
				opacity = 0.7,
				mode = "markers+lines",
				textposition = "bottom center",
				showlegend = False
				)
		chromosome_line_scatter.extend([tmp])


### Add repeat(19)
	if grch37_repeat_data is not None:
		repeat_x = []
		repeat_y = []
		name_list = []
		for each_bed in grch37_repeat_data:
			if each_bed.chrom not in chr_list_for_chain:
				continue
			chrm = each_bed.chrom.replace("chr", "")
			if chrm == "X":
				chrom = 23
			if chrm == "Y":
				chrom = 24
			offset = start_coordinate_37[chrm]
			repeat_x.append(offset + each_bed.start)
			repeat_x.append(offset + each_bed.end)
			repeat_x.append(None)
			repeat_y.append(-1 * y_axis["repeat"])
			repeat_y.append(-1 * y_axis["repeat"])
			repeat_y.append(None)
			name_list.append(each_bed.desc)
			name_list.append(each_bed.desc)
			name_list.append(None)
		tmp = go.Scattergl(
				x = repeat_x,
				y = repeat_y,
				line = dict(width = 5, color = "midnightblue"),
				#marker_line_color = "midnightblue", marker_color = "midnightblue", marker_line_width = 2, marker_size = 4,
				#marker_symbol = "line-ns",
				name = "grch37 repeat",
				text = name_list,
				opacity = 0.7,
				mode = "lines",
				textposition = "bottom center",
				showlegend = False
				)
		chromosome_line_scatter.extend([tmp])



### Add repeat(38)
	if grch38_repeat_data is not None:
		repeat_x = []
		repeat_y = []
		name_list = []
		for each_bed in grch38_repeat_data:
			if each_bed.chrom not in chr_list_for_chain:
				continue
			chrm = each_bed.chrom.replace("chr", "")
			if chrm == "X":
				chrom = 23
			if chrm == "Y":
				chrom = 24
			offset = start_coordinate_38[chrm]
			repeat_x.append(offset + each_bed.start)
			repeat_x.append(offset + each_bed.end)
			repeat_x.append(None)
			repeat_y.append(y_axis["repeat"])
			repeat_y.append(y_axis["repeat"])
			repeat_y.append(None)
			name_list.append(each_bed.desc)
			name_list.append(each_bed.desc)
			name_list.append(None)
		tmp = go.Scattergl(
				x = repeat_x,
				y = repeat_y,
				line = dict(width = 5, color = "midnightblue"),
				#marker_line_color = "midnightblue", marker_color = "midnightblue", marker_line_width = 2, marker_size = 4,
				#marker_symbol = "line-ns",
				name = "grch38 repeat",
				text = name_list,
				opacity = 0.7,
				mode = "lines",
				textposition = "bottom center",
				showlegend = False
				)
		chromosome_line_scatter.extend([tmp])




	fig = go.Figure(data = chromosome_line_scatter)

	fig.update_layout(height=600, width=2000, title_text="上が38")
	fig.update_yaxes(fixedrange=True, visible=False)

	#fig.update_xaxes(visible=False)


	fig.show()

if __name__ == '__main__':
	main()