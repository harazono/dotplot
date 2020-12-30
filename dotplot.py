#! /usr/bin/env python3

import sys
import csv
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
import argparse
import const
from enum import Enum
import pprint
pp = pprint.PrettyPrinter(indent=2)

class pafData():
	def __init__(self, ref_chrom, ref_start, ref_end, query_chrom, query_start, query_end, cigar, rev, dsc = None, color = None):
		self.ref_chrom = ref_chrom
		self.ref_start = ref_start
		self.ref_end = ref_end
		self.query_chrom = query_chrom
		self.query_start = query_start
		self.query_end = query_end
		self.cigar = cigar
		self.rev = rev
		self.dsc = dsc

	def __str__(self):
		return f"ref_chrom: {self.ref_chrom}, ref_start: {self.ref_start}, ref_end: {self.ref_end}, query_chrom: {self.query_chrom}, query_start: {self.query_start}, query_end: {self.query_end}, cigar: {self.cigar}, rev: {self.rev}, dsc: {self.dsc}"

class ref_or_query(Enum):
	R = 0
	Q = 1

class annotation_data():
	def __init__(self, chrom, start, end, name, annotation_type, ref_or_query):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.name = name
		self.annotation_type = annotation_type
		self.ref_or_query = ref_or_query
	def __str__(self):
		return f"{self.ref_or_query}\t{self.chrom}:{self.start}-{self.end}\t{self.annotation_type}\t{self.name}"

def chain_parser(filename, chrom):
	print(f"file {filename}, chrom {chrom}")
	ret_array = []
	score = 0
	tName = ""
	tSize = 0
	tStart = 0
	tEnd = 0
	qName = ""
	qSize = 0
	qStrand = ""
	qStart = 0
	qEnd = 0
	id_ = 0
	ref_current_pos = 0
	query_current_pos = 0
	with open(filename) as f:
		for line in f:
			if line.startswith("chain"):
				chain, score, tName, tSize, tStrand, tStart, tEnd, qName, qSize, qStrand, qStart, qEnd, id_ = line.strip().split()
				ref_current_pos = int(tStart)
				query_current_pos = int(qStart)
			elif line != "\n":
				block = line.strip().split()
				ref_start = int(ref_current_pos)
				query_start = int(query_current_pos)
				ref_end = ref_start + int(block[0])
				query_end = query_start + int(block[0])
				paf_alignment_data = pafData(
					ref_chrom = tName,
					ref_start = ref_start,
					ref_end = ref_end,
					query_chrom = qName,
					query_start = query_start,
					query_end = query_end,
					cigar = None,
					rev = qStrand,
					)
				ret_array.append(paf_alignment_data)
				if len(block) == 3:
					ref_current_pos += int(block[0]) + int(block[1])
					query_current_pos += int(block[0]) + int(block[2])
	return ret_array



def minimap2_paf_parser(filename:str, dsc = None):
	paf = []
	try:
		with open(filename) as f:
			reader = csv.reader(f, delimiter = "\t")
			l = [row for row in reader]
			for each_paf in l:
				tmp_paf = pafData(
						ref_chrom = each_paf[5],
						ref_start = int(each_paf[7]),
						ref_end = int(each_paf[8]),
						query_chrom = each_paf[0],
						query_start = int(each_paf[2]),
						query_end = int(each_paf[3]),
						cigar = None,
						rev = each_paf[4],
						dsc = dsc)
				paf.append(tmp_paf)
	except:
		pass
	return paf

def gtf_parser(gtf_file_name, ref_or_query):
	ret_array = []
	with open(gtf_file_name, "r") as f:
		for cols in csv.reader(f, delimiter='\t'):
			if cols[0].startswith("#"):
				continue
			else:
				chrom = cols[0]
				start = int(cols[3])
				end = int(cols[4])
				annotation_type = cols[2]
				name = cols[8].split(";")[0].split()[1].replace("\"", "")
				one_gtf = annotation_data(chrom, start, end, name, annotation_type, ref_or_query)
				ret_array.append(one_gtf)
	return ret_array



def draw_dotplot(
	PAFs, 
	query_centromere_breakpoint = None, 
	reference_centromere_breakpoint = None, 
	chrm = None, 
	query_annotation = None, 
	reference_annotation = None
	):
	counter = 0
	# # 63 6E FA -> rgba(99, 110, 250, 0.7)
	colorList = []
	for each_color in px.colors.qualitative.Dark24:
		r_str = each_color[1:3]
		g_str = each_color[3:5]
		b_str = each_color[5:7]
		r_dec = int(r_str, 16)
		g_dec = int(g_str, 16)
		b_dec = int(b_str, 16)
		colorList.append(f"rgba({r_dec}, {g_dec}, {b_dec}, 0.5)")

	main_line_scatter = []
	for paf in PAFs:
		x_points = []
		y_points = []
		for one_alignment in paf:
			if one_alignment.rev == "+":
				x_points.append(one_alignment.query_start)
				x_points.append(one_alignment.query_end)
				x_points.append(None)
				y_points.append(one_alignment.ref_start)
				y_points.append(one_alignment.ref_end)
				y_points.append(None)
			elif one_alignment.rev == "-":
				x_points.append(one_alignment.query_end)
				x_points.append(one_alignment.query_start)
				x_points.append(None)
				y_points.append(one_alignment.ref_start)
				y_points.append(one_alignment.ref_end)
				y_points.append(None)
		tmp = go.Scattergl(x = x_points, y = y_points, line=dict(width=3, color=colorList[counter]), mode='lines', name = paf[0].dsc) #, color=one_alignment.color
		main_line_scatter.append(tmp)
		counter += 1

	if chrm is not None:
		scale_end = const.GRCh38_chromosome_length[int(chrm)]
	else:
		scale_end = const.GRCh38_chromosome_length[1]

	query_annotation_scatter = []
	if query_annotation:
		x_points = []
		y_points = []
		name_list = []
		for each_query_anno in query_annotation:
			x_points = []
			y_points = []
			x_points.append(each_query_anno.start)
			x_points.append(each_query_anno.end)
			#x_points.append(None)
			y_points.append(0)
			y_points.append(0)
			#y_points.append(None)
			#name_list.append(each_query_anno.name)
			tmp = go.Scattergl(
				x = x_points,
				y = y_points,
				yaxis = "y2",
				line = dict(width = 40, color = "black"),
				name = each_query_anno.name,
				opacity = 0.25,
				showlegend = False
				#fill='toself'
				#fillcolor='rgba(0,100,80,0.2)',
				#line_color='rgba(255,255,255,0)'
				)
			query_annotation_scatter.append(tmp)
		#tmp = go.Scattergl(x = x_points, y = y_points, yaxis = "y2", line=dict(width = 30, color = "black"), name = "Gene", opacity = 1, showlegend = False) #, color=one_alignment.color
		#query_annotation_scatter.append(tmp)

	main_line_scatter.extend(query_annotation_scatter)
	main_line_figure = go.Figure(data = main_line_scatter)

	if query_centromere_breakpoint:
		main_line_figure.add_vrect(x0 = query_centromere_breakpoint[0], x1 = query_centromere_breakpoint[1], fillcolor = px.colors.qualitative.Pastel[0], opacity = 0.3, layer = "above", line_width=0)
	if reference_centromere_breakpoint:
		main_line_figure.add_vrect(y0 = reference_centromere_breakpoint[0], y1 = reference_centromere_breakpoint[1], fillcolor = px.colors.qualitative.Pastel[1], opacity = 0.3, layer = "below", line_width=0)


	main_line_figure.update_xaxes(title = {'text': "Query", "standoff": 1000},     zeroline = True,  range = [0, scale_end], rangemode = "tozero", showgrid = True,  gridwidth = 1, matches = 'x')
	main_line_figure.update_yaxes(title_text = 'Reference', zeroline = True,  range = [0, scale_end], rangemode = "tozero", showgrid = True,  gridwidth = 1, scaleanchor = "x", scaleratio = 1)
	#retObj.update_xaxes(title_text = 'Gene',      zeroline = True,  range = [0, scale_end], rangemode = "tozero", showgrid = True,  gridwidth = 1, matches = 'x')
	#retObj.update_yaxes(title_text = ''         , zeroline = False, range = [-0.5, 0.5],    rangemode = "tozero", showgrid = False, gridwidth = 0, matches = 'y2')

	main_line_figure.update_layout(
		title = {'text': f"Chromosome {chrm}", "y": 0.95, "x": 0.5},
		legend = {"yanchor": "top", "y": 0.99, "xanchor": "left" , "x": 0.01},
		autosize = False,
		width = 1200,
		height = 1200,
		yaxis  = {"domain": [0.05, 1]},
		yaxis2 = {"domain": [0, 0.05]}
		#xaxis2 = {"title": "Gene", "titlefont": {"color": "#ff7f0e"}, "tickfont": {"color": "#ff7f0e"}, "anchor": "free", "overlaying": "free", "side": "bottom", "position": 0.1}
		)



	return main_line_figure

def main():
	parser = argparse.ArgumentParser(description='Describe dot plot of alignments in PAF files. Alignments are grouped by PAF file name. This script will show dot plot on your browser and save a picture to the file which you specify.')
	parser.add_argument("Outputfilename", metavar='FileName', type=str, help='image file name.File name must be like fizz.png or fizz.svn. Please refer to https://plotly.com/python/static-image-export/')
	parser.add_argument("PAFfilename", metavar='PAF', type=str, nargs='+', help='PAF file(s)')
	#parser.add_argument("--qc", metavar='Query_centromere', type=str, nargs=2, help='query centromere coordinates')
	parser.add_argument("--chrm", metavar='Chromosome', type=int, help='chromosome')
	parser.add_argument("--ref_gene", metavar='ref_gene', type=str, help='gene annotation file for reference genome')
	parser.add_argument("--query_gene", metavar='query_gene', type=str, help='gene annotation file for query genome')
	parser.add_argument("--ref_repeat", metavar='ref_repeat', type=str, help='repeat annotation file for reference genome')
	parser.add_argument("--query_repeat", metavar='query_repeat', type=str, help='repeat annotation file for query genome')
	parser.add_argument("--chain", metavar='chain', type=str, help='hg19ToHg38')
	args = parser.parse_args()
	paf_file_names = args.PAFfilename
	out_file_name = args.Outputfilename
	chrm = args.chrm
	query_centromere_breakpoint = chrm
	ref_gene_file_name = args.ref_gene
	query_gene_file_name = args.query_gene

	ref_gene = None
	query_gene = None
	if ref_gene_file_name is not None:
		ref_gene = gtf_parser(ref_gene_file_name, ref_or_query.R)
		print(f"# of ref_anno: {len(ref_gene)}", file = sys.stderr)
	if query_gene_file_name is not None:
		query_gene = gtf_parser(query_gene_file_name, ref_or_query.Q)
		print(f"# of query_anno: {len(query_gene)}", file = sys.stderr)

	paf_instance_array = []
	for each_filename in paf_file_names:
		tmp_paf_instance = minimap2_paf_parser(each_filename, dsc = each_filename)
		paf_instance_array.append(tmp_paf_instance)

	chain_paf_format_data = chain_parser(args.chain, chrm)
	only_designated_paf = []
	for each_chain in chain_paf_format_data:
		if each_chain.ref_chrom == f"chr{chrm}" and each_chain.query_chrom == f"chr{chrm}":
			only_designated_paf.append(each_chain)
	print(len(only_designated_paf))
	fig = draw_dotplot([only_designated_paf])
	#fig = draw_dotplot(paf_instance_array, query_centromere_breakpoint = const.GRCh37_centromere_coordinates[chrm], reference_centromere_breakpoint = None, chrm = chrm, query_annotation = query_gene, reference_annotation = ref_gene)
	fig.show()
	pio.kaleido.scope.default_width = 2400
	pio.kaleido.scope.default_height = 2400
	#fig.write_image(out_file_name)

if __name__ == "__main__":
	main()