#! /usr/bin/env python3

import sys
import csv

import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.io as pio
from plotly.validators.scatter.marker import SymbolValidator
import argparse
from enum import Enum
import pprint
pp = pprint.PrettyPrinter(indent=2)
import re
import json
import gzip


class paf_data():
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

def chain_parser(filename, switchflag = False):
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
				ref_start = ref_current_pos
				query_start = query_current_pos
				ref_end = ref_start + int(block[0])
				query_end = query_start + int(block[0])
				paf_alignment_data = None
				if switchflag == False:
					paf_alignment_data = paf_data(
						query_chrom = tName,
						query_start = ref_start,
						query_end = ref_end,
						ref_chrom = qName,
						ref_start = query_start,
						ref_end = query_end,
						cigar = f"{query_end - query_start}M",
						rev = qStrand,
						dsc = filename
					)
				else:
					paf_alignment_data = paf_data(
						query_chrom = qName,
						query_start = query_start,
						query_end = query_end,
						ref_chrom = tName,
						ref_start = ref_start,
						ref_end = ref_end,
						cigar = f"{query_end - query_start}M",
						rev = tStrand,
						dsc = filename
					)
				
				assert(query_end - query_start == ref_end - ref_start), "chain file reports different length between reference and query"
				ret_array.append(paf_alignment_data)
				if len(block) == 3:
					ref_current_pos = ref_end + int(block[1])
					query_current_pos = query_end + int(block[2])
	return ret_array



def minimap2_paf_parser(filename:str, dsc = None):
	paf = []
	try:
		with open(filename) as f:
			reader = csv.reader(f, delimiter = "\t")
			l = [row for row in reader]
			for each_paf in l:
				cigar = str(int(each_paf[3]) - int(each_paf[2])) + "M"
				if each_paf[-1].split(":")[0].strip() == "cg":
					cigar = each_paf[-1].split(":")[-1]

				tmp_paf = paf_data(
						ref_chrom = each_paf[5],
						ref_start = int(each_paf[7]),
						ref_end = int(each_paf[8]),
						query_chrom = each_paf[0],
						query_start = int(each_paf[2]),
						query_end = int(each_paf[3]),
						cigar = cigar,
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
				one_gtf = annotation_data(chrom, int(start), int(end), name, annotation_type, ref_or_query)
				ret_array.append(one_gtf)
	return ret_array




def gff3_attribute_split(attribute):
	ret_dict = {}
	for each_item in attribute.split(";"):
		key, val = each_item.split("=")
		ret_dict[key] = val
	return ret_dict


def gff3_parser(gff3_file_name, ref_or_query):
	return_array = []
	try:
		with gzip.open(gff3_file_name, "rt") as f:
			for each_line in f:
				if each_line.startswith("#"):
					continue
				seqid, source, _type, start, end, score, strand, phase, attribute = each_line.split("\t")
				attr = gff3_attribute_split(attribute)
				if not _type in ["biological_region", "chromosome", "supercontig", "scaffold"]:
					name = attr.get("Name")
					_id = attr["ID"].split(":")[1]
					text = f"{_id}({str(name)})"
					return_array.append(annotation_data(seqid, int(start), int(end), text, None, ref_or_query))
		return return_array

	except gzip.BadGzipFile:
		with open(gff3_file_name, "r") as f:
			for each_line in f:
				if each_line.startswith("#"):
					continue
				seqid, source, _type, start, end, score, strand, phase, attribute = each_line.split("\t")
				attr = gff3_attribute_split(attribute)
				if not _type in ["biological_region", "chromosome", "supercontig", "scaffold"]:
					name = attr.get("Name")
					_id = attr["ID"].split(":")[1]
					text = f"{_id}({str(name)})"
					return_array.append(annotation_data(seqid, int(start), int(end), text, None, ref_or_query))
		return return_array






def cigar_parser(cigar_string):
	cigar_element = re.findall(r'([0-9]+[MIDNSHPX=])', cigar_string)
	return cigar_element

def cut_alignment_at_large_indel(paf, threthold = 50):
	rev = paf.rev
	current_ref_pos = paf.ref_start
	current_query_pos = paf.query_start

	points = []
	if rev == "+":
		points.append((paf.query_start, paf.ref_start))
	else:
		current_query_pos = paf.query_end
		points.append((paf.query_end, paf.ref_start))


	for each_chain in cigar_parser(paf.cigar):
		length = int(each_chain[:-1])
		if rev == "+":
			if each_chain[-1] == "D":# and length >= threthold:
				current_ref_pos += length
			elif each_chain[-1] == "I":# and length >= threthold:
				current_query_pos += length
			else:
				current_query_pos += length
				current_ref_pos += length
		else:
			if each_chain[-1] == "D":# and length >= threthold:
				current_ref_pos += length
			elif each_chain[-1] == "I":# and length >= threthold:
				current_query_pos -= length
			else:
				current_query_pos -= length
				current_ref_pos += length

		points.append((current_query_pos, current_ref_pos))
	points.append((None, None))
	#if rev == "+":
	#	points.append((paf.query_end, paf.ref_end))
	#else:
	#	points.append((paf.query_start, paf.ref_end))
	return points




def draw_dotplot(
	PAFs, 
	chrm,
	const,
	query_centromere_breakpoint = None, 
	reference_centromere_breakpoint = None,  
	query_annotation = None, 
	ref_annotation = None, 
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

	scale_end =10
	main_line_scatter = []
	for paf in PAFs:
		x_points = []
		y_points = []
		for one_alignment in paf:#ここでcigarのパーザーを呼ぶ。
			#scale_end = max(scale_end, one_alignment.ref_start, one_alignment.ref_end, one_alignment.query_start, one_alignment.query_end)
			points = cut_alignment_at_large_indel(one_alignment)
			for i in range(len(points)):
				x_points.append(points[i][0])
				y_points.append(points[i][1])
		tmp = go.Scattergl(
			x = x_points,
			y = y_points, 
			line=dict(width=3,color=colorList[counter % 24]),
			mode='lines',
			name = paf[0].dsc.split("/")[-1]
		) #, color=one_alignment.color
		main_line_scatter.append(tmp)
		counter += 1




	query_annotation_scatter = []
	if query_annotation is not None:
		x_points = []
		y_points = []
		name_list = []
		count = 0
		for each_query_anno in query_annotation:
			if str(each_query_anno.chrom) != str(chrm):
				continue
			length = each_query_anno.end - each_query_anno.start
			x_points.append(each_query_anno.start)
			x_points.append(each_query_anno.end)
			x_points.append(None)
			y_points.append(count % 12)
			y_points.append(count % 12)
			#y_points.append(0)
			#y_points.append(0)
			y_points.append(None)
			name_list.append(each_query_anno.name)
			name_list.append(each_query_anno.name)
			name_list.append("")
			count += 1
		tmp = go.Scattergl(
				x = x_points,
				y = y_points,
				yaxis = "y2",
				line = dict(width = 2, color = "midnightblue"),
				marker = dict(color = colorList),
				marker_line_color = "midnightblue", marker_color = "midnightblue", marker_line_width = 2, marker_size = 4,
				marker_symbol = 223,
				name = "query annotation",
				text = name_list,
				opacity = 0.7,
				mode = "markers+lines",
				textposition = "bottom center",
				
				showlegend = False
				)
		query_annotation_scatter.append(tmp)
	main_line_scatter.extend(query_annotation_scatter)


	ref_annotation_scatter = []
	if ref_annotation is not None:
		x_points = []
		y_points = []
		name_list = []
		count = 0
		for each_ref_anno in ref_annotation:
			if str(each_ref_anno.chrom) != str(chrm):
				continue
			length = each_ref_anno.end - each_ref_anno.start
			y_points.append(each_ref_anno.start)
			y_points.append(each_ref_anno.end)
			y_points.append(None)
			x_points.append(count % 12)
			x_points.append(count % 12)
			#x_points.append(0)
			#x_points.append(0)
			x_points.append(None)
			name_list.append(each_ref_anno.name)
			name_list.append(each_ref_anno.name)
			name_list.append("")
			count += 1
		tmp = go.Scattergl(
				x = x_points,
				y = y_points,
				xaxis = "x2",
				line = dict(width = 2, color = "midnightblue"),
				marker = dict(color = colorList),
				marker_line_color = "midnightblue", marker_color = "midnightblue", marker_line_width = 2, marker_size = 4,
				marker_symbol = 224,
				name = each_ref_anno.name,
				text = name_list,
				opacity = 0.7,
				mode = "markers+lines",
				textposition = "bottom center",
				
				showlegend = False
				)
		ref_annotation_scatter.append(tmp)
	main_line_scatter.extend(ref_annotation_scatter)


	main_line_figure = go.Figure(data = main_line_scatter)

	if query_centromere_breakpoint:
		main_line_figure.add_vrect(x0 = query_centromere_breakpoint[0], x1 = query_centromere_breakpoint[1], fillcolor = px.colors.qualitative.Pastel[0], opacity = 0.3, layer = "above", line_width=0)
	if reference_centromere_breakpoint:
		main_line_figure.add_vrect(y0 = reference_centromere_breakpoint[0], y1 = reference_centromere_breakpoint[1], fillcolor = px.colors.qualitative.Pastel[1], opacity = 0.3, layer = "below", line_width=0)

	if chrm:
		scale_end = const["GRCh38_chromosome_length"][str(chrm)]

	main_line_figure.update_xaxes(title = {'text': "Query", "standoff": 1100}, title_font = dict(size=18), zeroline = True,  range = [0, scale_end], rangemode = "tozero", showgrid = True,  gridwidth = 1, matches = 'x', anchor = "free", position = 1)
	main_line_figure.update_yaxes(title_text = 'Reference', zeroline = True,  range = [0, scale_end], rangemode = "tozero", showgrid = True,  gridwidth = 1, scaleanchor = "x", scaleratio = 1, autorange="reversed")

	main_line_figure.update_layout(
		title = {'text': f"Chromosome {chrm}", "y": 0.95, "x": 0.5},
		legend = {"yanchor": "top", "y": 0.98, "xanchor": "right" , "x": 1.0},
		autosize = False,
		width = 1200,
		height = 1200,
		yaxis  = {"domain": [0.05, 1]},
		yaxis2 = {"domain": [0, 0.05]},
		hovermode = 'x'
		#xaxis2 = {"title": "Gene", "titlefont": {"color": "#ff7f0e"}, "tickfont": {"color": "#ff7f0e"}, "anchor": "free", "overlaying": "free", "side": "bottom", "position": 0.1}
		)



	return main_line_figure

def main():
	parser = argparse.ArgumentParser(description='Describe dot plot of alignments in PAF files. Alignments are grouped by PAF file name. This script will show dot plot on your browser and save a picture to the file which you specify.')
	parser.add_argument("PAFfilename", metavar='PAF', type=str, nargs='+', help='PAF file(s)')
	parser.add_argument("-o", metavar='FileName', type=str, help='image file name.File name must be like fizz.png or fizz.svn. Please refer to https://plotly.com/python/static-image-export/')
	parser.add_argument("-c", metavar='Chromosome', type=int, help='chromosome')
	parser.add_argument("--ref_gene", metavar='ref_gene', type=str, help='gene annotation file for reference genome')
	parser.add_argument("--query_gene", metavar='query_gene', type=str, help='gene annotation file for query genome')
	parser.add_argument("--ref_repeat", metavar='ref_repeat', type=str, help='repeat annotation file for reference genome')
	parser.add_argument("--query_repeat", metavar='query_repeat', type=str, help='repeat annotation file for query genome')
	parser.add_argument("--chain", metavar='chain', type=str, help='hg19ToHg38')
	parser.add_argument("--sf", action='store_true', help='if you want to switch ref/query in chain file, set this flag')
	args = parser.parse_args()
	paf_file_names = args.PAFfilename
	out_file_name = args.o
	chrm = args.c
	query_centromere_breakpoint = chrm
	ref_gene_file_name = args.ref_gene
	query_gene_file_name = args.query_gene
	switchflag = args.sf
	ref_gene = None
	query_gene = None

	if ref_gene_file_name is not None:
		if ".gtf" in ref_gene_file_name:
			ref_gene = gtf_parser(ref_gene_file_name, ref_or_query.R)
		elif ".gff3" in ref_gene_file_name:
			ref_gene = gff3_parser(ref_gene_file_name, ref_or_query.R)
		print(f"# of ref_anno: {len(ref_gene)}", file = sys.stderr)

	if query_gene_file_name is not None:
		if ".gtf" in query_gene_file_name:
			query_gene = gtf_parser(query_gene_file_name, ref_or_query.Q)
		elif ".gff3" in query_gene_file_name:
			query_gene = gff3_parser(query_gene_file_name, ref_or_query.Q)
		print(f"# of query_anno: {len(query_gene)}", file = sys.stderr)

	paf_instance_array = []
	for each_filename in paf_file_names:
		tmp_paf_instance = minimap2_paf_parser(each_filename, dsc = each_filename)
		paf_instance_array.append(tmp_paf_instance)

	if args.chain is not None:
		chain_paf_format_data = chain_parser(args.chain, switchflag = switchflag)
		only_designated_paf = []
		for each_chain in chain_paf_format_data:
			if each_chain.ref_chrom == f"chr{chrm}" and each_chain.query_chrom == f"chr{chrm}":
				only_designated_paf.append(each_chain)
		paf_instance_array.extend([only_designated_paf])
		#print(str(only_designated_paf[0]))
		#print(str(only_designated_paf[1]))
		#sys.exit(0)

	const = {}
	with open("const.json", "r") as f:
		const = json.load(f)

	#fig = draw_dotplot(paf_instance_array, chrm, const, query_annotation = query_gene, ref_annotation = None)
	fig = draw_dotplot(paf_instance_array, chrm, const, reference_centromere_breakpoint = None, query_annotation = query_gene, reference_annotation = ref_gene)
	#pio.kaleido.scope.default_width = 2400
	#pio.kaleido.scope.default_height = 2400
	if out_file_name is not None:
		fig.write_image(out_file_name)
	else:
		fig.show()


if __name__ == "__main__":
	main()