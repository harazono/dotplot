#! /usr/bin/env python3

import csv
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
import argparse

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

	def print_paf_info(self):
		print(f"ref_chrom: {self.ref_chrom}")
		print(f"ref_start: {self.ref_start}")
		print(f"ref_end  : {self.ref_end}")
		print(f"query_chrom: {self.query_chrom}")
		print(f"query_start: {self.query_start}")
		print(f"query_end  : {self.query_end}")
		print(f"cigar: {self.cigar}")
		print(f"rev: {self.rev}")
		print(f"dsc: {self.dsc}")


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

GRCh38_chromosome_length = {
							1: 248956422,
							2: 242193529,
							3: 198295559,
							4: 190214555,
							5: 181538259,
							6: 170805979,
							7: 159345973,
							8: 145138636,
							9: 138394717,
							10: 133797422,
							11: 135086622,
							12: 133275309,
							13: 114364328,
							14: 107043718,
							15: 101991189,
							16: 90338345,
							17: 83257441,
							18: 80373285,
							19: 58617616,
							20: 64444167,
							21: 46709983,
							22: 50818468
}

def draw_dotplot(PAFs, query_centromere_breakpoint = None, reference_centromere_breakpoint = None, chrm = None):
	retObj = go.Figure()
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
			else:
				x_points.append(one_alignment.query_end)
				x_points.append(one_alignment.query_start)
				x_points.append(None)
				y_points.append(one_alignment.ref_start)
				y_points.append(one_alignment.ref_end)
				y_points.append(None)
		retObj.add_trace(go.Scattergl(x = x_points, y = y_points, line=dict(width=3, color=colorList[counter]), mode='lines', name = paf[0].dsc)) #, color=one_alignment.color
		if query_centromere_breakpoint:
			retObj.add_vrect(x0 = query_centromere_breakpoint[0], x1 = query_centromere_breakpoint[1], fillcolor = px.colors.qualitative.Pastel[0], opacity = 0.3, layer = "below", line_width=0)
		if reference_centromere_breakpoint:
			retObj.add_vrect(y0 = reference_centromere_breakpoint[0], y1 = reference_centromere_breakpoint[1], fillcolor = px.colors.qualitative.Pastel[1], opacity = 0.3, layer = "below", line_width=0)
		counter += 1
	if chrm is not None:
		scale_end = GRCh38_chromosome_length[int(chrm)]
	else:
		scale_end = GRCh38_chromosome_length[1]
	retObj.update_xaxes(title_text='Query', range = [0, scale_end], showgrid=True, gridwidth=1)
	retObj.update_yaxes(title_text='Reference', range = [0, scale_end], showgrid=True, gridwidth=1)
	retObj.update_layout(title={'text': f"Chromosome {chrm}"})
	return retObj

def main():
	parser = argparse.ArgumentParser(description='Describe dot plot of alignments in PAF files. Alignments are grouped by PAF file name. This script will show dot plot on your browser and save a picture to the file which you specify.')
	parser.add_argument("Outputfilename", metavar='FileName', type=str, help='image file name.File name must be like fizz.png or fizz.svn. Please refer to https://plotly.com/python/static-image-export/')
	parser.add_argument("PAFfilename", metavar='PAF', type=str, nargs='+', help='PAF file(s)')
	parser.add_argument("--qc", metavar='Query_centromere', type=str, nargs=2, help='query centromere coordinates')
	parser.add_argument("--chrm", metavar='Chromosome', type=str, help='chromosome')
	args = parser.parse_args()
	paf_file_names = args.PAFfilename
	out_file_name = args.Outputfilename
	query_centromere_breakpoint = args.qc
	chrm = args.chrm
	paf_instance_array = []
	for each_filename in paf_file_names:
		tmp_paf_instance = minimap2_paf_parser(each_filename, dsc = each_filename)
		paf_instance_array.append(tmp_paf_instance)

	fig = draw_dotplot(paf_instance_array, query_centromere_breakpoint = query_centromere_breakpoint, reference_centromere_breakpoint = None, chrm = chrm)
	#fig.show()
	pio.kaleido.scope.default_width = 2400
	pio.kaleido.scope.default_height = 1200
	fig.write_image(out_file_name)

if __name__ == "__main__":
	main()