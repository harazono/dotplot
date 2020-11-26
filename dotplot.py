#! /usr/bin/env python3

import csv
import plotly.graph_objects as go
import plotly.express as px
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


def draw_dotplot(PAFs):
	retObj = go.Figure()
	counter = 0
	# # 63 6E FA -> rgba(99, 110, 250, 0.7)
	colorList = []
	for each_color in px.colors.qualitative.Plotly:
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
		retObj.add_trace(go.Scattergl(x = x_points, y = y_points, line=dict(width=2, color=colorList[counter]), mode='lines', name = paf[0].dsc)) #, color=one_alignment.color
		counter += 1
	retObj.update_xaxes(title_text='Query')
	retObj.update_yaxes(title_text='Reference')
	return retObj

def main():
	parser = argparse.ArgumentParser(description='Describe dot plot of alignments in PAF files. Alignments are grouped by PAF file name. This script will show dot plot on your browser and save a picture to the file which you specify.')
	parser.add_argument("Outputfilename", metavar='FileName', type=str, help='image file name.File name must be like fizz.png or fizz.svn. Please refer to https://plotly.com/python/static-image-export/')
	parser.add_argument("PAFfilename", metavar='PAF', type=str, nargs='+', help='PAF file(s)')
	args = parser.parse_args()
	paf_file_names = args.PAFfilename
	out_file_name = args.Outputfilename
	paf_instance_array = []
	for each_filename in paf_file_names:
		tmp_paf_instance = minimap2_paf_parser(each_filename, dsc = each_filename)
		paf_instance_array.append(tmp_paf_instance)
	fig = draw_dotplot(paf_instance_array)
	fig.show()
	fig.write_image(out_file_name)

if __name__ == "__main__":
	main()