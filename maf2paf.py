#! /usr/bin/env python3
import dotplot
import sys

def last_maf_parser(filename:str, dsc = None, number = 0):
	tmp_array = []
	with open(filename) as f:
		for line in f:
			if line.startswith("#") or line == "\n":
				pass
				#print(line, end="")
			else:
				if line.startswith("s"):
					s = line.split()[1:-1]
					tmp_array.append(s)
	for i in range(0, len(tmp_array), 2):
		ref_idx = i
		qry_idx = i + 1
		ref = tmp_array[ref_idx]
		qry = tmp_array[qry_idx]
		ref_chrm, ref_start_pos, ref_aligned_bases, ref_strand, ref_seqlen = ref
		qry_chrm, qry_start_pos, qry_aligned_bases, qry_strand, qry_seqlen = qry
		ref_end_pos = int(ref_start_pos) + int(ref_aligned_bases)
		if qry_strand == '+':
			qry_end_pos = int(qry_start_pos) + int(qry_aligned_bases)
		else:
			qry_end_pos = int(qry_start_pos) - int(qry_aligned_bases)
		ref_end_pos = str(ref_end_pos)
		qry_end_pos = str(qry_end_pos)
		print_array = [qry_chrm, qry_seqlen, qry_start_pos, qry_end_pos, qry_strand, ref_chrm, ref_seqlen, ref_start_pos, ref_end_pos]
		print_str = "\t".join(print_array)
		print(print_str)


def main():
	if len(sys.argv) < 2:
		print(f"Usage: maf2paf.py yourAlignment.maf > yourAlignment.paf\nThis script will output result to stdout.")
	maf = last_maf_parser(sys.argv[1])

if __name__ == "__main__":
	main()