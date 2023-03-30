import argparse
import os
import logging
import pandas as pd # conda install pandas
import subprocess
import sys



def main(args):
	
	# GET BASE MODIFICATION
	mod = args.input_file.split('.')[1]


	# CREATE OUTPUT FILE
	try:
		outdir = '/'.join((args.output_path, args.output_directory))
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		prefix = '/'.join((outdir, args.savename))
		prefix = '_'.join((prefix, mod))
		outf = ".".join((prefix, 'gff'))
		f2 = open(outf, 'w')
		f2.write('##gff-version 3\n') # write header
	except:
		print('ERROR: Could not configure output GFA file. Skipping...')
		return

	# READ BED FILE
	f1 = open(args.input_file, 'r')
	count = 0 # counter for modification names
	for line in f1:
		l = line.split('\t')
		prct = float(l[10])
		# FILTER ON CONFIDENCE
		if prct >= args.threshold:
			seqname = l[0]
			source = 'Megalodon'
			feature = 'mod'
			start = l[1]
			end = str(int(l[2])-1)
			score = str(prct)
			strand = l[5]
			frame = '.'

			# CREATE ATTRIBUTE
			count = count + 1
			num = str(count).zfill(4)
			seqID = '_'.join(('mod', num))
			attributes = ('ID=', seqID, ';Name=', mod)
			attr = ''.join(attributes)

			record = (seqname, source, feature, start, end, score, strand, frame, attr)
			newRecord = '\t'.join(record)
			newRecord = ''.join((newRecord, '\n'))
			f2.write(newRecord)

	f1.close()
	f2.close()


if __name__ == "__main__":

	cwd = os.getcwd()

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='mega_bed2gff', add_help=False)
	parent_parser.add_argument('-i', '--input_file', help='Bed file (required)', required=True)
	parent_parser.add_argument('-o', '--output_directory', default='mega_bed2gff', help='Prefix of output directory', type=str)
	parent_parser.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parent_parser.add_argument('-s', '--savename', default='extract', help='Prefix for output.')
	parent_parser.add_argument('-t', '--threshold', default=80.0, help='Minimum confidence score to keep', type=float)
	subparsers = parent_parser.add_subparsers(help='sub-command help')

	args = parent_parser.parse_args()

	main(args)


