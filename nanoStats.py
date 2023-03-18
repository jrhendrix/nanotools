'''
FILE:	nanoStats.py
AUTHOR:	J. R. Hendrix
URL:	https://jrhendrix.github.io/
		http://stronglab.org
DESC:	This program gets the read lenths and qualities
		using a file of fastq files.
		Reports basic statistics

''' 

#IMPORT FROM PYTHON STANDARD LIBRARY
import argparse
import os
import statistics
import sys

def main(program):
	cwd = os.getcwd()

	# TODO add exact match

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='nanoLength', add_help=False)
	parent_parser.add_argument('-i', '--input_file', help='Path to input directory')
	#parent_parser.add_argument('-o', '--output_directory', default='filtered_reads', help='Name of output directory')
	parent_parser.add_argument('-m', '--min_qscore', default=0, help='Minimum qscore to consider read', type=float)
	parent_parser.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parent_parser.add_argument('-s', '--savename', default='nano', help='Name of output file', type=str)
	subparsers = parent_parser.add_subparsers(help='sub-command help')

	args = parent_parser.parse_args()


	# CREATE TABLE OUTPUT
	try:
		outdir = args.output_path
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		prefix = '/'.join((outdir, args.savename))

		outf = "_".join((prefix, 'reads.tsv'))

		f2 = open(outf, 'w')
		header = 'Seq_name\tlength\tquality\n'
		f2.write(header)

	except:
		print('ERROR: Could not configure output file. Exit.')
		exit()

	# OPEN INPUT FILE
	try:
		f1 = open(args.input_file, 'r')
	except:
		print('ERROR: Could not open input file. Exit.')
		exit()

	lens = []
	quals = []
	count = 0
	# LOOP OVER READS IN FILE
	while True:
		# GET READ
		r_id = f1.readline().rstrip().split(' ')[0]	# line1:	ID
		if not r_id: break #EOF

		seq = f1.readline().rstrip()	# line2:	SEQUENCE
		if not seq: break  #EOF

		ignore = f1.readline().rstrip() # line3:	+
		if not ignore: break  #EOF
								

		q = f1.readline().rstrip()		# line4:	QUALITY
		if not q: break    #EOF

		length = len(seq)
		quality = sum([ord(n)-33 for n in q]) / len(q)

		if quality < args.min_qscore:
			continue
		lens.append(length)
		quals.append(quality)

		entry = (r_id, str(length), str(quality))
		record = '\t'.join(entry) + '\n'
		count = count + 1

		f2.write(record)

	f1.close()
	f2.close()

	

	# CREATE REPORT OUTPUT
	try:
		outdir = args.output_path
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		prefix = '/'.join((outdir, args.savename))

		outf = "_".join((prefix, 'stats.tsv'))

		f3 = open(outf, 'w')

	except:
		print('ERROR: Could not configure output report file. Exit.')
		exit()

	# CALCULATE READ STATS
	lAvg = statistics.mean(lens)
	lMin = min(lens)
	lMax = max(lens)
	lMed = statistics.median(lens)
	lStd = statistics.stdev(lens)
	lNum = sum(lens)

	qAvg = statistics.mean(quals)
	qMin = min(quals)
	qMax = max(quals)
	qMed = statistics.median(quals)
	qStd = statistics.stdev(quals)

	data = (f'Total reads: {str(count)}',
			f'Total bases: {str(lNum)},'
			'',
			f'length average: {str(lAvg)}',
			f'length median: {str(lMed)}',
			f'length minimum: {str(lMin)}',
			f'length maximum: {str(lMax)}',
			f'length stan dev: {str(lStd)}',
			'',
			f'quality average: {str(qAvg)}',
			f'quality median: {str(qMed)}',
			f'quality minimum: {str(qMin)}',
			f'quality maximum: {str(qMax)}',
			f'quality stan dev: {str(qStd)}')
	entry = '\n'.join(data)
	f3.write(entry)

	f3.close()

if __name__== "__main__":
	main(sys.argv[1])











