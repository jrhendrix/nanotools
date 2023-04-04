'''
FILE:	runFaidx.py
AUTHOR:	J.R. Hendrix
URL: 	http://stronglab.org
DESC:	This script that wraps samtools faidx
		using more flexible input formats

'''

# IMPORT FROM PYTHON STANDARD LIBRARY
import argparse
import glob
import logging
import os
import subprocess
import sys


class Dir:
	""" Base class for system directories """

	def __init__(self, path):
		self._path = None
		self.path = path

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isdir(value):
			self._path = value
		else:
			raise NotADirectoryError(value)

	@property
	def dirname(self):
		return self.path.strip("/").split("/")[-1]

	@property
	def children(self):
		children = [Dir(os.path.join(self.path, subdir)) 
			for subdir in os.listdir(self.path) 
			if os.path.isdir(os.path.join(self.path, subdir))]
		if len(children) > 0:
			return children
		else:
			return None

	@property
	def files(self):
		files = [File(os.path.join(self.path, file))
			for file in os.listdir(self.path)
			if os.path.isfile(os.path.join(self.path, file))]
		if len(files) > 0:
			return files
		else:
			return None

	def join(self, *args):
		return os.path.join(self.path, *args)

	def make_subdir(self, *args):
		""" Makes recursive subdirectories from 'os.path.join' like arguments """
		subdir = self.join(*args)
		return self.make(subdir)

	@classmethod
	def make(cls, path):
		try:
			os.makedirs(path)
			return cls(path)
		except FileExistsError:
			return cls(path)

	def __repr__(self):
		return self.path
	


class File:
	""" Base class for all file-types """

	def __init__(self, path, file_type=None):
		self._path = None
		self.path = path
		self.file_type = file_type

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isfile(value):
			self._path = value
		else:
			raise FileNotFoundError(value)

	@property
	def dir(self):
		return Dir(os.path.dirname(self.path))

	@property
	def filename(self):
		return os.path.basename(self.path)

	@property
	def file_prefix(self):
		return self.filename.split(".")[0]

	@property
	def extension(self):
		return self.filename.split(".")[-1]
	


def get_keys_from_gene_name(args, ifile):
	''' 
	Ientifies keys from gene name and Prokka output (.tsv)
	'''

	keys = []

	gene = args.gene

	# LOOP THROUGH GENE NAMES
	command = ['grep', gene, ifile.path]
	process = subprocess.Popen(command, shell=False, stdout=subprocess.PIPE)
	result = process.stdout.readlines()
	process.stdout.close()
	process.wait()

	# if no result, continue
	if len(result) < 1:
		print()
		return keys

	# extract sequence ID
	for line in result:
		line = line.decode("utf-8")
		line = line.split('\t')

		if line[1] != 'CDS':
			continue
		gName = line[3]
		if args.match == 'exact':
			if gName != gene:
				continue
		elif args.match == 'gene':
			if gName.split('_')[0] != gene:
				continue

		keys.append(line[0])

	return keys


def get_keys_from_table(args, ifile, suffix):

	f = open(ifile.path, 'r')
	next(f) # skip header line
	for line in f:
		keys = []
		l = line.strip().split('\t')
		fname = l[0]
		for i in range(1, len(l)):
			keys.append(l[i])

		extract_seqs(args, keys, suffix, fname)
	return

def get_keys_from_roary(args, ifile, suffix):
	f = open(ifile.path, 'r')

	next(f) # Skip header line
	for line in f:
		l = line.split(',')
		fname = l[0].strip('"')
		print('\n')
		try:
			keyList = l[14:]
		except:
			continue

		print(keyList)
		# REMOVE EMPTY ELEMENTS
		keys = []
		for e in keyList:
			e = e.strip('"')
			if ' ' in e or ':' in e:
				continue
			if e:
				keys.append(e)

		# EXTRACT SEQUENCES
		extract_seqs(args, keys, suffix, fname)

	f.close()
	return


def get_keys_for_pangenome(args, ifile):

	f = open(ifile.path, 'r')
	accKeys = []
	coreKeys = []
	#next(f) # Skip header line
	# GET NUMBER OF SAMPLES
	header = f.readline()
	h = header.split(',')[14:]
	header.pop(-1)
	numSamps = len(header)

	for line in f:
		try:
			l = line.split(',')[14:]
			l.pop(-1)
		except:
			continue

		# REMOVE EMPTY ELEMENTS
		myIDs = []
		n = 0
		for e in l:
			e = e.strip('"')
			if ' ' in e or ':' in e:
				continue
			if e:
				myIDs.append(e)
				n = n + 1

		if len(myIDs) == numSamps:
			coreKeys.append(myIDs[0])
		elif len(myIDs) < numSamps:
			accKeys.append(myIDs[0])
		else:
			print('ERROR: Something went wrong.')
		
	f.close()

	return accKeys, coreKeys

def get_keys_from_file(args, ifile):
	''' Extract keys from a file. Return list of keys '''

	# try to open file
	f = open(ifile.path, 'r')
	keys = []
	for line in f:
		#if line.startswith('>'):	# identify header
		e = str(line.split()[0])	# split the header and get unique ID
		keys.append(e)

	print(keys)
	f.close()
	return keys


def extract_seqs(args, keys, suffix, tag=None):
	''' Accepts a list of keys and uses samtools faidx to extract fasta elements matching key 
	Writes hit(s) to file '''

	# use input file for output
		# preserves origin of data
		# preserves data type in extension

	TOP_DIR = Dir(args.output_path)
	OUTDIR = TOP_DIR.make_subdir(args.output_directory)

	
	if tag is not None:
		outname = '_'.join((args.savename, tag))
	else:
		outname = args.savename
	outname = '.'.join((outname, suffix))
	outfile = '/'.join((OUTDIR.path, outname))

	try:
		f = open(outfile, 'w')

	except IOError:
		print("ERROR: Failed to open and write to file")
		return

	# loop through keys and extract corresponding sequence
	for i in keys:
		#print(i)
		cmd = ['samtools', 'faidx', args.fasta, i]
		process = subprocess.run(cmd, capture_output=True)
		#print('new entry')
		entry = process.stdout.decode("utf-8").strip()
		#print(entry)
		#seq = entry.split('\n')[1]
		#print(entry)
		#if len(seq.strip()) > 0: 
		#	print('stuff enough')
			#entry = entry.strip()
		entry = ''.join((entry, '\n'))
		if entry == '>\n':
			continue
		f.write(entry)
		#else:
		#	print(entry)
		#	print(seq)

	f.close()

	



def main(program):

	cwd = os.getcwd()

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='run_faidx', add_help=False)
	parent_parser.add_argument('-f', '--fasta', help='FASTA file of sequences (required)')
	parent_parser.add_argument('-i', '--input_file', help='Path to input file.')
	#parent_parser.add_argument('-k', '--keys', 'List of keys', type=list)
	#parent_parser.add_argument('-m', '--method', help='Dictate source of keys: file or list', choices=['roary', 'gene'])
	parent_parser.add_argument('-o', '--output_directory', default='subset_faidx', help='Prefix of output directory', type=str)
	parent_parser.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parent_parser.add_argument('-s', '--savename', default='subset', help = 'Name for output file.', type=str)
	parent_parser.add_argument('--version', action='version', version='%(prog)s 0.0.4')
	subparsers = parent_parser.add_subparsers(help='sub-command help')

	from_roary = subparsers.add_parser('roary', help='When using Roary output', parents=[parent_parser])
	#from_roary.add_argument('-n', '--num_samps', help='Number of samples to look at.', type = int)

	from_pangenome = subparsers.add_parser('pangenome', help='Separate core and accessory sequences', parents=[parent_parser])

	from_gene = subparsers.add_parser('by_gene', help='When using gene names', parents=[parent_parser])
	from_gene.add_argument('-m', '--match', default='gene', help='Select match level', choices=['exact', 'gene', 'close'])
	from_gene.add_argument('-g', '--gene', help='Gene name to search for.')

	from_file = subparsers.add_parser('by_id', help='When using a file of keys', parents=[parent_parser])
	from_table = subparsers.add_parser('table', help='When using a table of keys', parents=[parent_parser])


	args = parent_parser.parse_args()

	# TODO: Add functionality to take input as list
	try:
		ifile = File(args.input_file)
		fa = File(args.fasta)
		suffix = fa.extension
	except:
		print('ERROR: Could not find input file.')

	if program == 'pangenome':
		accKeys, coreKeys = get_keys_for_pangenome(args, ifile)

		extract_seqs(args, accKeys, suffix, 'access')
		extract_seqs(args, coreKeys, suffix, 'core')

	elif program == 'by_gene':
		keys = get_keys_from_gene_name(args, ifile)
		print(keys)
		extract_seqs(args, keys, suffix, str(args.gene))

	elif program == 'file':
		keys = get_keys_from_file(args, ifile)
		extract_seqs(args, keys, suffix, args.o)

	elif program == 'roary':
		get_keys_from_roary(args, ifile, suffix)

	elif program == 'table':
		get_keys_from_table(args, ifile, suffix)


	# TODO: Add functionality to find samtools

	#extract_seqs(args, keys)
	




if __name__ == "__main__":
	main(sys.argv[1])