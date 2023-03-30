# nanotools
A repository of tools for analyzing nanopore data

# nanoStats
nanoStats.py calculates the read length and average quality score for each read in a Fastq file.

Fastq files contain four lines for each read. Line 1 specifies the read ID and any meta data. Line 2 contains the nucleotide sequence. Line 3 contains a “+” to indicate a quality score identifier. Finally, line 4 encodes the quality sequence. The quality sequence consists of ASCII characters from 33 to 126. In Illumina reads, these characters have an offset of 33 and represent Phred scores of 0 to 93. For example, the character ‘!’ has an ASCII value of 33 and represents the worst Phred score of 0 while the character ‘]’ has an ASCII value of 126 and represented the highest possible Phred score of 93. Though not explicity stated, Nanopore appears to use the same system for calculating per base quality scores.

nanoStats.py provides two outputs to the user. 
* Table file (.tsv format): A table where each record contains the read ID, read length, and average quality score. 
* Report file: Provides basic statistics on the length and quality metrics. 
    * Total read number
    * Total number of bases
    * Length: mean, median, minimum, maximum, and standard deviation
    * Average read quality: mean, median, minimum, maximum, and standard deviation

Finally, nanoStats.py offers the option to filter the reads based on quality so that only reads with an average quality score at or above this threshold are displayed in the table and assessed for the report.

## Requirements
* Python and Python standard libraries

## Usage
```
# Example: Basic usage
python nanoStats.py -i input.fastq

# Example: Send output to another directory
python nanoStats.py -i input.fastq -p /other/path/

# Example: Change name of output
python nanoStats.py -i input.fastq -s newname

# Example: Calculate stats on reads with average qscore > 10
python nanoStats.py -i input.fastq -m 10
```



# mega_bed2gff
mega_bed2gff.py is a tool to filter and reformat Megalodon modification calling output.

The program Megalodon examines the raw ONT reads at specific motifs in a given genome and determines the probability that each site is modified. It outputs these probabilities as a confidence score for each position in the genome in .bed format. However, loading this track into a genome viewer gives the visual appearance that every site is modified. For example, if using a model trained to identify 5mC and 5hmC base modifications at CG motifs, then Megalodon will generate the files modified_bases.5mC.bed and modified_bases.5hmC.bed. These files will have a confidence score for 5mC or 5hmC, respectively, modifications for each instance of CG in the genome sequence. 

mega_bed2gff.py was developed to (1) filter the modification calls based on a specified confidence threshold and (2) output the result in .gff format. As input, mega_bed2gff.py only requires a .bed file generated by Megalodon and it will output a .gff file. The .gff files are compatible for downstream analysis and visualization in the genome browser IGV.

By default, mega_bed2gff.py uses a 80% confidence score when filtering (per ONT recommendation), but the threshold can be changed using the `-t` flag. 

## Requirements
Python and Python standard libraries

## Usage
```
# Example: Basic usage
Python mega_bed2gff.py -i input.bed

# Example: Send output to another directory
python mega_bed2gff.py -i input.fastq -p /other/path/

# Example: Change name of output
python mega_bed2gff.py -i input.fastq -s newname

# Example: Filter for modifications with a 90% confidence
python mega_bed2gff.py -i input.fastq -t 90
```
