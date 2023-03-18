# nanotools
A repository of tools for analyzing nanopore data

# nanoStats
nanoStats.py calculates the read length and average quality score for each read in a Fastq file.

Fastq files contain four lines for each read. Line 1 specifies the read ID and any meta data. Line 2 contains the nucleotide sequence. Line 3 contains a “+” to indicate a quality score identifier. Finally, line 4 encodes the quality sequence. The quality sequence consists of ASCII characters from 33 to 126. In Illumina reads, these characters have an offset of 33 and represent Phred scores of 0 to 93. For example, the character ‘!’ has an ASCII value of 33 and represents the worst Phred score of 0 while the character ‘]’ has an ASCII value of 126 and represented the highest possible Phred score of 93. Though not explicity stated, Nanopore appears to use the same system for calculating per base quality scores.

nanoStats.py provides two outputs to the user. 
* Table file (.tsv format): A table where each record contains the read ID, read length, and average quality score. 
* Report file: Provides basic statistics on the length and quality metrics. 
** Total read number
** Length: mean, median, minimum, maximum, standard deviation, and total bases
** Average read quality: mean, median, minimum, maximum, and standard deviation

Finally, nanoStats.py offers the option to filter the reads based on quality so that only reads with an average quality score at or above this threshold are displayed in the table and assessed for the report.
