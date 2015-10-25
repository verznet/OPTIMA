OPTIMA

=== Software Information ===

OPTIMA v.f-1.3 -- 6 October 2015
Copyright (C) Davide Verzotto and Niranjan Nagarajan


OPTIMA: Sensitive and accurate whole-genome alignment of error-prone genomic
maps by combinatorial indexing and technology-agnostic statistical analysis
--
Index-based map-to-sequence alignment in large eukaryotic genomes


Resolution of complex repeat structures and rearrangements in the assembly
and analysis of large eukaryotic genomes is often aided by a combination
of high-throughput sequencing and mapping technologies (e.g. optical
restriction mapping). In particular, mapping technologies can generate
sparse maps of large DNA fragments (150 kbp--2 Mbp) and thus provide a
unique source of information for disambiguating complex rearrangements in
cancer genomes. Despite their utility, combining high-throughput
sequencing and mapping technologies has been challenging due to the lack
of efficient and freely available software for robustly aligning maps to
sequences.

Here, we introduce two new map-to-sequence alignment algorithms that
efficiently and accurately align high-throughput mapping datasets to
large, eukaryotic genomes while accounting for high error rates. In order
to do so, these methods (OPTIMA for glocal and OPTIMA-Overlap for overlap
alignment) exploit the ability to create efficient data structures that
index continuous-valued mapping data while accounting for errors. We also
introduce an approach for evaluating the significance of alignments that
avoids expensive permutation-based tests while being agnostic to
technology-dependent error rates.
Our benchmarking results suggest that OPTIMA and OPTIMA-Overlap outperform
state-of-the-art approaches in sensitivity (1.6--2X improvement) while
simultaneously being more efficient (170--200%) and precise in their
alignments (99% precision). These advantages are independent of the
quality of the data, suggesting that our indexing approach and statistical
evaluation are robust and provide improved sensitivity while guaranteeing
high precision. 



=== License and Disclaimer ===

OPTIMA
Copyright (C) Davide Verzotto and Niranjan Nagarajan

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License 2.1 as published by the Free Software Foundation.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License 2.1 along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.



=== References ===

Please cite the following papers:

Davide Verzotto*, Audrey S.M. Teo, Axel M. Hillmer, Niranjan Nagarajan*:
"OPTIMA: Sensitive and accurate whole-genome alignment of error-prone genomic
maps by combinatorial indexing and technology-agnostic statistical analysis."
Submitted to GigaScience 2015.

Davide Verzotto*, Audrey S.M. Teo, Axel M. Hillmer, Niranjan Nagarajan*:
"Index-based map-to-sequence alignment in large eukaryotic genomes."
In the Proceedings of the Fifth RECOMB Satellite Workshop on Massively
Parallel Sequencing, RECOMB-Seq 2015, Warsaw (Poland).



=== Feedback ===

For questions and suggestions about the project, please contact Davide Verzotto
(d.verzotto@gmail.com) or Niranjan Nagarajan (nagarajann@gis.a-star.edu.sg).



----



== Requirements and Dependencies ==

- Oracle Java SE Development Kit 7

- Apache Commons Math 3.2 (http://commons.apache.org/proper/commons-math/download_math.cgi,
  library needed: commons-math3-3.2.jar)

- CERN Colt 1.2.0 (https://dst.lbl.gov/ACSSoftware/colt/install.html,
  library needed: cern.jar)



=== Execution ===

Install Java 7+ SDK, and place the jar libraries "commons-math3-3.2.jar" and
"cern.jar" in the same directory containing OPTIMA (along with the source code).


Compile:

	javac  OPTIMA/Align.java


Execute:

	java  OPTIMA/Align  OpticalMaps.maps  OutputFileName  InSilicoMaps.silico 
	      [pvalue|score]  [allMaps|select]  (OPTIONAL:  firstIndex  lastIndex)
	      > OutputFileName.log  2> OutputFileName.err

Where:
- the .maps file contains genomic/optical maps in the format currently accepted
by OpGen (i.e. default output of OpGen's Argus System),
- the .silico file represents scaffold/chromosome in silico maps in SOMA v.2 format
(http://www.cbcb.umd.edu/finishing/) with the additional '>' FASTA-like symbol
to distinguish multiple scaffolds/chromosomes in the same file.

Maps and in silico files need to be generated using the same genomic or
restriction pattern(s) (e.g. with a restriction enzyme).

Options:
- You can decide between using OPTIMA's original p-value-based function ("pvalue")
or the Nagarajan et al. SOMA scoring function ("score").
- By choosing "allMaps" you will run the computation for the entire set of maps,
and therefore do not need to specify the map index range.  With "select", you
must subsequently specify the range of map indexes to be processed (note that
the range "-1 -1" will actually act as "allMaps").



Output files:

- .ok				best alignments found (significant and unique)
  					(to be further analyzed with the q-value FDR estimation);
  					
- .notOk			non-significant and/or non-unique best alignments found;

- .notFound			maps with no alignment found;

- .discarded		discarded maps (thresholds currently set to minimum 10
					fragments and 50 kbp -- for optical maps only);
					
- .otherSolutions	all candidate alignment results in the format: Z-score (to
                    be converted into p-value), map size.  (These values can be
                    used later in the q-value analysis);
  					
- other supporting files, output log file and error log file.



Explanation of the main output file in .ok format:

	Map_index	Map_name	Significant	Unique	Map_length	Map_size	P-value	Map_on_reverse_strand	Scaffold_first_location	Length_in_scaffold	Used_size_comparison_map/scaffold	Map_used_size	Scaffold_used_size	Length_comparison_scaffold/map	Non-end_OM_frags_found/map_length	Matches	Missing_cuts	False_cuts	Missing_frags	Perc_missing_cuts	Ratio_false_cuts	Ratio_small_frags_found/total_small_frags	WHT_chi_square	Chi_square_per_match	Chi_square	DP_score	P_missing_cuts	P_false_cuts	P_WHT_chi_square	Z-score	Z-score_error_cuts	Z-score_no_of_matches	Z-score_WHT_chi_square	AFS_map	AFS_map_non-end_frags	AFS_used_scaffold	Abs_fragment_sizing_error	Relative_fragment_sizing_error	Map_length2	Map_length_non-end_frags	Used_scaffold_length	Used_scaffold_length2	Small_frags	Used_scaffold_size-missing_frags_found	Map_seed_location	Map_seed_first_size	Scaffold_seed_first_size	Abs_WHT_comparison	Abs_WHT_squared_comparison	Relative_WHT_comparison	Seeds_in_same_region	Alignment	Other_seeds	Very_small_frags_not_detected_as_missing	Very_small_frags_detected_as_missing	Second_best_p-value	Second_best_location



Indexes/coordinates start from 0 and are all-inclusive.  Orientation and
coordinates are based on the input scaffold/chromosome in silico maps
(concatenated to each other with a virtual separation fragment).



In order to get the left-most location of an alignment (i.e. 'Scaffold_first_location',
say L), please run the following awk UNIX/Linux script by setting the location L
(e.g., L = 259466) and the in silico file name (e.g. Human_reference_hg19-KpnI.silico):

	awk 'BEGIN{inputFragLoc=L; inputFragLoc += 1; allChrFragsSum=0; nextSeq =0; remainingInputFrags=0; lastChrFragSize=0} {header=substr($1,1,1); if (header==">"){allChrFragsSum += ($3 + 2); if (allChrFragsSum == inputFragLoc) {print "ERROR: End of scaffold/chromosome!!!"; exit;} if (allChrFragsSum >= inputFragLoc) {print substr($1,2,length($1)-1); remainingInputFrags=(inputFragLoc - (allChrFragsSum - $3 - 2)); lastChrFragSize=$2; nextSeq++;}} else if (nextSeq == 1) {len=split($0,fragSizesArray," "); for (countFrags=1; countFrags <= len; countFrags++) {if (countFrags == remainingInputFrags) {print (countFrags >= 2 ? fragSizesArray[countFrags] : "0"); exit;}} print lastChrFragSize; exit; }}' ./Human_reference_hg19-KpnI.silico     # consider as number of scaffold/chromosome fragments:  no. of cuts + 1 + 1 (the latter refers a separation fragment)   # delta-like indexes



A q-value false discovery rate (FDR) estimation can be performed by using e.g.
the following Bioconductor Q-value v.1.40 software package in R (input: p-values
from the output file .otherSolutions):

	Storey, J.D., Tibshirani, R.: Statistical significance for genomewide
	studies.  Proceedings of the National Academy of Science 100, 9440-9445
	(2003).
	Bioconductor URL:
	https://www.bioconductor.org/packages/release/bioc/html/qvalue.html



=== OPTIMA-Overlap ===

Overlap alignment can be achieved by splitting the input optical maps into
sliding windows of the desired length, and separately aligning all of them to
the input scaffold/chromosome in silico maps with OPTIMA.
