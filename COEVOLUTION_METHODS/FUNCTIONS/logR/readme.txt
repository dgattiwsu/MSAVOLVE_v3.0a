This directory contains all the code that is necessary to calculate coevolution between pairs of residues
in a multiple sequence alignment (described in Burger and van Nimwegen, PLoS Comp Biol 2010). This version of the code does NOT allow the incorporation of prior information.
The directory contains a (very simple) makefile to compile the C++ file that does most of the computations. So in order to be able to use
the code, type make and then use the perl wrapper runContactPredictions.pl as follows:
perl runContactPredictions.pl multifasta_file
The multifasta_file should end in .fasta. The directory contains one example of such a file: FKB_binding.fasta
If you execute the wrapper with FKB_binding.fasta,  it automatically generates a number of different files:
FKB_binding.aln: different format of the alignment, used by the C++ program
FKB_binding_ungapped.aln: alignment with all columns removed that contain more than 20 percent gaps (this is set by the variable threshold in runContactPredictions.pl and can be easily changed there)
FKB_binding.mapping: contains the mapping from gapped to ungapped alignment
FKB_binding_ungapped.post: contains posterior probabilities of direct coevolution between pairs of residues in the coordinate system of the ungapped alignment
FKB_binding.post: This is the main file of interest: It contains posterior probabilities of direct coevolution between pairs of residues in the original coordinates system of the input multifasta file. The format of this file is:
position1 position2 log(posterior_probability)
Note that the positions are counted starting at 0.

for questions or comments, contact me: lukas.burger@fmi.ch
