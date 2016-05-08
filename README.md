All code is now found in the daf16_mmp.zip
The final report is called daf16_FinalReport_2016.pdf

Some files are missing the import: import numpy as np
Any file names that are not found need the file path corrected for the directory the files are in.

The "Python_Code_Files" folder contains all the python scripts developed.
Descriptions for these can be found below.

All files should be called using "python 'filename.py'" on the command line
when in the same directory as the file.
For example: "Python ANGeneRank_final.py" 

All input files for the following should be "(geneName)_100nn.tsv" These can be
found in the "DataFiles" foulder. For the protein-protein interactions network use
the files found in the "AmmendedDataFile" foulder. Some of these files have an extra line at
the end as required for the files to run. For all files, The 1st element should be the geneID,
the 5th element should be the LogFc value, the 7th element should be the gene name the 14th 
element should be the Function GO ID's, the 15th element should be the Process GO ID's and the 
16th element should be the Component GO ID's 
For protein-protein interactions the network input files should be "(geneName)_100nn_StringNet.csv" 

ANGeneRank_Prototype is the prototype for the algorithm and final use file. 
This is intended to show that the key aspects of the algorithm run and work as expected.
This uses the Gene_ontology_annotations.txt file. 

ANGeneRank_Final is the final version of the above which outputs the ranking
for all genes. This is the file intended for normal use of the algorithm. This only 
deals with GO Ontology networks. As is this file uses the full GO network, if only one
of the three GO Ontology networks is to be used then un-comment and comment out the 
required readFile method as well as the call to runRanking in the main method.

ANGeneRank_For_PToP is the same as the above but for the protein-protein
interactions network.

Python_GeneRank_Method contains only the GeneRank method and could be used as a library.
This returns the ranking value for each gene, in the order in which they appear in the lists which
are entered.

KO_Ranking_For_All_D calculates the ranking for all genes for all values of d and then
writes the rank of the KO gene for all values of d to a specified file. This file asks
for an index for the GO terms of each gene in the input file. If all 3 networks are to be 
included use the commented put version of readFile. This uses hard coded file paths and Ko gene names 
as it is only needed for evaluation.

PToP_KO_Ranking_For_All_D does the above but for the protein-protein interactions
network.

Calculate_ROC_Over_all_40_Files calculates the roc value over all 40 KO gene files
for each value of d. Some code also exists for the KO_Ranking_For_All_D method, 
with some tweaking and use of the commented out sections it may be possible to get rid of the 
KO_Ranking_For_All_D file and use Calculate_ROC_Over_all_40_Files for all evaluation measures.
This file asks for an index for the GO terms of each gene in the input file. If all 3 networks 
are to be included use the commented out version of readFile.
This uses hard coded file paths and Ko gene names as it is only needed for evaluation.

PToP_ROC_Over_All_40_files does the same as the above but for the protein-proteins
interactions network. Much like in the above, the ROC and Ranking_For_all_D files
could potentially be put into one file.

Calculate_Evaluation_Measures takes in the output file for the 
and calculates the three evaluation methods, not including the roc score, for
all values of d. The evaluation scores for all values of d are then written to
a specified file. Use the same file for reading and writing. For this to run all 
instances of "[", "]" and "," had to be removed using find and replace in NotePad++

Rank_On_Expression_Value creates the ranking for each KO gene based only on their
expression values. Then outputs these rankings to a file. This uses hard coded file
paths and Ko gene names as it is only needed for evaluation

The "DataFiles" folder contains the expression change data and GO Ontology data for each gene.
The "AmmendedDataFiles" folder contains the same as the above but some files have an added blank line
at the end as required for the PToP files to run.
The "PToPDataFiles folder" contains the network files for the protein-protein interactions networks.
The "Graphs" folder contains 5 folders, one for each network used. In these 5 folders there is a
network graph for each of the 40 KO gene files.
The "SampleRankingOutput" folder contains sample files for the Abca1 gene for each of the 5 networks.
Using only d = 0.65. It also contains a sample ranking for the prototype.
The "KORankings" folder contains the ranking and evaluation scores for all values of d for each
of the 5 networks.

For all previous version of code go to https://github.com/D4nfire/Dissertation
which is currently private and has to be requested. Email daf16@aber.ac.uk if this
is required.
