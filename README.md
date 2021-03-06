# AGILE: Assembled Genome minIng pipeLinE 
Methodologies for genome mining and annotation
---------
Synopsis
---------
The goal of AGILE is to mine and annotate genes from a target genome  using a set of reference genes from a closely related taxon. This script combines a number of bioinformatics tools, in addition to a number of extra steps to annotate genomic regions while overcoming limitations in genome annotation software. AGILE is wirtten in Perl for command line execution on a linux operating system.
 
------------------------
Pre-constructed pipeline
------------------------
We have 2 pre-constructed versions of AGILE available, with all necessary dependancies, binaries and files installed already. One is a Virtual Box disk image (VDI), and the other is a Singularity image. 
Both can be downloaded here: https://figshare.com/s/a0004bf93dc43484b0c0

Virtual Box (password: guest) can be downloaded here: https://www.virtualbox.org/

Singularity can be downloaded here: https://singularity.lbl.gov/install-linux

-------------
Prerequisites
-------------
AGILE combines binaries from BLAST, MAKER, AUGUSTUS, SNAP, EMBOSS and Framebot
  
The following are locations for acquiring each binary (sites active as of writing this readme):

BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ 

MAKER: http://www.yandell-lab.org/software/maker.html

EMBOSS: ftp://emboss.open-bio.org/pub/EMBOSS/   

Framebot: https://github.com/rdpstaff/

Note that MAKER contains binaries for AUGUSTUS and SNAP as part of its installation. MAKER also contains BLAST binaries, and therefore contains almost all external software required. Only the matcher binary is required from EMBOSS.

==MAKER Instillation==

MAKER can be found at the above address. This pipeline has been written using version 2.31-8. 
Once MAKER has been unpacked, go to the /maker/src/ directory. In here,type 'perl Buil.PL'. This will install a number of dependancies for maker to run (you can choose No for installing MPI, this will not affect maker. From here, you can install all depednables by typing:
'./Build installdeps' and './Build install exes'. You can see thet build status by typing './Build status'. This will tell you what files are needed for maker to run.

Using './Build', you can make the binaries for SNAP, augustus and exonerate with './Build snap', './Build augustus' and './Build exonerate'. The download locations for these binaries can be found (or modified) in the 'downloads' file in the src directory in /maker/. 

==EMBOSS Instillation==

To install emboss once downloaded, type the following commands:
'./configure --prefix=/usr/local/emboss --without-x'
'make'
'sudo make install'

Once installed, you can place the "matcher" binary into the maker/bin/ directory or add the path to matcher's location to the ctrlfile 


==FrameBot Instillation==
you can downlaod RDPTools from GitHub. 
Once downloaded and unzipped (or cloned from github), you may have to install 'ant' via the command 'sudo apt install ant'
When you have ant installed, you can use 'ant jar' to install. 

Alternatively, you can download pre-compiled versions in our github directory (lib.zip and Framebotfiles.zip) and move them into the maker/bin/ directory

Once installed, you can copy all "jar" files and the directory 'lib/' to /maker/bin/ or add the path to the ctrl file
   
----------------------
Installation/File list
----------------------
AGILE is written in Perl, and has been tested using version 5.10.1. The execution of agile does not require additional, non-standard perl libraries or BioPerl. 

Installation of maker will result in a maker/bin/ directory. The AGILE script, and all associated files should be placed within this 'maker/bin/' directory. Additionally, all blast binaries and EMBOSS matcher can be placed within this directory, however their location can be specified in the ctrlfile..

Files required for the execution of this script are:

agile.pl: The perl file containing the code

ctrlfile:Control file containing modifieable parameters used by agile.pl

maker_bopts.ctl, maker_opts.ctl: standard maker control files, modified for agile.pl


-----------------------
File format/Preparation
-----------------------
=>Gene list file: The gene list, named "target_genes" must contain a list of genes that the user wishes to identify and mine/annotate in the target genome. See "target_genes" in the example directory. The following is an example of the "target_genes" list:
       
       TPP1
       
       ELOVL6
       
       TNMD
       
       SCYL3
       
       FGR

=>The file containing reference sequences must be in fasta format, with ">" delimiting the gene information. This fasta header must contain the gene name in the format "(Gene name)", where 'Gene name' matches the gene name specified in the aforementioned gene list file: 
      
      ">NM_000391.3 Homo sapiens tripeptidyl peptidase I (TPP1), mRNA" would match "TPP1" in the gene list

=>The genome file must contain fasta headers, however these headers will be modifed during the "shearing" step, for ease of parsing downstream. 

-----------------------
Parameters and Examples
-----------------------
AGILE can be executed by calling Perl: "perl agile.pl". Alternatively, the script can be made executable by typing "chmod u+x "agile.pl" and calling it as "./agile.pl". 

The ctrlfile contains a number of parameters and locations for executables. The current ctrlfile expects a file of reference gene names called "target_genes",which will be used to pull reference sequences from "Homo_sapiens_longest_transcript.fa" to mine "GCF_000001405.35_GRCh38.p9_genomic.fna" and output results to a directory called "HumanvHumanTest":

Target_genes:"target_genes";
GenomeFile:"GCF_000001405.35_GRCh38.p9_genomic.fna";
SisterTaxa:"Homo_sapiens_longest_transcript.fa";
TestName:"HumanvHumanTest";  

The parameters used are:

=>Overhang:"0"; The number of bp overlapping between sheared contig files.
=>ShearLength:"4000"; The length of each sheared contig. 
=>OverlapDiff:"20";The minimum number of overlapping bp between matching regions to be considered putative paralogs .
=>Evalue:"1e-5"; Blast hit threshold.
=>Score:"40";#Minimum score to be observed in BLAST results before a result can be considered.
=>Ref_Coverage:"70"; The minimum percentage of coverage between the  annotated/mined gene and initial reference gene, used to remove poor quality results.  

The example directory contains the list of genes required from the human genome in "target_genes". The reference genes used to mine the human genome are contained in "Homo_sapiens_longest_transcript.fa" (part of benchmarking data). The final sequences output by this pipeline are contained in "HumanvHuman_Final_Genes.fa". These data were mined from the human genome, version GRCh38.p9.   


-------
Outputs
-------
The script will create a directory whose name is specified in the ctrlfile. All outputs, including maker data, alignmnets and blast results files will be transfered to this directory. 
A sheared copy of the genome file is created in the initial steps of the workflow. This can be removed after completion. 

Final gene sequences can be found in the "XXX_Final_Genes.fa" file. Each sequence contains Identity shared with the reference sequence, The percentage of the reference covered by the final sequence and the percentage of the maker output used. A good output will have a high percentage of coverage and identity, ideally representing 100% of the original maker gene: 

">HumanvHumanTest_TPP1_IdentwithRefspec_99.9_CoverageofRef_99.17_percentage_of_original_maker_gene_100.00"

Note that identity and coverage will decrease using more divergent taxa, with varying number of exons and sequence diversity


-------
License
-------
See file 'LICENSE.md' for license. If you found our software useful in your analysis, we ask that you cite us as:

