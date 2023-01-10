## RNA-seq (bulk) of sleep deprivation in wildtype mice at 4 different ages
This repository contains the code and analyses for the analysis of bulk RNA-seq analysis of sleep deprivation in wildtype mutant mice at different developmental stages.
Paper: "Ontogenesis of the molecular response to sleep loss", 2023
Goal: investigate transcriptional changes in the prefrontal cortex of mice after sleep deprivation at different developmental ages. 

### Authors:
Christine Muheim (christine.muheim@wsu.edu)
Katie Ford (kaitlyn.ford@wsu.edu)
Lucia Peixoto (lucia.peixoto@wsu.edu)

### Data
bulk RNA-seq from prefrontal cortex, the FASTQ files for WT mice at postnatal day 16, 24, 30 and 90 can be found in GEO project GSCE 211301. 
Download GEO metadata.

### Quantification
To quantify the raw fastq files, we used salmon quant (v1.8, Patro et al., 2017) using an index file build with GENCODE. Below are the steps to build the index and how to run Salmon.

download GENCODE files
Run the shell script 01_quantification/download-gencode-files.R. This downloads GENCODE files and creates the files necessary for salmon. Reference files were obtained for mouse release version M28.
Files needed:
GRCm38.primary_assembly.genome.fa.gz - nucleotide (DNA) sequences of the GRCm38 primary genome assembly.
gencode.vM25.transcripts.fa.gz - nucleotide (DNA) sequences of all transcripts on reference chromosomes.
gencode.vM25.annotation.gtf.gz - gene annotation on the reference chromosomes (i.e. for humans, these are chromosomes 1 to 22, X, and Y), i.e. locations of genes and other information about the genes, gene structure
Gene transfer format (GTF) is a file format used to hold information about gene structure. It is a tab-delimited text format based on the general feature format (GFF), but contains some additional conventions specific to gene information.

The specific locations of where the files were pulled from are from here:
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/GRCm39.primary_assembly.genome.fa.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.transcripts.fa.gz
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M28/gencode.vM28.annotation.gtf.gz

Next create decoys for salmon index
The decoy sequence is going to be the whole genome sequence (GRCm39.primary_assembly.genome.fa.gz). You can read more about decoy sequences in Salmon below:
https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
https://github.com/COMBINE-lab/SalmonTools/blob/master/README.md
Source for code: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

In general:
After creating the fasta file with transcript and intron sequences as above, we can index it using Salmon. Here, we add the full genome sequence as decoy sequences (see Srivastava et al., 2019 for more details). This is recommended in order to avoid reads truly originating from intergenic regions being assigned to a suboptimal transcriptome location. However, the effect of including decoys is typically smaller when both transcripts and introns are being quantified than for ‘regular’ gene expression quantification, since in the former case a larger fraction of the genome is already covered by the features of interest.

We use decoys for bulk RNA-seq analyses.
To use a decoy, we need to create two files:
decoys.txt is the names of the genome targets (decoys), will be used in the -d parameter in build-index-salmon.sh
gentrome_transcripts_mouse.fa.gz is a concatenated FASTA transcriptome, will be used in the -t parameter in build-index-salmon.sh (see below). Note that you need to recreate this once per time you set up salmon to for quantification.
These two files are created in the 01_quantification/create-decoys-salmon.sh and will both be used 01_quantification/build-index-salmon.sh file to build the salmon index (see next section).

#### Install and build salmon index

This part will have to be done for each user, we suggest however to use a high performance computing cluster. 
To install salmon v1.8.0, use terminal:

cd /users/christinemuheim/data/
wget https://github.com/COMBINE-lab/salmon/releases/tag/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz
tar xzvf salmon-1.38.0_linux_x86_64.tar.gz
rm salmon-1.8.0_linux_x86_64.tar.gz
Also, make sure this is in your .bash_profile file

PATH=$PATH:/users/christinemuheim/data/salmon-latest_linux_x86_64/bin
You can check to make sure salmon has been upgraded correctly using salmon -h inside terminal (or help with specific parts of using salmon using e.g. salmon index -h for help with the index step).

The -t argument is the input transcripts file.
The -i argument is the index file to create.
The -d argument is the decoy sequence.
The --keepDuplicates argument forces all duplicate transcripts (for example, multiple unspliced transcript of the same gene that are identical for example) that appear in the input will be retained and quantified separately. If you keep the duplicates they will be assigned identical expression levels since salmon can’t tell them apart. When you aggregate on the gene level, this will not make a difference any more. Therefore, I do not keep the duplicates as we are interested in gene level aggregation.
The --gencode flag will handle the composite fasta headers in GENCODE transcript fasta files and split the transcript name at the first '|' character.
The --threads argument says how many threads to use when building the index.
The salmon index is built using the 01_quantification/build-index-salmon.sh file (used 4 cores). The bulk RNA-seq index uses the decoys.txt and is built from the combined FASTA file (gentrome_transcripts_mouse.fa.gz). The snRNA-seq index uses the GRCm39.primary_assembly.genome.chrnames.txt and is built from the combined FASTA file.

#### Running salmon quant
We use the index (gencode.vM28-salmon-index-v1.0.0-mouse-withdecoy) created by Build-SalmonIndex.sh for the bulk RNA-seq analysis. See the 01_quantification/SalmonQuant.sh file. This script quantifies reads at the transcript level, which can be summarized later at the gene level using the tximeta::summarizeToGene() function.

#### Create SummarizedExperiment object
Here we import the quant.sf files produced by salmon quant into R/Bioconductor using the tximeta package with the 01-quantification/Tximeta.R script. We save the SummarizedExperiment object at the transcript level and gene level as a .RDS file: data/se_mouse_sleep_complete.rds (transcript level) and data/gse_mouse_sleep_complete.rds (gene level).
For the gene level, we remove inferential repeats because the uncertainty information is not used as the  suggested by Mike Love (https://mikelove.github.io/counts-model/quantification.html).

Helpful vignette: https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html

## RUVs

The vignette for RUVs can be found here: https://bioconductor.org/packages/release/bioc/manuals/RUVSeq/man/RUVSeq.pdf. 
We run RUVs with biological replicas and a negative gene list (02-analysis/RUV/Data/NegControl_Genes.txt) and tested the efficacy with a list of genes (02-analysis/RUV/Data/PosControl_Genes.txt) known to be changed by sleep deprivation.

## Visualization
After RUVs we use Venn-Euler.R to get venn (or Euler) diagrams and the list of interest for enrichment analysis (DAVID).
The output from DAVID was visualized using BubblePlots_EnrichedTerms.R. HeatMaps_ExpressionAccrossAge.R was run on a selected set of genes.

