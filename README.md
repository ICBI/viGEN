# viGEN
A bioinformatics pipeline for the detection and quantification of viral RNA in human NGS data

## Steps in general
![Image](https://github.com/ICBI/viGEN/blob/master/workflow.png)

## About the data

*	For this tutorial, we used a sample RNA-seq file from liver cancer from public study SRA http://www.ncbi.nlm.nih.gov/bioproject/PRJNA279878 or GEO http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65485 (also available at EBI SRA here: http://www.ebi.ac.uk/ena/data/view/SRR1946637). 
* About the dataset: Data consists of Liver cancer patients in China. The authors look for HBV-human integration sites in the dataset of 50 Liver cancer patients in China and 5 adjacent normal tissues. We downloaded the raw reads for one sample SRR1946637 and unzipped them. This sample is from an adjacent normal liver tissue. 

## Workflow steps in detail
* Alignment to human reference
  -	In this tutorial, we use tool RSEM for alignment
  -	Install RSEM as shown in https://github.com/bli25ucb/RSEM_tutorial 
  -	RSEM uses aligner Bowtie1, Bowtie2 or STAR. We chose to use Bowtie1 which is available as pre-compiled executable. For Bowtie1 installation, follow instructions in this manual: http://bowtie-bio.sourceforge.net/manual.shtml#obtaining-bowtie 
  - Test RSEM successful installation
  `./rsem-calculate-expression --version`
	`./rsem-calculate-expression --help`
  - Generate transcriptome reference for RSEM:
    
    RSEM uses a transcriptome reference for alignment. RSEM has a command called `rsem-prepare-reference` where it uses the genome fasta reference file and the GTF file to create the transcriptome reference. We have generated and provided these files for human hg19 reference genome in the [**RSEM_reference/Rsem_ref_hg19**!](https://drive.google.com/drive/folders/0B3-883ME4sP3dDF3ZllrN1JSWmM?usp=sharing) folder via google drive. 
   
  
### Citation
* NCBI SRA http://www.ncbi.nlm.nih.gov/bioproject/PRJNA279878
* EBI SRA http://www.ebi.ac.uk/ena/data/view/SRR1946637
