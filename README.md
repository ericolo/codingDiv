# potential-garbanzo
garbanzo means pois-chiche, it was suggested by github

Exploring metagenome microdiversity to better find protein coding genes

# Usage with docker container
You need to have installed the docker app:

https://docs.docker.com/get-docker/ 


Once your docker ins running, clone the repository and move into the `potential-garbanzo` directory, then run:

```diff
- docker build --tag codingdiv .
`

This will take around 20 minutes the first time but will be cached if you need to rebuild it in case of an update.

You have now built your docker image named `codingdiv`, run the example to test if it's working:

```diff
- docker run -v /Users/ONE/Downloads/codingdiv:/data codingdiv codingDiv.sh tylcv.fna blast_hits_90.fna 90 1 2 1 3 N
```

The -v option allows you to read and write files on your OS from your container, just specify your working directory, mine was `/Users/ONE/Downloads/codingdiv` and add `:/data` as it is the working directory in the container. 

Normally you should get three files in the `final_results` directory, two tables and an SVG plot as the one showed below.


# Required tools and versions
**Python v3.6.4 & R v3.6.2**

## Protein predictors
**prodigal v2.6.2**
https://github.com/hyattpd/Prodigal

**phanotate v1.5.0**
https://github.com/deprekate/PHANOTATE

**getorf from the EMBOSS suite v6.6.0.0**
http://emboss.open-bio.org/html/adm/ch01s01.html 

## Python packages
**Biopython v1.79**
https://biopython.org/wiki/Download 

**matplotlib v3.3.4**
https://matplotlib.org/stable/users/installing/index.html

**dna_features_viewer v3.1.1**
https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer

**svg_stack v0.1.0**
https://github.com/astraw/svg_stack

## R packages
**tidyverse v1.3.0**
https://www.tidyverse.org/packages/#:~:text=Installation%20and%20use,in%20your%20current%20R%20session.

**ggplot2 v3.3.5**
https://ggplot2.tidyverse.org/

## Software needed for mapping reads and SNP calling
**samtools v1.10**
http://www.sthda.com/english/wiki/install-samtools-on-unix-system#google_vignette

Keep only samtools from this repo:
https://sourceforge.net/projects/samtools/files/samtools/1.10/

**bcftools v1.14**

You can follow the same process as done for samtools, and you can find the right version here:
https://sourceforge.net/projects/samtools/files/samtools/1.14/

**bwa v0.7.17-r1198-dirty**

Same installation process as for samtools
https://sourceforge.net/projects/bio-bwa/files/

Once you have all needed software, download and add the garbanzo directory to your path (in your bash config file) 
```bash
export PATH=$PATH:"~/your_own_dir/garbanzo"
```

Also check that the correct versions of the dependencies are the ones called by default, confusion can happen if you have several versions of the same tool.

You are now ready to use it!

# Usage

Example command: 

```bash
pipeline.sh reference_genome.fasta reads_or_contigs.fasta 60 1 2 1 3 N
```

pipeline.sh v1.0 

Positional arguments: 
1- Reference genome / Studied genome (FASTA)

2- Reads or contigs to map (FASTA or FASTQ)

3- Minimal ORF size (in nucleotides) [integer]

4- Translation table number used by EMBOSS getorf - https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi [integer 1-23]

5- Minimal number of reads required to evaluate an SNP [integer]

6- Minimal % of total mapped reads required to evaluate an SNP [double]

7- Number of CPU allowed for mapping [integer]

8- Force SVG for a very large genome, over 100 kilobases [Y|N]

This last option is not recomended as it will generate a very large SVG file.
A better option would be splitting your genomes in several regions.

Cite us:

blablabla et al. 2036
Laboratoire Microorganismes Genome & Environnement (LMGE)
Clermont-Auvergne University (UCA)

# Output & errors

A directory `final_results` will contain both an SVG plot and a TSV file with all results.
The other created directories contain intermediate files such as mapping results, protein prediction files... 

If your final files are empty, you can check the `stdout.err` file which is a raw output file.

Expected plot:
![alt text](https://drive.google.com/uc?export=view&id=1tjHziIe0J7N43GqA1VKhoe1I2qfTlvsa)
