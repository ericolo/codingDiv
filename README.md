# potential-garbanzo
garbanzo means pois-chiche, it was suggested by github

Exploring metagenome microdiversity to better find protein coding genes

## Usage with docker container
You need to have installed the docker app:

https://docs.docker.com/get-docker/ 


Once your docker is running, clone the repository and move into the `potential-garbanzo` directory, then run:

```diff
docker build --tag codingdiv .
```

This will take around 20 minutes the first time but will be cached if you need to rebuild it in case of an update.

You have now built your docker image named `codingdiv`.

In the cloned repository you have an example dataset with a reference genome in the `tylcv.fna` file (GenBank accession #) and 1,106 genomes exhibiting at least 90% nucleotide identity with the reference in the `blast_hits_90.fna`. 

To work with this example data, you can start by creating a custom directory on your OS, here we created a `test` directory in our home repertoire:

```bash
mkdir /home/ericolo/test/
```

Then copy the two example files into this new directory:

```bash
cp tylcv.fna /home/ericolo/test/
cp blast_hits_90.fna /home/ericolo/test/
```
You a re now ready to go ! 

**Note that your working directory can be an existing one, or even the one you just cloned, it just needs to contain the data files.**

To launch codingDiv on the example run the following command:

```diff
docker run -v /home/ericolo/test:/data codingdiv codingDiv.sh tylcv.fna blast_hits_90.fna 90 1 2 1 3 N
```

The `-v` flag tells the docker that you are working in your custom directory, which is mirrored inside the container in a directory named `/data`. This allows docker to read and write files in your custom directory on your OS. Thus you should replace `/home/ericolo/test` with the path to your own directory.

Then you are calling the `codingdiv` image, that you just built a few minutes ago with the `docker build` command, and finally comes the actual **codingDiv** script with all the positional arguments.


## Positional arguments

```bash
codingDiv.sh tylcv.fna blast_hits_90.fna 90 1 2 1 3 N
```

This is the command launched inside the container on our example dataset, the fiver numbers and the final letter are positional arguments for which you will get the detailed explanation if you run the script with no options:

```diff
docker run -v /home/ericolo/test:/data codingdiv codingDiv.sh 
```

There are 8 arguments read by codingDiv:

1- Reference genome / Studied genome (FASTA)

2- Reads or contigs to map (FASTA or FASTQ)

3- Minimal ORF size (in nucleotides) [integer]

4- Translation table number used by EMBOSS getorf - https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi [integer 1-23]

5- Minimal number of reads required to evaluate an SNP [integer]

6- Minimal % of total mapped reads required to evaluate an SNP [double]

7- Number of CPU allowed for mapping [integer]

8- Force SVG for a very large genome, over 100 kilobases [Y|N]

This last option is not recomended as it will generate a very large SVG file.
A better option would be splitting your genomes in several regions of 100Kb.


## Output & errors

If you want to jump to the results just open the `final_results` directory, but if you are interested in the intermediate files, four different directories will be created:

`final_results` will contain an SVG plot summarizing protein prediction and mapping results on genomic maps, and two TSV files which will contain the plotted information as a tab-separated table, one summarizing at the ORF level and the last one detailing the SNPs at each position of these ORFs.

`output_files` will have all the intermediate files used to generate the SVG plot and the final tables. 

`mapping_results` stores all the mapping files produced by the **BWA** aligner.

`prot_prediction` keeps track of protein and ORF prediction keeping all the fasta files. 

## Example SVG

This is the SVG plot you should expect, here is the example run on TYLCV:
![alt text](https://drive.google.com/uc?export=view&id=1tjHziIe0J7N43GqA1VKhoe1I2qfTlvsa)

## Citation

For a detailed description of the pipeline, and to cite us:

CodingDiv : visualize SNP-level microdiversity to discriminate between coding and noncoding regions.
Olo Ndela, Sasha, Yoann, Johannes Enault 2023

Laboratoire Microorganismes Genome & Environnement (LMGE)
Clermont-Auvergne University (UCA)
