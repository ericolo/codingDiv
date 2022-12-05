# potential-garbanzo
garbanzo means pois-chiche, it was suggested by github

Exploring metagenome microdiversity to better find protein coding genes

# Usage with docker container
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


# Positional arguments

```bash
codingDiv.sh tylcv.fna blast_hits_90.fna 90 1 2 1 3 N
```

codingDiv.sh v1.0 

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

# Output & errors

A directory `final_results` will contain both an SVG plot and a TSV file with all results.
The other created directories contain intermediate files such as mapping results, protein prediction files... 

If your final files are empty, you can check the `stdout.err` file which is a raw output file.

Expected plot:
![alt text](https://drive.google.com/uc?export=view&id=1tjHziIe0J7N43GqA1VKhoe1I2qfTlvsa)

# Citation

Cite us:

Olo Ndela & Enault 2022

Laboratoire Microorganismes Genome & Environnement (LMGE)
Clermont-Auvergne University (UCA)
