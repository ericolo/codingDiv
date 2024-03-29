#!/usr/bin/env python3

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from collections import defaultdict
import re
from dna_features_viewer import GraphicFeature, GraphicRecord,CircularGraphicRecord
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import numpy as np

genomes=[sys.argv[1]]

sizes=[int(sys.argv[2])]

orf_file=sys.argv[3]

if sizes[0]<=20000:
	plot_size=12
elif sizes[0]<=50000:
	plot_size=22
elif sizes[0]>50000:
	plot_size=42

for i in range(len(genomes)):

	genome=genomes[i]
	size=sizes[i]
	print(genome)

	fig, (ax, ax2) = plt.subplots(
	    2, 1, figsize=(plot_size, 3), sharex=True, gridspec_kw={"height_ratios": [4, 1]}
	)


	dprot_coord=defaultdict(list)
	dprot_clu=defaultdict(str)
	dprot_color=defaultdict(str)

	with open(orf_file,"r") as f1:
		for li in f1:
			li=li.rstrip()
			lp=li.split()

			if li.startswith(">") and re.sub("_[0-9]+$","",lp[0].lstrip(">"))==genome and "REVERSE" in li:
				dprot_coord[lp[0].lstrip(">")]=[int(lp[3].rstrip("]")),int(lp[1].lstrip("["))]



	##########################


	#To have each prot in its frame
	dprot_frame=defaultdict(str)

	for prot in dprot_coord:
		if size-(dprot_coord[prot][1]-1) in list(range(1,size,3)):
			dprot_frame[prot]="1.4"

		elif size-(dprot_coord[prot][1]-1) in list(range(2,size,3)):
			dprot_frame[prot]="0.7"

		elif size-(dprot_coord[prot][1]-1) in list(range(3,size,3)):
			dprot_frame[prot]="0"	


	#############################

	features=[]
	sjr=[]

	#Every gene needs a name, unless you remove the label tag
	for prot in dprot_coord:
		sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][1], strand=-1, color="white", label=prot, linecolor="green"))


	record = GraphicRecord(sequence_length=size, features=features)


	genome_plot, _=record.plot(figure_width=plot_size,with_ruler=True,figure_height=2)

	for elt in sjr:
		sjr_plot=record.plot_feature(genome_plot,elt,level=float(dprot_frame[elt.label]))

	genome_plot.figure.savefig("./"+genome+'_getorf_neg.svg', bbox_inches='tight', dpi=300)