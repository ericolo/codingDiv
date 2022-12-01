#!/usr/bin/env python3

from collections import defaultdict
import re
from dna_features_viewer import GraphicFeature, GraphicRecord,CircularGraphicRecord
from Bio import SeqIO
import matplotlib.pyplot as plt
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

	dprot_coord=defaultdict(list)
	dprot_clu=defaultdict(str)
	dprot_color=defaultdict(str)

	with open(orf_file,"r") as f1:
		next(f1)
		for li in f1:
			li=li.rstrip()
			lp=li.split()

			if li.startswith(">") and re.sub("_[0-9]+$","",lp[0].lstrip(">"))==genome and "REVERSE" not in li:
				dprot_coord[lp[0].lstrip(">")]=[int(lp[1].lstrip("[")),int(lp[3].rstrip("]"))]


	##########################


	#To have each prot in its frame
	dprot_frame=defaultdict(str)

	for prot in dprot_coord:
		if dprot_coord[prot][0] in list(range(1,size,3)):
			dprot_frame[prot]="0"

		elif dprot_coord[prot][0] in list(range(2,size,3)):
			dprot_frame[prot]="0.7"

		elif dprot_coord[prot][0] in list(range(3,size,3)):
			dprot_frame[prot]="1.4"	


	#############################

	#pNpS coloring

	#PNPS data
	 
	dprot_pnps=defaultdict(float)

	dprot_fb=defaultdict(float)

	with open("summary_table.tsv","r") as f2:
		next(f2)
		for li in f2:
			li=li.rstrip()
			lp=li.split("\t")

			if lp[-3]!="NA":
				dprot_pnps[lp[0].rstrip("+-")]=float(lp[-3])
			elif lp[-3]=="NA":
				dprot_pnps[lp[-3]]="NA"

			#For the false beginnings (50 nucleotides tested if >100, else 25)
			if lp[-2]=="yes":
				dprot_fb[lp[0].rstrip("+-")]=float(lp[-1])

	for prot in dprot_coord:
		if prot not in dprot_pnps:
			dprot_pnps[prot]="NA"

	features=[]

	#separate zone for sjr
	sjr=[]

	#Every gene needs a name, unless you remove the label tag

	for prot in dprot_coord:

		if dprot_pnps[prot]=="NA":
			color="gray"

		elif dprot_pnps[prot]<=0.2:
			color="#2c4099"

		elif dprot_pnps[prot]<=0.4:
			color="#5f62ae"

		elif dprot_pnps[prot]<=0.6:
			color="#8887c2"

		elif dprot_pnps[prot]<=0.8:
			color="#b0add6"

		elif dprot_pnps[prot]<1:
			color="#d7d5eb"

		elif dprot_pnps[prot]<=1.2:
			color="#ffffff"

		elif dprot_pnps[prot]<=1.4:
			color="#eecac6"

		elif dprot_pnps[prot]<=1.6:
			color="#d89791"

		elif dprot_pnps[prot]<=1.8:
			color="#be645e"

		elif dprot_pnps[prot]>1.8:
			color="#a02d30"

		
		sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][1], strand=1, color=color, label=prot ))

		if dprot_fb[prot]=="NA":
			color="gray"

		elif dprot_fb[prot]<=0.2:
			color="#2c4099"

		elif dprot_fb[prot]<=0.4:
			color="#5f62ae"

		elif dprot_fb[prot]<=0.6:
			color="#8887c2"

		elif dprot_fb[prot]<=0.8:
			color="#b0add6"

		elif dprot_fb[prot]<1:
			color="#d7d5eb"

		elif dprot_fb[prot]<=1.2:
			color="#ffffff"

		elif dprot_fb[prot]<=1.4:
			color="#eecac6"

		elif dprot_fb[prot]<=1.6:
			color="#d89791"

		elif dprot_fb[prot]<=1.8:
			color="#be645e"

		elif dprot_fb[prot]>1.8:
			color="#a02d30"

		if prot in dprot_fb and dprot_fb[prot]!=0 and dprot_coord[prot][1]-dprot_coord[prot][0] > 100:
			sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][0]+50, strand=0, color=color, label=prot ))
		else:
			sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][0]+25, strand=0, color=color, label=prot ))


	record = GraphicRecord(sequence_length=size, features=features)


	genome_plot, _=record.plot(figure_width=plot_size,figure_height=2,with_ruler=True)
	
	
	for elt in sjr:
		sjr_plot=record.plot_feature(genome_plot,elt,level=float(dprot_frame[elt.label]))

	genome_plot.set_xlabel("pNeg/pS ratio | positive strand",loc="left", weight='bold', color="black",size=12)
	genome_plot.xaxis.set_label_position('top')


	genome_plot.figure.savefig("./"+genome+'_pnps.svg', bbox_inches='tight', dpi=300)

