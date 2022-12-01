#!/usr/bin/env python3

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

depth_file=sys.argv[4]

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
		next(f1)
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
			print(print)

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

		
		#print(prot, dprot_pnps[prot], color)

		sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][1], strand=-1, color=color, label=prot ))

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

		if prot in dprot_fb and dprot_fb[prot]!=0 and dprot_coord[prot][1]-dprot_coord[prot][0] >= 100:
			sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][0]+50, strand=0, color=color, label=prot ))
		elif dprot_fb[prot]!=0 and dprot_coord[prot][1]-dprot_coord[prot][0] < 100:
			sjr.append(GraphicFeature(start=dprot_coord[prot][0], end=dprot_coord[prot][0]+25, strand=0, color=color, label=prot ))

	record = GraphicRecord(sequence_length=size, features=features)


	genome_plot=record.plot(ax=ax,figure_width=plot_size, figure_height=2,with_ruler=True)

	for elt in sjr:
		sjr_plot=record.plot_feature(ax,elt,level=float(dprot_frame[elt.label]))

	ax.set_xlabel("pNeg/pS ratio | negative strand",loc="left", weight='bold', color="black",size=12)
	ax.xaxis.set_label_position('top')

	dgenome_gc=defaultdict(list)

	dpos_gc=defaultdict(int)

	with open(depth_file) as f4:
	        for li in f4:
	                li=li.rstrip()
	                lp=li.split("\t")

	                dpos_gc[int(lp[1])]=int(lp[-1])

	for pos in range(1,size+1):
		if pos in dpos_gc:
			dgenome_gc[genome].append(dpos_gc[pos])

		else:
			dgenome_gc[genome].append(0) 

	ax2.plot(np.arange(0,size), dgenome_gc[genome], alpha=1,color="green")

	ax2.set_xlabel("#mapped reads per position (max="+str(max(dgenome_gc[genome]))+")")

	#ax2.yaxis.set_label_position("right")
	#ax2.yaxis.tick_right()

	ax2.axes.get_yaxis().set_ticks([])

	fig.savefig("./"+genome+'_neg_strand_pnps.svg', bbox_inches='tight', dpi=300)



#######################################################
#### pNpS legend
fig, ax = plt.subplots(figsize=(4.5, 0.5))
fig.subplots_adjust(bottom=0.5)

cmap = mpl.colors.ListedColormap(["#2c4099", "#5f62ae",
                                  "#8887c2", "#b0add6", "#d7d5eb", "#ffffff", "#eecac6", "#d89791", "#be645e", "#a02d30"])
cmap.set_over('red')
cmap.set_under('blue')

bounds = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8,2]

norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb3 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                norm=norm,
                                boundaries=bounds,
                                #extend='max',
                                #extendfrac='auto',
                                ticks=bounds,
                                spacing='uniform',
                                orientation='horizontal')
cb3.set_label('pNeg/pS gradient',labelpad=-26, y=1.5, x=-0.22, weight='bold')
fig.savefig("pnps_legend.svg", bbox_inches='tight', dpi=300)