#!/usr/bin/env python3

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'

from collections import defaultdict
import re
import sys


#Script trying to find the mirror effect of SNPs between reading frames

genomes=[sys.argv[1]]

sizes=[int(sys.argv[2])]

orf_file=sys.argv[3]

for i in range(len(genomes)):

	genome=genomes[i]
	size=sizes[i]

	dprot_coord=defaultdict(list)
	dprot_clu=defaultdict(str)
	dprot_color=defaultdict(str)

	with open(orf_file,"r") as f1:
		for li in f1:
			li=li.rstrip()
			lp=li.split()

			if li.startswith(">") and re.sub("_[0-9]+$","",lp[0].lstrip(">"))==genome and "REVERSE" not in li:
				dprot_coord[lp[0].lstrip(">")+"+"]=[int(lp[1].lstrip("[")),int(lp[3].rstrip("]"))]
	
			elif li.startswith(">") and re.sub("_[0-9]+$","",lp[0].lstrip(">"))==genome and "REVERSE" in li:
				dprot_coord[lp[0].lstrip(">")+"-"]=[int(lp[3].rstrip("]")),int(lp[1].lstrip("["))]


	##########################


	#To have each prot in its frame
	dprot_frame=defaultdict(str)

	for prot in dprot_coord:
		if re.search("\+$",prot):
			if dprot_coord[prot][0] in list(range(1,size,3)):	
					dprot_frame[prot]=1

			elif dprot_coord[prot][0] in list(range(2,size,3)):
					dprot_frame[prot]=2

			elif dprot_coord[prot][0] in list(range(3,size,3)):
					dprot_frame[prot]=3	

		else:
			if size-(dprot_coord[prot][1]-1) in list(range(1,size,3)):
				dprot_frame[prot]=-1

			elif size-(dprot_coord[prot][1]-1) in list(range(2,size,3)):
				dprot_frame[prot]=-2

			elif size-(dprot_coord[prot][1]-1) in list(range(3,size,3)):
				dprot_frame[prot]=-3	

	#############################

	#pNpS coloring

	#PNPS data
	 
	dprot_pnps=defaultdict(float)

	dprot_fb=defaultdict(float)

	les_bleus=list()

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

			if float(lp[-3]) <= 1.2:
				les_bleus.append(lp[0])

	for prot in dprot_coord:
		if prot not in dprot_pnps:
			dprot_pnps[prot]="NA"

	#Finding the mirror strand:
	if size%3 == 0:
		mirror=[(2,-1),(1,-2),(3,-3)]

	elif size%3 == 1:
		mirror=[(3,-1),(2,-2),(1,-3)]

	elif size%3 == 2:
		mirror=[(3,-2),(2,-3),(1,-1)]

	rouges=set()

	for orf1 in les_bleus:
		a=re.search("\+|-$",orf1)
		if a:
			strand1=a.group(0)

		for orf2 in les_bleus:
			b=re.search("\+|-$",orf2)
			if b:
				strand2=b.group(0)
		
			if strand1!=strand2:
				if (dprot_frame[orf1],dprot_frame[orf2]) in mirror or (dprot_frame[orf2],dprot_frame[orf1]) in mirror:

					start1=dprot_coord[orf1][0]
					end1=dprot_coord[orf1][1]
					size1=end1-start1

					start2=dprot_coord[orf2][0]
					end2=dprot_coord[orf2][1]
					size2=end2-start2

					if start2 <= start1 <= end2 and start2 <= end1 <= end2:
						if (end1-start1) / size1 >= 0.75:
							
							if size1 > size2:
								rouges.add(orf2.rstrip("+-"))

							else:
								rouges.add(orf1.rstrip("+-"))						

					elif  start2 <= start1 <= end2:
						if (end2-start1) / size1 >= 0.75 or (end2-start1) / size2 >= 0.75 :
							
							if size1 > size2:
								rouges.add(orf2.rstrip("+-"))

							else:
								rouges.add(orf1.rstrip("+-"))

					elif start2 <= end1 <= end2:
						if (end1-start2) / size1 >= 0.75 or (end1-start2) / size2 >= 0.75:
							
							if size1 > size2:
								rouges.add(orf2.rstrip("+-"))

							else:
								rouges.add(orf1.rstrip("+-"))

					if start1 <= start2 <= end1 and start1 <= end2 <= end1:
						if (end2-start2) / size2 >= 0.75:
							
							if size1 > size2:
								rouges.add(orf2.rstrip("+-"))

							else:
								rouges.add(orf1.rstrip("+-"))						

					elif  start1 <= start2 <= end1:
						if (end1-start2) / size2 >= 0.75 or (end1-start2) / size1 >= 0.75:
							
							if size1 > size2:
								rouges.add(orf2.rstrip("+-"))

							else:
								rouges.add(orf1.rstrip("+-"))

					elif start1 <= end2 <= end1:
						if (end2-start1) / size2 >= 0.75 or (end2-start1) / size1 >= 0.75:
							
							if size1 > size2:
								rouges.add(orf2.rstrip("+-"))

							else:
								rouges.add(orf1.rstrip("+-"))
