#!/usr/bin/env python3

from collections import defaultdict
import re 
from Bio import SeqIO
from Bio import Seq
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
import glob
import sys
from Bio.Seq import Seq

#Script checking if mapped snps changed amino acids
#SPECIFIC FOR MINNOW THAT HAVE A REVERSE CONSERVED PROT


###### This script takes into account several subs per position

#clu glob, use \* on the terminal
prots=sys.argv[1]

snp=sys.argv[2]

genomes=sys.argv[3]

genetic_code=sys.argv[4]

######### genomes file

genome_parser1 = SeqIO.parse(genomes,"fasta")

dgenome_size=defaultdict(int)

for record in genome_parser1:
	dgenome_size[str(record.id)]=len(str(record.seq))


######### SNP file

dgenome_snp=defaultdict(list)

with open(snp,"r") as f2:
	for li in f2:
		li=li.rstrip()
		lp=li.split("\t")

		current_genome=lp[0]

		depths=list(map(int,lp[-1].split(",")))
		total_depth=sum(depths)

		if lp[4]!=".": #key=position, values=[ref_nuc,alt_nuc,total_depth,[depth1,depth2,depth3,depth4],score,total_snp_depth]

			if len(lp[4])==1:		

				#forward
				dgenome_snp[(int(lp[1]),lp[4])]=[lp[3],lp[4],total_depth,depths[1],lp[5],sum(depths[1::]),int(lp[1]),"forward"]

				#reverse
				dgenome_snp[( dgenome_size[current_genome]-(int(lp[1])-1) , str(Seq(lp[4]).complement()) )]=[str(Seq(lp[3]).complement()),str(Seq(lp[4]).complement()),total_depth,depths[1],lp[5],sum(depths[1::]),int(lp[1]),"reverse"]

			#If several SNPs at the same position
			else:

				count=0
				for alt_nuc in re.sub(",","",lp[4]):

					count+=1

					dgenome_snp[(int(lp[1]),alt_nuc)]=[lp[3],alt_nuc,total_depth,depths[count],lp[5],sum(depths[1::]),int(lp[1]),"forward"]

					dgenome_snp[( dgenome_size[current_genome]-(int(lp[1])-1) , str(Seq(alt_nuc).complement()) )]=[str(Seq(lp[3]).complement()),str(Seq(alt_nuc).complement()),total_depth,depths[count],lp[5],sum(depths[1::]),int(lp[1]),"reverse"]



######### genomes file

current_genome_seq=str()

genome_parser = SeqIO.parse(genomes,"fasta")

for record in genome_parser:
	if str(record.id)==current_genome:
		
		genome_seq=">"+str(record.seq) #added a > so I don't have to worry about indexes starting at 0 when extracting a substring
		genome_size=len(str(record.seq))

		rev_genome_seq=">"+str(record.seq.reverse_complement())


######### proteins file

#Be careful when orfs cross the breakpoint

dprot_coord=defaultdict(list)

dprot_seq=defaultdict(str)

dprot_clu=defaultdict(str)

lcoord_prot_list = [[] for x in range(genome_size+1)]

for file_ in glob.glob(prots):
	prot_parser = SeqIO.parse(file_,"fasta")
	
	for record in prot_parser:

		if re.sub("_[0-9]+$","",str(record.id))==current_genome:

			if "REVERSE" not in str(record.description):

				#adding a tag to the prot name to know the strand
				dprot_seq[str(record.id)+"+"]=">"+str(record.seq) #added a > so I don't have to worry about indexes starting at 0 when extracting a substring

				start=str(record.description).split()[1].strip("[]")
				end=str(record.description).split()[3].strip("[]")

				dprot_coord[str(record.id)+"+"]=[int(start),int(end)]

				for pos in range(int(start),int(end)):
					lcoord_prot_list[pos].append(str(record.id)+"+")

			else:

				dprot_seq[str(record.id)+"-"]=">"+str(record.seq)

				start=str(record.description).split()[1].strip("[]")
				end=str(record.description).split()[3].strip("[]")
				
				dprot_coord[str(record.id)+"-"]=[genome_size-(int(start)-1),genome_size-(int(end)-1)]

				for pos in range(genome_size-(int(start)-1),genome_size-(int(end)-1)):
					lcoord_prot_list[pos].append(str(record.id)+"-")



######### creating blosum62 dict and properties dict

blosum_dict=defaultdict(float)

dprop=defaultdict(str)

for key in matlist.blosum62:
	blosum_dict[(key[0],key[1])]=matlist.blosum62[key]
	blosum_dict[(key[1],key[0])]=matlist.blosum62[key]

	if key[0] in "RHK":
		dprop[key[0]]="Charged+"
	elif key[0] in "DE":
		dprop[key[0]]="Charged-"
	elif key[0] in "STNQ":
		dprop[key[0]]="Polar"
	elif key[0] in "AVILMFYW":
		dprop[key[0]]="Hydrophobic"
	elif key[0] in "CGP":
		dprop[key[0]]="Special"

	if key[1] in "RHK":
		dprop[key[1]]="Charged+"
	elif key[1] in "DE":
		dprop[key[1]]="Charged-"
	elif key[1] in "STNQ":
		dprop[key[1]]="Polar"
	elif key[1] in "AVILMFYW":
		dprop[key[1]]="Hydrophobic"
	elif key[1] in "CGP":
		dprop[key[1]]="Special"

dprop["*"]="STOP"

#Frames list:

dpos_frame=defaultdict(str)

for i in range(1,genome_size,3):
	dpos_frame[i]=1

for i in range(2,genome_size,3):
	dpos_frame[i]=2

for i in range(3,genome_size,3):
	dpos_frame[i]=3

######### Check

print("prot","clu","position_on_genome","corrected_position","frame","codon_position","ref_codon","ref_aa","ref_type","####","alt_codon","alt_aa","alt_type","blosum62","total_depth","snp_depth","total_snp_depth","qual_score",sep="\t")

for snp in dgenome_snp:
	#print(snp)

	for prot in lcoord_prot_list[snp[0]]:
		count=0


		#Impacted nuc
		if (snp[0] % 3) == 0:
			nuc_index=3

		elif (snp[0] % 3) == 2:
			nuc_index=2

		elif (snp[0] % 3) == 1: 
			nuc_index=1

		#to get strand of snp:
		a=re.search("[-|+]$",prot)
		if a:
			current_strand=a.group(0)

		if dgenome_snp[snp][-1]=="forward" and current_strand=="+":
			current_genome_seq=genome_seq


			if (dpos_frame[dprot_coord[prot][0]]==1 and nuc_index==1) or (dpos_frame[dprot_coord[prot][0]]==2 and nuc_index==2) or (dpos_frame[dprot_coord[prot][0]]==3 and nuc_index==3):
				ref_codon=current_genome_seq[snp[0]]+current_genome_seq[snp[0]+1]+current_genome_seq[snp[0]+2]

				alt_codon=dgenome_snp[snp][1]+current_genome_seq[snp[0]+1]+current_genome_seq[snp[0]+2]

				ref_aa=Seq.translate(ref_codon,table=genetic_code)

				alt_aa=Seq.translate(alt_codon,table=genetic_code)

				codon_position=snp[0]


			elif (dpos_frame[dprot_coord[prot][0]]==2 and nuc_index==1) or (dpos_frame[dprot_coord[prot][0]]==1 and nuc_index==3) or (dpos_frame[dprot_coord[prot][0]]==3 and nuc_index==2):
				ref_codon=current_genome_seq[snp[0]-2]+current_genome_seq[snp[0]-1]+current_genome_seq[snp[0]]

				alt_codon=current_genome_seq[snp[0]-2]+current_genome_seq[snp[0]-1]+dgenome_snp[snp][1]

				ref_aa=Seq.translate(ref_codon,table=genetic_code)

				alt_aa=Seq.translate(alt_codon,table=genetic_code)

				codon_position=snp[0]-2


			elif (dpos_frame[dprot_coord[prot][0]]==3 and nuc_index==1) or (dpos_frame[dprot_coord[prot][0]]==2 and nuc_index==3) or (dpos_frame[dprot_coord[prot][0]]==1 and nuc_index==2):
				ref_codon=current_genome_seq[snp[0]-1] + current_genome_seq[snp[0]] + current_genome_seq[snp[0]+1]

				alt_codon=current_genome_seq[snp[0]-1] + dgenome_snp[snp][1] + current_genome_seq[snp[0]+1]

				ref_aa=Seq.translate(ref_codon,table=genetic_code)

				alt_aa=Seq.translate(alt_codon,table=genetic_code)

				codon_position=snp[0]-1


			if alt_aa!="*" and ref_aa!="*":
				print(prot,dprot_clu[prot],dgenome_snp[snp][-2],snp[0],str(dpos_frame[dprot_coord[prot][0]]),codon_position,ref_codon,ref_aa,dprop[ref_aa],"####",alt_codon,alt_aa,dprop[alt_aa],blosum_dict[(ref_aa,alt_aa)],dgenome_snp[snp][2],dgenome_snp[snp][3],dgenome_snp[snp][5],dgenome_snp[snp][4],sep="\t")

			else:
				print(prot,dprot_clu[prot],dgenome_snp[snp][-2],snp[0],str(dpos_frame[dprot_coord[prot][0]]),codon_position,ref_codon,ref_aa,dprop[ref_aa],"####",alt_codon,alt_aa,dprop[alt_aa],"NA",dgenome_snp[snp][2],dgenome_snp[snp][3],dgenome_snp[snp][5],dgenome_snp[snp][4],sep="\t")


		elif dgenome_snp[snp][-1]=="reverse" and current_strand=="-":
			current_genome_seq=rev_genome_seq

			if (dpos_frame[dprot_coord[prot][0]]==1 and nuc_index==1) or (dpos_frame[dprot_coord[prot][0]]==2 and nuc_index==2) or (dpos_frame[dprot_coord[prot][0]]==3 and nuc_index==3):
				ref_codon=current_genome_seq[snp[0]]+current_genome_seq[snp[0]+1]+current_genome_seq[snp[0]+2]

				alt_codon=dgenome_snp[snp][1]+current_genome_seq[snp[0]+1]+current_genome_seq[snp[0]+2]

				ref_aa=Seq.translate(ref_codon,table=genetic_code)

				alt_aa=Seq.translate(alt_codon,table=genetic_code)

				codon_position=snp[0]


			elif (dpos_frame[dprot_coord[prot][0]]==2 and nuc_index==1) or (dpos_frame[dprot_coord[prot][0]]==1 and nuc_index==3) or (dpos_frame[dprot_coord[prot][0]]==3 and nuc_index==2):
				ref_codon=current_genome_seq[snp[0]-2]+current_genome_seq[snp[0]-1]+current_genome_seq[snp[0]]

				alt_codon=current_genome_seq[snp[0]-2]+current_genome_seq[snp[0]-1]+dgenome_snp[snp][1]

				ref_aa=Seq.translate(ref_codon,table=genetic_code)

				alt_aa=Seq.translate(alt_codon,table=genetic_code)

				codon_position=snp[0]-2


			elif (dpos_frame[dprot_coord[prot][0]]==3 and nuc_index==1) or (dpos_frame[dprot_coord[prot][0]]==2 and nuc_index==3) or (dpos_frame[dprot_coord[prot][0]]==1 and nuc_index==2):
				ref_codon=current_genome_seq[snp[0]-1] + current_genome_seq[snp[0]] + current_genome_seq[snp[0]+1]

				alt_codon=current_genome_seq[snp[0]-1] + dgenome_snp[snp][1] + current_genome_seq[snp[0]+1]

				ref_aa=Seq.translate(ref_codon,table=genetic_code)

				alt_aa=Seq.translate(alt_codon,table=genetic_code)

				codon_position=snp[0]-1


			if alt_aa!="*" and ref_aa!="*":
				print(prot,dprot_clu[prot],dgenome_snp[snp][-2],snp[0],current_strand+str(dpos_frame[dprot_coord[prot][0]]),codon_position,ref_codon,ref_aa,dprop[ref_aa],"####",alt_codon,alt_aa,dprop[alt_aa],blosum_dict[(ref_aa,alt_aa)],dgenome_snp[snp][2],dgenome_snp[snp][3],dgenome_snp[snp][5],dgenome_snp[snp][4],sep="\t")

			else:
				print(prot,dprot_clu[prot],dgenome_snp[snp][-2],snp[0],current_strand+str(dpos_frame[dprot_coord[prot][0]]),codon_position,ref_codon,ref_aa,dprop[ref_aa],"####",alt_codon,alt_aa,dprop[alt_aa],"NA",dgenome_snp[snp][2],dgenome_snp[snp][3],dgenome_snp[snp][5],dgenome_snp[snp][4],sep="\t")

