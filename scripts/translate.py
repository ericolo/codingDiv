#!/usr/bin/env python3

import os 
import re
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

#transalate

seqs=SeqIO.parse(sys.argv[1],"fasta")


#with open("my_output.fasta","w") as f1:

for record in seqs:

	frame1=Seq(str(record.seq)[0::])
	frame2=Seq(str(record.seq)[1::])
	frame3=Seq(str(record.seq)[2::])

	print(">"+str(record.id)+"_"+"frame_1"+"_"+"5'3'")
	print(frame1.translate(table=1))


	print(">"+str(record.id)+"_"+"frame_2"+"_"+"5'3'")
	print(frame2.translate(table=1))


	print(">"+str(record.id)+"_"+"frame_3"+"_"+"5'3'")
	print(frame3.translate(table=1))

	
	revcomp=record.seq.reverse_complement()

	frame_1=Seq(str(revcomp)[0::])
	frame_2=Seq(str(revcomp)[1::])
	frame_3=Seq(str(revcomp)[2::])


	print(">"+str(record.id)+"_"+"frame_-1"+"_"+"3'5'")
	print(frame_1.translate(table=1))


	print(">"+str(record.id)+"_"+"frame_-2"+"_"+"3'5'")
	print(frame_2.translate(table=1))


	print(">"+str(record.id)+"_"+"frame_-3"+"_"+"3'5'")
	print(frame_3.translate(table=1))

		#os.system('curl -s -d "dna_sequence='+str(record.seq)+'&output_format=fasta" https://web.expasy.org/cgi-bin/translate/dna2aa.cgi >> my_output.fasta' )


