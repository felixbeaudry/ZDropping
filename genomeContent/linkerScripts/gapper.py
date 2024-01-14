#python3
#march 23, 2021
#Felix

import fileinput, argparse
import sys
import getopt
from collections import Counter

def arguments():
        parser  = argparse.ArgumentParser(description="window results from fasta assembly summary")
        parser.add_argument("-i","--input",help="input file",required=True)
        parser.add_argument("-w","--window",help="window size in kb",required=True)
        args = parser.parse_args()
        return(args)

args = arguments()


reader = open(args.input)
win = int(args.window) * 1000

scaf = None
scafName = ""

for record in reader:
	#on first loop, just record scaffold names
	if scaf == None:
		scaf = ""
		scafName = record.rstrip()
		print("scaffold bp N G C A T")
	elif record.startswith(">"):
		#subset previous scaffold
		#chunks,chunk_size=len(scaf),len(scaf)//win
		subset_scaf = [ scaf[i:i+win] for i in range(0,len(scaf),win) ]
		#print(chunks,chunk_size,len(subset_scaf))
		for i in range(0,len(subset_scaf)):
			letterCount = Counter(subset_scaf[i])
			print(scafName+" ",end='',flush=True)
			print((i*win),letterCount['N'],letterCount['G'],letterCount['C'],letterCount['T'],letterCount['A'])
		#reset for next scaffold
		scaf = ""
		scafName = record.rstrip()
	else:
		scaf=scaf+record.rstrip()


