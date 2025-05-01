
import warnings
warnings.filterwarnings('ignore')

import sys
import argparse
import plotnine as p9
from plotnine import ggplot, aes, geom_line
from plotnine import *
import pandas as pd


def get_frac(args, infilename):

	infile = open(infilename, "r")
	lines = infile.readlines()

	nots = [0,1,2]
	dictionary = {}
	outlist = []
	total = 0

	nraround = args.nr_around

	for line in lines:
		line = line.strip()
		plotindex = int(line.split("\t")[1])
		nr = int(line.split("\t")[-1])

		for i in range((nraround)*-3,(nraround*3)+3):
			if plotindex == i:
				if plotindex not in nots:
					if plotindex < 0:
						total = total + nr
						dictionary[plotindex] = nr
					else:
						plotindex = plotindex-2
						total = total + nr
						dictionary[plotindex] = nr

	indexlist = list(dictionary.keys())

	for k,v in dictionary.items():
		outlist.append(v/total)

	return outlist, indexlist


def allparsers(parser):

	parser.add_argument("-is","--input_substitution_file", type = str, help = "Inputfile", required = True)
	parser.add_argument("-iwt","--input_wt_file", type = str, help = "Reference file", required = True)
	parser.add_argument("-nra","--nr_around", type = int, help= "nr around the codon of interest (default 10)", required = False, default = 10)
	parser.add_argument("-o","--outfileprefix", type =str, help= "name of outfile prefix", required = True)
	parser.add_argument("-wo","--writeout", help="Flag to use if you want the dataframe to be written into a textfile" , action='store_true')

def get_diff(args):

	list1 = get_frac(args, args.input_substitution_file)[0]
	list2 = get_frac(args, args.input_wt_file)[0]
	indexlist = get_frac(args, args.input_substitution_file)[1]
	difflist = []

	for i in range(len(list1)):
		difflist.append(str(list1[i] - list2[i]))

	df = pd.DataFrame({"index":indexlist, "f1":list1, "f2":list2,"diff":difflist}, dtype = float)

	return df


def plot_diff(diff_dataframe, args):

	outnamepre = args.outfileprefix
	plotoutname = outnamepre + ".pdf"

	nraround = args.nr_around

	xlabname = "pos"
	ylabname = "frequency difference" 

	plotp = (p9.ggplot(data=diff_dataframe, mapping=p9.aes(x='index', y='diff'))) + geom_line(color = "red", size = 2) + p9.scale_x_continuous(breaks=range((nraround)*-3,(nraround*3)+1),minor_breaks= [1]) + geom_vline(xintercept=0, linetype = "dashed") + geom_hline(yintercept=0, linetype ="dashed") + ylab(ylabname) + xlab(xlabname) + ylim(-0.02,0.02)
	plotp.save(plotoutname,  width = 20, height =8, dpi=1000)


def write_out(diff_dataframe, args):

	outnamepre = args.outfileprefix
	textoutname = outnamepre + ".txt"

	diff_dataframe.to_csv(textoutname, sep='\t', index=False,)

def main():

	parser = argparse.ArgumentParser()
	allparsers(parser)
	args = parser.parse_args()

	diff_dataframe = get_diff(args)

	plot_diff(diff_dataframe, args)

	if args.writeout:
		write_out(diff_dataframe, args)


main()




