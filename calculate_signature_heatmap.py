
import sys
import re
import argparse
import os, fnmatch
import pandas as pd

def allparsers(parser):

	parser.add_argument("-i1","--inputfolder1", type = str, help = "Name of inputfolder1", required = True)
	parser.add_argument("-i2","--inputfolder2", type = str, help = "Name of inputfolder2", required = True)
	parser.add_argument("-o","--outfolder", type =str, help= "name of outfolder", required = True)


def runx(inputfoldername, all_ACGT, infilename):

	infilename = inputfoldername + infilename
	infile = open(infilename, "r")

	lines = infile.readlines()

	for line in lines:
		line = line.strip()
		pos = line.split("\t")[1]
		value = line.split("\t")[-1]
		nuc = line.split("\t")[2]

		if int(pos) in range(-3,6):
			all_ACGT[str(pos)][nuc] = value

	return all_ACGT


def empty_ACGT():

	nucleotides = ["A","C","G","T"]
	all_ACGT = {}

	for pos in range(-3,6):
		all_ACGT[str(pos)] = {}
		for nuc in nucleotides:
			all_ACGT[str(pos)][nuc] = 0

	return all_ACGT

def calc_pd_df(dict):

	df = pd.DataFrame.from_dict(dict)
	df.astype('int').dtypes

	cols = df.columns
	df[cols] = pd.to_numeric(df[cols].stack(), errors='coerce').unstack()

	for pos in range(-3,6):
		df_fracstr = str(pos) +  "_frac"
		df[df_fracstr] = df[str(pos)]/df[str(pos)].sum()

	return df


def calc_difference(df1, df2):

	diffdf = pd.DataFrame()

	for pos in range(-3,6):
		df_fracstr = str(pos) +  "_frac"
		diff_fracstr = "diff_" + df_fracstr
		diffdf[diff_fracstr] = df1[df_fracstr]-df2[df_fracstr]

	diffdf = diffdf.drop("diff_0_frac", axis=1) # drop columns 0,1,2 as tehse are of the codon itself
	diffdf = diffdf.drop("diff_1_frac", axis=1)
	diffdf = diffdf.drop("diff_2_frac", axis=1)

	diffdf.columns = ["-3", "-2", "-1", "1", "2","3"]

	diffdf["nucleotide"] = diffdf.index

	return diffdf


def sep_by_codon(inputfolder, codonfiles):

	codonseps = {}

	for codonfile in codonfiles:
		codon = codonfile.split("codon.txt")[0].split("_")[1][:3]
		nuc =  codonfile.split("codon.txt")[0].split("_")[1][-1]

		if codon not in codonseps:
			codonseps[codon] = [codonfile]
		else:
			codonseps[codon].append(codonfile)

	maindict = {}

	for codon in codonseps.keys():

		cod_ACGT = empty_ACGT()
		maindict[codon] = cod_ACGT

		for codonfile in codonseps[codon]:
                        #print(codon,codonfile)
			maindict[codon] = runx(inputfolder, cod_ACGT, codonfile)

	return maindict

#	return codonseps


def main():

#	print("") #

	parser = argparse.ArgumentParser()
	allparsers(parser)
	args = parser.parse_args()

	if not os.path.isdir(args.inputfolder1):
		print("Folder given with -i1 ", args.inputfolder1, "does not exist")
		sys.exit()

	if not os.path.isdir(args.inputfolder2):
		print("Folder given with -i2 ", args.inputfolder2, "does not exist")
		sys.exit()

	nuc_ACGT1 = empty_ACGT()
	nuc_ACGT2 = empty_ACGT()

	nucfiles1 = fnmatch.filter(os.listdir(args.inputfolder1), '*nuc.txt')
	nucfiles2 = fnmatch.filter(os.listdir(args.inputfolder2), '*nuc.txt')

	for nucfile in nucfiles1:
		nuc_ACGT1 = runx(args.inputfolder1, nuc_ACGT1, nucfile)

	for nucfile in nucfiles2:
		nuc_ACGT2 = runx(args.inputfolder2, nuc_ACGT2, nucfile)

	nuc_df1 = calc_pd_df(nuc_ACGT1)
	nuc_df2 = calc_pd_df(nuc_ACGT2)
	nuc_diffdf = calc_difference(nuc_df1, nuc_df2)

	nuc_outfilename = args.outfolder + "/Nucleotide_frequency_difference.txt"
	nuc_diffdf.to_csv(nuc_outfilename, sep ='\t', index = False)
	print("nucleotidefile", nuc_outfilename, " finished")

##############################################

	codonfiles1 = fnmatch.filter(os.listdir(args.inputfolder1), '*codon.txt')
	codon_dict1 = sep_by_codon(args.inputfolder1, codonfiles1)

	codonfiles2 = fnmatch.filter(os.listdir(args.inputfolder2), '*codon.txt')
	codon_dict2 = sep_by_codon(args.inputfolder2, codonfiles2)

	for codon in codon_dict1.keys():
		cod_df1 = calc_pd_df(codon_dict1[codon])
		cod_df2 = calc_pd_df(codon_dict2[codon])
		cod_diffdf = calc_difference(cod_df1, cod_df2)
		cod_outfilename = args.outfolder + "/" + str(codon) + "_Codon_frequency_difference.txt"
		cod_diffdf.to_csv(cod_outfilename, sep ='\t', index = False)

		print("codonfile", cod_outfilename, " finished")

main()
