
import sys
import re
import argparse
import os

def allparsers(parser):

	parser.add_argument("-i","--inputfile", type = str, help = "Name of inputfile, with IDs in first column", required = True)
	parser.add_argument("-r","--reference_file", type = str, help = "Reference mRNA fasta file", required = True)
	parser.add_argument("-s","--substitution", type =str, help= "substitution to perform e.g HQ", required = True)
	parser.add_argument("-o","--outfolder", type =str, help= "name of outfolder (will also be used as pre for file names)", required = True)
	parser.add_argument("-t","--type", type = str, help = "is it Substitutant or WT peptides (S or WT) ?",required = True)

def return_codondict():

	codon2AA = {
		  'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
		  'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
		  'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
		  'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
		  'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
		  'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
		  'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
		  'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
		  'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
		  'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
		  'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
		  'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
		  'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
		  'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
		  'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
		  'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

	return codon2AA

def codon_to_AA(codon2AA, codon_seq):

	AAseq = []

	for codon in codon_seq:

		if codon not in codon2AA.keys(): break

		if codon2AA[codon] == "*" : break

		AAseq.append(codon2AA[codon])

	AAseq = "".join(AAseq)

	return AAseq

##################3

#stopcodons = ["TAG","TAA","TGA"]

def dnatocodon(codon2AA, seq):

	codons = []

	for i in range(0,len(seq),3):

		if codon2AA[seq[i:i+3]] == "*": break

		codons.append(seq[i:i+3])

	return codons


##############

def load_REFfile(args, codon2AA):


	REFfile = open(args.reference_file, "r")

	lines = REFfile.readlines()

	REF = {}

	for line in lines:
		line = line.strip()

		if line.startswith(">"):
			header = line
			ENST = header.split("|")[0][1:]
		else:
			seq = line
			codonseq = dnatocodon(codon2AA, seq)
			AAseq = codon_to_AA(codon2AA, codonseq)
			if ENST not in REF:
				REF[ENST] = {}
				REF[ENST]["rna"] = codonseq
				REF[ENST]["prot"] = AAseq

	return REF

######################

def invert_dict(dict):

	newdict = {}

	for k,v in dict.items():
		if v not in newdict.keys():
			newdict[v] = [k]
		else:
			newdict[v].append(k)

	return newdict

def getnucmatrix(args, codon2AA):

	allmatrix = {}
	nraround = 10

	AA2codon = invert_dict(codon2AA)
	usedcodons = AA2codon[args.substitution[0]]
	usedcodons.append("all")

	for codon in usedcodons:
		allmatrix[codon] = {}
		for i in range(((nraround*2)+1)*3):
			allmatrix[codon] [i] = {"A":0, "C":0, "T":0, "G":0}

	AAmatrix = {}
	for i in range((nraround*2)+1):
		AAmatrix[i] = {"A":0, "C":0, "D":0, "E":0, "F":0, "G":0, "H":0, "I":0, "K":0, "L":0, "M":0, "N":0, "P":0, "Q":0, "R":0, "S":0, "T":0, "V":0, "W":0, "Y":0}
	allmatrix["AA"] = AAmatrix

	return allmatrix

###########################################################

def mappeptides(args, allmatrix, REF, codon2AA):

	bef_aft = []
	nraround = 10

	infile = open(args.inputfile, "r")
	lines = infile.readlines()

	AA2codon = invert_dict(codon2AA)
	usedcodons = AA2codon[args.substitution[0]]


	for line in lines:
		line = line.strip()

		if not line.startswith("ID"):
			ID = line.split("\t")[0]
			seq = ID.split("_")[-1]
			ENST = ID.split("_")[0]
			Gene = ID.split("_")[1]
			codon = ID.split("_")[2]

			xstring = "all" + args.substitution

			if args.type == "WT":
				group = ID.split("_")[3]
				start = int(ID.split("_")[4])
				end = int(ID.split("_")[5])

				if ENST in REF:
					AAseq = REF[ENST]["prot"]
					codonseq = REF[ENST]["rna"]

					peptide_AA = AAseq[start:end]
					peptide_codonseq = codonseq[start:end]

					for codon in usedcodons:
						multiples = []
						multiples.extend([i for i, x in enumerate(peptide_codonseq) if x == codon])

						for subcodon in multiples:
							relativepos = int(start) + int(subcodon)
							aroundrna = "".join(REF[ENST]["rna"][relativepos-10:relativepos+11])
							aroundseq = AAseq[relativepos-10:relativepos+11]

							targetcodon = REF[ENST]["rna"][relativepos]
							followcodon = REF[ENST]["rna"][relativepos+1]
							beforecodon = REF[ENST]["rna"][relativepos-1]

							if len(aroundseq) == 21:
								bef_aft.append([beforecodon, targetcodon, followcodon])

								for i in range(((nraround*2)+1)*3):
									nuc = aroundrna[i]
									allmatrix[codon][i][nuc] = allmatrix[codon][i][nuc] + 1

								for i in range(((nraround*2)+1)*3):
									nuc = aroundrna[i]
									allmatrix["all"][i][nuc] = allmatrix["all"][i][nuc] + 1

								for i in range((nraround*2)+1):
									AA = aroundseq[i]
									allmatrix["AA"][i][AA] = allmatrix["AA"][i][AA] + 1
				
			else:

				#print(codon)

				if not codon == xstring:
					pos = int(ID.split("_")[3])
					start = int(ID.split("_")[5])
					end = int(ID.split("_")[6])

					if ENST in REF:
						AAseq = REF[ENST]["prot"]

						aroundseq = AAseq[pos-10:pos+11]
						if len(aroundseq) == 21:

							targetcodon = REF[ENST]["rna"][pos]
							followcodon = REF[ENST]["rna"][pos+1]
							beforecodon = REF[ENST]["rna"][pos-1]
							bef_aft.append([beforecodon, targetcodon, followcodon])

							aroundrna = "".join(REF[ENST]["rna"][pos-10:pos+11])

							for key in allmatrix.keys():

								if key != "AA":

									if key == "all":

										for i in range(((nraround*2)+1)*3):
											nuc = aroundrna[i]
											allmatrix[key][i][nuc] = allmatrix[key][i][nuc] + 1
									else:
										for i in range(((nraround*2)+1)*3):
											nuc = aroundrna[i]
											allmatrix[codon][i][nuc] = allmatrix[codon][i][nuc] + 1	
#
								else:

									for i in range((nraround*2)+1):
										AA = aroundseq[i]
										allmatrix["AA"][i][AA] = allmatrix["AA"][i][AA] + 1

	return allmatrix

##################################################################


def writeouts(args, allmatrix, codon2AA):

	nucs = ["A","C","G","T"]
	nraround = 10

	for nuc in nucs:
		outname = args.outfolder + "/" + args.substitution + "_" + nuc + "nuc.txt"
		nucfile = open(outname, "w")

		for pos, rest in allmatrix["all"].items():
			for restnuc in rest.keys():
				if restnuc == nuc:
					outstrng = str(pos) + "\t" + str(pos-(nraround*3)) + "\t" + str(nuc) + "\t" + str(rest[nuc]) + "\n"
					nucfile.write(outstrng)
		nucfile.close()
	##

	AA2codon = invert_dict(codon2AA)
	usedcodons = AA2codon[args.substitution[0]]

	for codon in usedcodons:

		for nuc in nucs:
			outname = args.outfolder + "/" + args.substitution + "_" + codon + nuc + "codon.txt"
			nucfile = open(outname, "w")

			for pos, rest in allmatrix[codon].items():
				for restnuc in rest.keys():
					if restnuc == nuc:
						outstrng = str(pos) + "\t" + str(pos-(nraround*3)) + "\t" + str(nuc) + "\t" + str(rest[nuc]) + "\n"
						nucfile.write(outstrng)
			nucfile.close()

	##

	for AA in AA2codon.keys():

		alreadypos = []

		if AA != "*":
			outname = args.outfolder + "/" + args.substitution + "_" + AA + "_AA.txt"
			AAfile = open(outname, "w")

			for pos, rest in allmatrix["AA"].items():
				for restAA in rest.keys():

					if pos not in alreadypos:
						outstrng = str(pos) + "\t" + str(pos-nraround) + "\t" + str(AA) + "\t" + str(rest[AA]) + "\n"
						AAfile.write(outstrng)
						alreadypos.append(pos)
			AAfile.close()

#############################

def main():

	parser = argparse.ArgumentParser()
	allparsers(parser)
	args = parser.parse_args()

	if not os.path.isdir(args.outfolder):
		print("Folder given with -o ", args.outfolder, "does not exist") 
		sys.exit()

	codon2AA = return_codondict()

	REF = load_REFfile(args, codon2AA)

	print("loaded reference file")

	allmatrix = getnucmatrix(args, codon2AA)
	allmatrix = mappeptides(args, allmatrix, REF, codon2AA)

	print("created and filled matrices")

	writeouts(args, allmatrix, codon2AA)

	print("write all matrices into outfiles")

main()

