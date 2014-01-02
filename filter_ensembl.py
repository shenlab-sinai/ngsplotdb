import os
import sys
import string
import re

def main():
	inFileName = sys.argv[1]
	dbTxt = sys.argv[2]
	outFileName = sys.argv[3]
	species = sys.argv[4]

	annoDict = {}
	inFile = file(inFileName)
	allowedType = ["lincRNA", "miRNA", "pseudogene", "protein_coding"]
	p = re.compile(r'.*\S*; gene_name "(\S*;\S*)"; transcript_name.*')
	while True:
		line = inFile.readline()
		if len(line) == 0:
			break;
		line = string.strip(line)
		annoPart = string.split(line, sep="\t")[-1][:-1].replace("\"", "")
		annoPartL = string.split(annoPart, sep=";")
		source = line.split("\t")[1]  # LS.
		try:
			# geneDict = dict(item.split(" ")[1:3] for item in annoPart.split(";"))
			geneDict = dict(item.strip().split(" ")[:2] for item in annoPart.split(";"))  # LS.
		except:
			raw_name = p.match(line).group(1)
			target_name = raw_name.replace(";", "-")
			line = string.replace(line, raw_name, target_name)
			annoPart = string.split(line, sep="\t")[-1][:-1].replace("\"", "")
			annoPartL = string.split(annoPart, sep=";")
			geneDict = dict(item.split(" ")[1:3] for item in annoPart.split(";"))
		
		if "gene_biotype" not in geneDict:
			biotype = source
		else:
			biotype = geneDict["gene_biotype"]
		if biotype in allowedType:
			annoDict[geneDict["gene_id"]] = biotype
		else:
			annoDict[geneDict["gene_id"]] = "misc"

	inFile.close()

	inFile = file(dbTxt)
	outFile = file(outFileName, "w")
	# pattern = re.compile(r'[0-9]|X|Y|x|y|MT|Mt|M')
	while True:
		line = inFile.readline()
		if len(line) == 0:
			break
		line = string.strip(line)
		lineL = string.split(line, sep="\t")
		if lineL[0].find("_") > 0 or lineL[0].find(".") > 0:
			continue
		#if (pattern.match(lineL[0]) is not None):
		lineL[0] = "chr" + lineL[0]
		if (lineL[0] == "chrMt" or lineL[0] == "chrMT"):
			lineL[0] = "chrM"
		if lineL[3] in annoDict:
			lineL.append(annoDict[lineL[3]])
		else:
			lineL.append("unknown")
		newLine = "\t".join(lineL)
		outFile.write(newLine + "\n")
	inFile.close()
	outFile.close()

	if os.path.isfile("./tmp/" + species + ".cgi.ensembl.txt"):
		inFile = file("./tmp/" + species + ".cgi.ensembl.txt")
		outFile = file("./tmp/" + species + ".cgi.ensembl.biotype.txt", "w")
		while True:
			line = inFile.readline()
			if len(line) == 0:
				break;
			line = string.strip(line)
			lineL = string.split(line, sep="\t")
			if lineL[0].find("_") > 0 or lineL[0].find(".") > 0:
				continue
			if lineL[3] in annoDict:
				lineL.append(annoDict[lineL[3]])
			else:
				lineL.append("No_anno")
			newLine = "\t".join(lineL)
			outFile.write(newLine + "\n")
		inFile.close()
		outFile.close()


if __name__ == '__main__':
	main()
