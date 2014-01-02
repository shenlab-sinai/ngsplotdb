import os
import sys
import string

def main():
	dbTxt = sys.argv[1]
	outFileName = sys.argv[2]
	species = sys.argv[3]

	inFile = file(dbTxt)
	outFile = file(outFileName, "w")
	while True:
		line = inFile.readline()
		if len(line) == 0:
			break;
		line = string.strip(line)
		lineL = string.split(line, sep="\t")
		if lineL[0].find("_") > 0 or lineL[0].find(".") > 0:
			continue
		gname = lineL[4]
		gid = lineL[5]
		if gname.startswith("LINC") and gid.startswith("NR"):
			lineL.append("lincRNA")
		elif gname.startswith("MIR") and gid.startswith("NR"):
			lineL.append("miRNA")
		elif gid.startswith("NM"):
			lineL.append("protein_coding")
		else:
			lineL.append("misc")
		newLine = "\t".join(lineL)
		outFile.write(newLine + "\n")
	inFile.close()
	outFile.close()

	if os.path.isfile("./tmp/" + species + ".cgi.refseq.txt"):
		inFile = file("./tmp/" + species + ".cgi.refseq.txt")
		outFile = file("./tmp/" + species + ".cgi.refseq.biotype.txt", "w")
		while True:
			line = inFile.readline()
			if len(line) == 0:
				break;
			line = string.strip(line)
			lineL = string.split(line, sep="\t")
			if lineL[0].find("_") > 0 or lineL[0].find(".") > 0:
				continue
			gname = lineL[3]
			gid = lineL[4]
			if gname.startswith("LINC") and gid.startswith("NR"):
				lineL.append("lincRNA")
			elif gname.startswith("MIR") and gid.startswith("NR"):
				lineL.append("miRNA")
			elif gid.startswith("NM"):
				lineL.append("protein_coding")
			elif len(gid)!=0:
				lineL.append("misc")
			else:
				lineL.append("No_anno")
			newLine = "\t".join(lineL)
			outFile.write(newLine + "\n")
		inFile.close()
		outFile.close()


if __name__ == '__main__':
	main()