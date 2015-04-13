#! /usr/bin/env python
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
    p1 = re.compile(r'.*\S*; gene_name "(\S*;\S*)"; transcript_name.*')
    p2 = re.compile(r'.*\S*; transcript_name "(\S*;\S*)"; *')
    while True:
        line = inFile.readline()
        if len(line) == 0:
            break;
        line = string.strip(line)
        if line.startswith("#"):
            continue;
        annoPart = string.split(line, sep="\t")[-1][:-1].replace("\"", "")
        # annoPartL = string.split(annoPart, sep=";")
        source = line.split("\t")[1]  # LS.
        try:
            # geneDict = dict(item.split(" ")[1:3] for item in annoPart.split(";"))
            geneDict = dict(item.strip().split(" ")[:2] for item in annoPart.split(";"))  # LS.
        except:
            try:
                # to replace some ";" in gene_name!
                raw_name = p1.match(line).group(1)
                target_name = raw_name.replace(";", "-")
                line = string.replace(line, raw_name, target_name)
                annoPart = string.split(line, sep="\t")[-1][:-1].replace("\"", "")
                # annoPartL = string.split(annoPart, sep=";")
                geneDict = dict(item.split(" ")[1:3] for item in annoPart.split(";"))
            except:
                # replace some ";" in transcript_id!
                # like 'transcript_name "PHT4;7-202";' in Oryza_sativa, IRGSP-1
                raw_name = p2.match(line).group(1)
                target_name = raw_name.replace(";", "-")
                line = string.replace(line, raw_name, target_name)
                annoPart = string.split(line, sep="\t")[-1][:-1].replace("\"", "")
                # annoPartL = string.split(annoPart, sep=";")
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


if __name__ == '__main__':
    main()
