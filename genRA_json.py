#! /usr/bin/env python

import os
import sys
import shutil
import regionanalysis
import json

def loadJSON(json_file):
    fp = open(json_file)
    json_info = json.load(fp)
    fp.close()
    return json_info

def main():
    genome = sys.argv[1]
    NPVer = sys.argv[2]
    species = sys.argv[3]
    assembly = sys.argv[4]
    ENSVer = sys.argv[5]
    RA_path = sys.argv[6]
    temp = loadJSON("json/RA_template.json")
    temp["genome"] = genome
    temp["version"] = NPVer
    temp["species"] = species
    temp["assembly"] = assembly
    for anno_db in temp["databases"]:
        if anno_db["database"] == "refseq":
            anno_db["version"] = genome
        if anno_db["database"] == "ensembl":
            anno_db["version"] = ENSVer
    with open(os.path.join(RA_path, genome+".json"), "w") as genome_json:
        json.dump(temp, genome_json, sort_keys=True, indent=2)

#-------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#-------------------------------------------------------------------------
# EOF
