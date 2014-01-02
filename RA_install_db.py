#! /usr/bin/env python

import os
import sys
import shutil
import regionanalysis

def main():
    genome = sys.argv[1]
    module_dir = os.path.dirname(os.path.realpath(regionanalysis.__file__))
    db_path = os.path.join(module_dir, "database/")
    origin_path = "./spec_anno/%s"%genome
    files = [".genome", "_pericentromere.bed", "_subtelomere.bed", "_geneDesert.bed", \
    ".ensembl.biotype_region_ext.bed", ".refseq.biotype_region_ext.bed"]
    for i in files:
        try:
            src_file = os.path.join(origin_path, genome+i)
            db_file = os.path.join(db_path, genome+i)
            shutil.copyfile(src_file, db_file)
        except:
            ## if no annotation in UCSC, just skip it.
            if i == ".refseq.biotype_region_ext.bed":
                continue
            else:
                sys.stderr.write("No %s%s file in annotation folder!\n"%(genome, i))


#-------------------------------------------------------------------------
if __name__ == '__main__':
    main()
#-------------------------------------------------------------------------
# EOF
