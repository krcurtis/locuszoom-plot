# Copyright 2019 Fred Hutchinson Cancer Research Center
################################################################################
### Helper functions to generate LD info using PLINK and the 1000G
### data provided with LocusZoom 1.4



import os
import sys

from parallel_task_database.liteweight_worker import invoke_system as invoke_system


################################################################################
### MAGIC ENVIRONMENT MODULE STUFF
### Automatically load plink environment module when this module is loaded

sys.path.insert(0,"/app/Lmod/lmod/lmod/init")
from env_modules_python import module
module('load','plink')



################################################################################


LOCUSZOOM_1000G_DIR = "/fh/fast/lampe_j/Gut_Bugs/keith-working/locuszoom-1.4/locuszoom/data/1000G/genotypes/2014-10-14" # EUR

# ancestries_to_check = [ "AFR", "ASN", "AMR", "EUR"]

################################################################################





def generate_plink_ld_file(output_file, ancestry, chromosome_text, target_variant, window_kb=1000):
    known_sets = set(["AFR", "ASN", "AMR", "EUR", "SAN"])

    if ancestry not in known_sets:
        raise Exception("ERROR unknown ancestry: {a}".format(a=ancestry))

    params = [ "plink",
               "--bfile", os.path.join(LOCUSZOOM_1000G_DIR, ancestry, chromosome_text),
               "--r2",
               "--ld-snp", target_variant,
               "--ld-window-kb", str(window_kb),
               "--ld-window", "99999",
               "--ld-window-r2", "0",
               "--out", output_file
    ]


    invoke_system(params)

    plink_out = output_file + ".ld"
    plink_log = output_file + ".log"
    plink_nosex = output_file + ".nosex"

    os.unlink(plink_log)
    os.unlink(plink_nosex)
    os.rename(plink_out, output_file)


if "__main__" == __name__:
    generate_plink_ld_file("test.ld", "EUR", "chr11", "chr11:74898804")
