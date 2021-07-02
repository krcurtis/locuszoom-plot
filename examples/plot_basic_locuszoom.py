################################################################################
### Plot a basic locuszoom plot using the example data



import os

import pandas as pd
import numpy as np

import locuszoom_plot as lzp


################################################################################

OUTPUT_PLOT = "my_example_locuszoom.png"


ANCESTRY =  "EUR"

LOCUSZOOM_DIR = os.environ["LOCUSZOOM_DIR"]  # or set to where your locuszoom databases are located
LOCUSZOOM_GENE_DB = os.path.join(LOCUSZOOM_DIR, "locuszoom/data/database/locuszoom_hg19.db")
LOCUSZOOM_GENOTYPES_TEMPLATE = os.path.join(LOCUSZOOM_DIR, "locuszoom/data/1000G/genotypes/2014-10-14/{ancestry}/{chrom}")

SCRATCH_DIR = "/tmp"

PVALUE_FILE = "random_example_data.csv"

################################################################################

rs_name = "rs3197999"
target_variant = "chr3:49721532"  # from dbsnp search results
target_pos = 49721532
chromosome = 3

pvalue_frame = pd.read_csv(PVALUE_FILE)
ancestry="EUR"
scratch_ld_file = os.path.join(SCRATCH_DIR, target_variant + "_" + ancestry + ".ld")


lzp.generate_plink_ld_file(scratch_ld_file, ancestry, "chr{chrom}".format(chrom=chromosome), target_variant,
                           locuszoom_template=LOCUSZOOM_GENOTYPES_TEMPLATE)

lzp.basic_locuszoom(pvalue_frame, scratch_ld_file, target_variant, target_pos, rs_name, target_window_size=500000,
                    output_plot=OUTPUT_PLOT, output_pdf=None, title="Example with random data", locuszoom_gene_db=LOCUSZOOM_GENE_DB)





