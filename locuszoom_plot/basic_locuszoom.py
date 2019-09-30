# Copyright 2019 Fred Hutchinson Cancer Research Center
################################################################################
### Combine gene region and pvalue/R2 plot to make a basic locus zoom plot





import numpy as np
import matplotlib.pyplot as pyplot

from .plot_gene_region import load_gene_region_info
from .plot_gene_region import sort_gene_locations
from .plot_gene_region import plot_gene_region_worker


from .plot_r2_region import plot_r2_region_worker
from .plot_r2_region import colorbar_magic
from .plot_r2_region import load_and_format_pvalue_file_custom
from .plot_r2_region import load_plink_r2_results_file
from .plot_r2_region import merge_pvalue_ld




################################################################################


################################################################################


def basic_locuszoom(pvalue_file, plink_file, target_variant, target_pos, fancy_name, output_plot=None, output_pdf=None, title=None):

    target_chromosome = target_variant.split(":")[0]
    position_min = target_pos - 1000000
    position_max = target_pos + 1000000
    region_info = load_gene_region_info(target_chromosome, position_min, position_max)
    gene_rows = sort_gene_locations(region_info, position_min, position_max)



    pvalue_frame = load_and_format_pvalue_file_custom(pvalue_file, target_variant, position_min, position_max)
    ld_frame = load_plink_r2_results_file(plink_file, target_variant, position_min, position_max)
    pvalue_ld_result = merge_pvalue_ld(pvalue_frame, ld_frame)

    n_gene_rows = len(gene_rows)

    pyplot.clf()
    mainfig, (r2_axes, gene_axes) = pyplot.subplots(2,1, gridspec_kw={'height_ratios': [10, n_gene_rows]}, num=None, figsize=(8,6), dpi=150)    

    plot_r2_region_worker(r2_axes, pvalue_ld_result, target_variant, position_min, position_max, fancy_name)

    plot_gene_region_worker(gene_axes, gene_rows, position_min, position_max)

    colorbar_magic(mainfig)

    if None != title:
        mainfig.suptitle(title)

    if None != output_plot:
        pyplot.savefig(output_plot, bbox_inches='tight');
        pyplot.draw()
    elif None != output_pdf:
        output_pdf.savefig(mainfig)
    else:
        raise Exception("ERROR no output plot or output pdf object given")


