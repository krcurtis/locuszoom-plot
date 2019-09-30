# Copyright 2019 Fred Hutchinson Cancer Research Center
################################################################################
### Combine gene region and pvalue/R2 plots for different ancestries






import numpy as np
import matplotlib.pyplot as pyplot



from .plot_gene_region import load_gene_region_info
from .plot_gene_region import sort_gene_locations
from .plot_gene_region import plot_gene_region_worker


from .plot_r2_region import plot_r2_region_worker
from .plot_r2_region import colorbar_magic
from .plot_r2_region import load_plink_r2_results_file
from .plot_r2_region import merge_pvalue_ld
from .plot_r2_region import window_pvalue



################################################################################


################################################################################


def multi_ancestry_locuszoom(pvalue_frame, ancestry_file_set, target_variant, target_pos, fancy_name, output_plot=None, output_pdf=None, title=None):

    target_chromosome = target_variant.split(":")[0]
    position_min = target_pos - 1000000
    position_max = target_pos + 1000000
    region_info = load_gene_region_info(target_chromosome, position_min, position_max)
    gene_rows = sort_gene_locations(region_info, position_min, position_max)


    pvalue_frame = window_pvalue(pvalue_frame, position_min, position_max)

    n_gene_rows = len(gene_rows)

    n_ancestry = len(ancestry_file_set)
    height_ratios = {'height_ratios': [10 for i in range(n_ancestry)] + [n_gene_rows]}
    pyplot.clf()


    mainfig, axes_objects = pyplot.subplots(n_ancestry+1,1, gridspec_kw=height_ratios, num=None, figsize=(8,11), dpi=150)    

    gene_axes = axes_objects[n_ancestry]
    for i, ancestry_group in enumerate(ancestry_file_set):
        label, plink_file = ancestry_group

        ld_frame = load_plink_r2_results_file(plink_file, target_variant, position_min, position_max)
        pvalue_ld_result = merge_pvalue_ld(pvalue_frame, ld_frame)
        r2_axes = axes_objects[i]

        plot_r2_region_worker(r2_axes, pvalue_ld_result, target_variant, position_min, position_max, fancy_name)
        r2_axes.set_title(label)

    plot_gene_region_worker(gene_axes, gene_rows, position_min, position_max)

    colorbar_magic(mainfig, n_ancestry)

    #if None != title:
    #    mainfig.suptitle(title)

    if None != output_plot:
        pyplot.savefig(output_plot, bbox_inches='tight');
        pyplot.draw()
    elif None != output_pdf:
        output_pdf.savefig(mainfig)
    else:
        raise Exception("ERROR no output plot or output pdf object given")


