# Copyright 2019 Fred Hutchinson Cancer Research Center
################################################################################
### Plot p-value and R2 of region




import re

import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.cm as cmx
from matplotlib.patches import Rectangle


################################################################################

LD_REGIMES =  [  (0,   0.2, (0,0,0x80)),
                 (0.2, 0.4, (0x87,0xce, 0xfa)),
                 (0.4, 0.6, (0x00,0xff, 0x00)),
                 (0.6, 0.8, (0xff,0xa5, 0x00)),
                 (0.8, 1.0, (0xff,0x00, 0x00)),
]



################################################################################



def load_custom_pvalue_file(filename, threshold=1):
    frame = pd.read_csv(filename)
    S = frame['simple_pvalue'] < threshold
    reduced_frame = frame[S]
    reduced_frame = reduced_frame[ ["chromo", "position", "simple_pvalue"]]
    replace = { "chromo" : "chrom",
                "simple_pvalue" : "pvalue" }

    reduced_frame = reduced_frame.rename(columns=replace)
    reduced_frame['variant'] = [ 'chr{chrom}:{pos}'.format(chrom=c, pos=p) for c,p in zip(reduced_frame['chrom'], reduced_frame["position"])]

    return reduced_frame
    

def load_and_format_pvalue_file_custom(filename, target_variant, position_min, position_max, threshold=1):
    frame = pd.read_csv(filename)
    S = frame['simple_pvalue'] < threshold
    reduced_frame = frame[S]
    reduced_frame = reduced_frame[ ["chromo", "position", "simple_pvalue"]]
    replace = { "chromo" : "chrom",
                "simple_pvalue" : "pvalue" }

    reduced_frame = reduced_frame.rename(columns=replace)
    reduced_frame['variant'] = [ 'chr{chrom}:{pos}'.format(chrom=c, pos=p) for c,p in zip(reduced_frame['chrom'], reduced_frame["position"])]

    exclude_select = (reduced_frame['position'] < position_min) | (position_max < reduced_frame['position'])
    results = reduced_frame[~exclude_select]
    return results
    

def window_pvalue(pvalue_frame, position_min, position_max):
    exclude_select = (pvalue_frame['position'] < position_min) | (position_max < pvalue_frame['position'])
    results = pvalue_frame[~exclude_select]
    return results
    




def load_plink_r2_results_file(filename, target_variant, position_min, position_max):
    # PLINK has strange formatting for its output files
    lines = [line.strip(" \t\n") for line in open(filename)]
    records = [ re.split("[ \t]+", line) for line in lines]
    header = records[0]
    rest = records[1:]
    initial = pd.DataFrame.from_records(rest, columns=header)
    replace = { "SNP_A" : "target_variant",
                "CHR_B" : "chrom",
                "BP_B"  : "position",
                "SNP_B" : "variant",
                "R2"    : "ld_r2" }
    initial = initial.rename(columns=replace)
    initial['ld_r2'] = [ float(x) for x in initial['ld_r2']]
    results = initial[["target_variant", "chrom", "position", "variant", "ld_r2"]]
    return results



def merge_pvalue_ld(pvalue_frame, ld_frame):
    results = pd.merge(left=pvalue_frame, right=ld_frame, on="variant", how="left")
    results = results.rename(columns={"position_x": "position"})
    return results



def colorbar_magic(figure, n_plots):
    colors = []
    ncolors = len(LD_REGIMES)
    for _, _, c in LD_REGIMES:
        red, green, blue = c
        colors.append( [float(red)/255, float(green)/255, float(blue)/255, 1])

    custom_color_map = matplotlib.colors.ListedColormap(colors)

    cNorm  = matplotlib.colors.Normalize(vmin=0, vmax=1)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=custom_color_map)
    scalarMap.set_array([])

    #pyplot.gca().add_patch(Rectangle((0.1, 45), 40, 55, edgecolor='gray',
    #                                        linewidth=3, fill=False))

    if 1 == n_plots:
        positioning_top = 0.67 + 0.15
    else:
        positioning_top = 0.67 + 0.15+ 0.04
    positioning_magic = [0.15, positioning_top - 0.15/n_plots, 0.015, 0.16/n_plots] ## this position seems a bit magical and is probably liable to break
    fontsize_magic = 7 # 8 works for one plot but seems perhaps a little big

    cax = figure.add_axes(positioning_magic)  
    cax.set_title("r2", fontsize=fontsize_magic)  

    pyplot.colorbar(scalarMap, cax = cax, ticks=[float(i) / ncolors for i in range(0, ncolors+1)],
                                            orientation='vertical')
    for text_obj in cax.get_yticklabels():
        text_obj.set_fontsize(fontsize_magic)





def plot_r2_region_worker(r2_axes, pvalue_ld_frame, target_variant, position_min, position_max, fancy_variant_name):
    M = 1e6
    position_min = position_min/M
    position_max = position_max/M

    ld_color_groups = {} ## will include "gray" for unknown ld

    colors = ["grey"]
    for ld_min, ld_max, color_tuple in LD_REGIMES:
        red, green, blue = color_tuple
        color = "#{r:02X}{g:02X}{b:02X}".format(r=red, g=green, b=blue) 
        colors.append(color)
        x = []
        y = []
        for pos,pvalue,ld,variant in zip(pvalue_ld_frame['position'], pvalue_ld_frame['pvalue'], pvalue_ld_frame['ld_r2'], pvalue_ld_frame['variant']):
            if pd.isnull(ld):
                continue

            elif variant == target_variant:
                continue

            elif ld_min <= ld and ld < ld_max:
                x.append(pos)
                y.append(-np.log10(pvalue))

        ld_color_groups[color] = (x,y)

    x = []
    y = []
    for pos,pvalue,ld,variant in zip(pvalue_ld_frame['position'], pvalue_ld_frame['pvalue'], pvalue_ld_frame['ld_r2'], pvalue_ld_frame['variant']):
        if pd.isnull(ld) and variant != target_variant:
            x.append(pos)
            y.append(-np.log10(pvalue))

    ld_color_groups["grey"] = (x,y)

    variant_pos = None
    variant_y = None
    for pos,pvalue,variant in zip(pvalue_ld_frame['position'], pvalue_ld_frame['pvalue'], pvalue_ld_frame['variant']):
        if variant == target_variant:
            variant_pos = pos
            variant_y = -np.log10(pvalue)



    #scale and plot
    variant_pos /=M
    for color in colors:
        x,y = ld_color_groups[color]
        r2_axes.plot(np.array(x)/M, y,'.', color=color)

    r2_axes.plot([variant_pos], [variant_y], 'D', color="purple")
    #pyplot.text(variant_pos, variant_y, fancy_variant_name, horizontalalignment='center', verticalalignment='bottom')

    # annotate on side
    r2_axes.annotate(fancy_variant_name, (variant_pos, variant_y),
                     horizontalalignment='left', verticalalignment='center', xytext=(6, 0), textcoords="offset points")

    # annotate to top
    #r2_axes.annotate(fancy_variant_name, (variant_pos, variant_y),
    #                 horizontalalignment='center', verticalalignment='bottom', xytext=(0, 6), textcoords="offset points")


    r2_axes.set_xlim(position_min, position_max)
    # r2_axes.set_xlabel('position (nt)')
    r2_axes.set_xticks([])
    r2_axes.set_ylabel('-log10(p-value)')

    #colorbar_magic(mainfig)
    #colorbar_magic(r2_axes)



def plot_r2_region(output_plot, pvalue_ld_frame, target_variant, position_min, position_max, fancy_variant_name):
    pyplot.clf()
    mainfig, r2_axes = pyplot.subplots(1,1, num=None, figsize=(8,6), dpi=150)

    plot_r2_region_worker(r2_axes, pvalue_ld_frame, target_variant, position_min, position_max, fancy_variant_name)

    #pyplot.show()
    if None != output_plot:
        pyplot.savefig(output_plot, bbox_inches='tight');
        pyplot.draw()
    else:
        pyplot.show()
    pyplot.close(mainfig)



