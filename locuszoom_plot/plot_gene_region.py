# -*- coding: utf-8 -*-
# Copyright 2019 Fred Hutchinson Cancer Research Center
################################################################################
### test gridspec stuff for locuszoom plotting


import sqlite3

import pandas as pd
import matplotlib.pyplot as pyplot
import matplotlib.gridspec as gridspec

################################################################################



# LOCUSZOOM_GENE = "locuszoom1.4_data/locuszoom/data/database/locuszoom_hg19.db"


################################################################################


def overlap_interval(a_min, a_max, b_min, b_max):
    if a_max < b_min:
        return False
    elif b_max < a_min:
        return False
    return True

def overlap_region(a, b):
    return overlap_interval(a['txStart'], a['txEnd'], b['txStart'], b['txEnd'])


def overlap_text(a, b, position_min, position_max):
    text_scaling_factor = .015  ## this is a pure guesstimate here

    position_range = (position_max - position_min)

    center_a = (a['txStart'] + a['txEnd'])/2
    center_a_x = (center_a - position_min) / position_range

    a_width = text_scaling_factor * (1 + len(a['geneName']))

    center_b = (b['txStart'] + b['txEnd'])/2
    center_b_x = (center_b - position_min) / position_range
    b_width = text_scaling_factor * (1 + len(b['geneName']))

    return overlap_interval(center_a_x - a_width/2, center_a_x + a_width/2, 
                            center_b_x - b_width/2, center_b_x + b_width/2)


def sort_gene_locations(gene_frame, position_min, position_max):
    records = gene_frame.to_dict('records')

    for gene in records:
        gene['exonStarts'] = [ int(s) for s in gene['exonStarts'].split(',') if len(s) > 0]
        gene['exonEnds'] =  [ int(s) for s in gene['exonEnds'].split(',') if len(s) > 0]



    rows = []
    while len(records) > 0:

        current_row = []
        last_item = None
        remainder = []
        
        for r in records:
            if None == last_item:
                current_row.append(r)
                last_item = r
            elif overlap_region(last_item, r):
                remainder.append(r)
            elif overlap_text(last_item, r, position_min, position_max):
                remainder.append(r)
            else:
                current_row.append(r)
                last_item = r
        rows.append(current_row)
        records = remainder


    return rows






def load_gene_region_info(chromosome, position_min, position_max, locuszoom_gene_db):
    conn = sqlite3.connect(locuszoom_gene_db)

    cursor = conn.cursor()
    cursor.execute("SELECT * FROM refFlat where refFlat.chrom='{chrom}';".format(chrom=chromosome))
    records = [ r for r in cursor.fetchall()]
    columns = [ t[0] for t in cursor.description]
    frame = pd.DataFrame.from_records(records, columns=columns)
    exclude_select = (frame['txEnd'] < position_min) | (position_max < frame['txStart'])
    frame = frame[~exclude_select]
    frame = frame.drop_duplicates('geneName')
    # drop duplicate gene names
    return frame.sort_values(by="txStart")


def scale_gene_rows(gene_rows):
    """Scale gene rows to use megabases or mega-nucleotides rather than plain nucleotides"""

    M = 1e6
    for row in gene_rows:
        for gene in row:
            gene['txStart'] = float(gene['txStart']) / M
            gene['txEnd'] = float(gene['txEnd']) / M
            gene['exonStarts'] = [ float(s)/M for s in gene['exonStarts']]
            gene['exonEnds'] = [ float(s)/M for s in gene['exonEnds']]

    return gene_rows


def plot_gene_region_worker(gene_axes, gene_rows, position_min, position_max):
    M = 1e6
    position_min = position_min/ M
    position_max = position_max/M

    scale_gene_rows(gene_rows)

    text_yoffset = 0.4
    fontsize_magic = 5
    nrows = len(gene_rows)
    for i,row in enumerate(gene_rows):
        y_coord = -i
        for gene in row:
            center_x = (gene['txStart'] +  gene['txEnd'])/2
            x = [float(gene['txStart']), float(gene['txEnd'])]
            y = [y_coord, y_coord]
            gene_axes.plot(x,y, 'b-|', linewidth=1)   ## force vertical line at beginning and end for short genes that might otherwis get dropped by the plot rendering engine

            
            exonStarts = gene['exonStarts']
            exonEnds = gene['exonEnds']


            for es, ee in zip(exonStarts, exonEnds):
                x = [es, ee]
                y = [y_coord, y_coord]
                gene_axes.plot(x,y, 'b-|', linewidth=5)

            if '+' == gene['strand']:
                gene_axes.text(center_x, y_coord + text_yoffset, gene['geneName'] + "→",
                               horizontalalignment='center', verticalalignment='center',
                               fontsize=fontsize_magic)
            elif '-' == gene['strand']:
                gene_axes.text(center_x, y_coord + text_yoffset, "←" + gene['geneName'],
                               horizontalalignment='center', verticalalignment='center',
                               fontsize=fontsize_magic)


    gene_axes.set_ylim(-nrows+0.5, 1)
    gene_axes.set_yticks([])
    gene_axes.set_xlim(position_min, position_max)
    gene_axes.set_xlabel('position (Mb)')  ## or mega-nucelotides?

    gene_axes.spines['top'].set_visible(False)
    gene_axes.spines['right'].set_visible(False)
    gene_axes.spines['left'].set_visible(False)
    #gene_axes.spines['bottom'].set_visible(False)



def plot_gene_region(output_plot, gene_rows, position_min, position_max):
    pyplot.clf()
    mainfig, gene_axes = pyplot.subplots(1,1, num=None, figsize=(8,6), dpi=150)

    plot_gene_region_worker(gene_axes, gene_rows, position_min, position_max)

    #pyplot.show()
    if None != output_plot:
        pyplot.savefig(output_plot, bbox_inches='tight');
        pyplot.draw()
    else:
        pyplot.show()
    pyplot.close(mainfig)

