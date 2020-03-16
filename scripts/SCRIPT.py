#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import cairo
import colorsys
import random
import re
import numpy as np

parser = argparse.ArgumentParser(description="This is a prfrom_fileram for determining binding motifs for transcription factors RBFOX and MBNL (format is intron-exon-intron) and outputs to-scale images of motif locations.")

parser.add_argument("-f", "--fasta", help="Input FASTA file containing inrons and exons.", required=True, type=str)
parser.add_argument("-m", "--motifs", help="Text file containing potential motifs.", required=True, type=str)

args = parser.parse_args()

fasta = args.fasta
motif_file = args.motifs


def GET_EXON(length, exon_start, exon_length, y):
    '''Draws an exon based off of its starting position and length'''

    start = 100
    exon = start + exon_start
    exon2 = start + exon_length
    context.set_line_width(3)
    context.move_to(start,y)
    context.line_to(length+start,y)
    context.stroke()
    context.set_line_width(12)
    context.move_to(exon,y)
    context.line_to(exon2,y)
    context.stroke()

def GET_REGEX(file):
    '''Converts a motif file to a regex format'''
    motifs_regex_list = []
    with open(file, "r") as fh:
        for line in fh:
            motif = ""
            line = line.strip("\n").upper()
            for char in line:
                motif += IUPAC[char]
            motifs_regex_list.append(motif)
    return motifs_regex_list

def MOTIFS(file):
    '''Makes a motif dictionary'''
    motifs = {}
    regex = []
    with open(file, "r") as fh:
        for line in fh:
            motif = ""
            line = line.strip("\n").upper()
            for char in line:
                motif += IUPAC[char]
            regex.append(motif)
            motifs[motif] = line
    return motifs

def GET_LEN(file):
    '''Outputs sequence length from a given FASTA record'''
    length = {}
    with open(fasta) as file:
        for record in SeqIO.parse(file, "fasta"):
            record = ("{},{}".format(record.id, record.seq))
            gene_name = record.split(",")[0]
            seq = record.split(",")[1]
            length[gene_name] = len(seq)
    return length

# motif_file = '/Users/abbiefayeolson/MOTIF_MARK/Fig_1_motifs.txt'
# fasta = '/Users/abbiefayeolson/MOTIF_MARK/Figure_1.fasta'

# MAIN SCRIPT

IUPAC = {"A": "[Aa]",
        "T":"[TtUu]",
        "C":"[Cc]",
        "G":"[Gg]",
        "U":"[UuTt]",
        "R": "[AaGg]",
        "Y":"[TtCcUu]",
        "S":"CcGg",
        "W":"AaTtUu",
        "K":"[GgTtUu]",
        "M":"[AaCc]",
        "B":"[CcGgTtUu]",
        "D":"[AaGgTtUu]",
        "H":"[AaCcTtUu]",
        "V":"[AaCcGg]",
        "N":"[AaTtCcGgUu]"}

COLOR = (138,43,226), (210,105,30), (124,252,0), (255,99,71)

# Iterates through a given FASTA file and initializes dictionaries
gene_length = GET_LEN(MOTIFS)
record_loc = {}
exon_pos = {}
with open(fasta) as file:
    for record in SeqIO.parse(file, "fasta"):
        record = ("{},{}".format(record.id, record.seq))
        gene_name = record.split(",")[0]
        seq = record.split(",")[1]
        seq_split = []
        exon = re.search('[A-Z]', seq)
        seq_split.append(seq[0:exon.start()])
        intron2 = re.search('[a-z]', seq[exon.start():])
        seq_split.append(seq[exon.start():exon.start()+intron2.start()])
        seq_split.append(seq[exon.start()+intron2.start():])
        record_loc[gene_name] = seq_split
        exon_pos[gene_name] = exon.start(), exon.start()+intron2.start()

# Records exon positions
dict = {}
motifs_regex_list = GET_REGEX(motif_file)
motifs_from_file = MOTIFS(motif_file)

for gene in record_loc.keys():
    seq = "".join(record_loc[gene])
    for motif in motifs_regex_list:
        loc = re.finditer(motif, seq)
        dict[gene+"_"+motifs_from_file[motif]] = []
        for l in loc:
            dict[gene+"_"+motifs_from_file[motif]].append(l.start())

# Assigns a color to each motif
col_dict = {}
with open(motif_file, "r") as fh:
        no = 0
        for line in fh:
            motif = ""
            line = line.strip("\n").upper()
            col_dict[line] = COLOR[no]
            no += 1

# Create a variable for scaling
max_gene = gene_length[max(gene_length, key=lambda i: gene_length[i])]
surface = cairo.SVGSurface("motifs.svg", max_gene*1.2, len(record_loc)*60+50)   ## width x height, use max_gene gene as setting for width, use no of genes as setting for height
context = cairo.Context(surface)

# Draw the gene in proportion to the length and exon
pixels= 25
right = 100
for i in gene_length.keys():
    context.set_source_rgb(0,0,0)
    context.move_to(50, pixels + 4)
    context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
    context.set_font_size(15)
    context.show_text(i)
    GET_EXON(gene_length[i], exon_pos[i][0], exon_pos[i][1], pixels)

    for from_file in motifs_from_file.values():
        m = i+"_"+from_file
        context.set_source_rgb(col_dict[from_file][0]/265, col_dict[from_file][1]/265, col_dict[from_file][2]/265)
        for a in range(len(dict[m])):
            context.set_line_width(12)
            context.move_to(dict[m][a]+100, pixels)
            context.line_to(dict[m][a]+100+len(from_file),pixels) ## make motifs proportional
            context.stroke()
    pixels+= 65

# Draw legend
for from_file in motifs_from_file.values():
        m = i+"_"+from_file
        context.set_source_rgb(col_dict[from_file][0]/265, col_dict[from_file][1]/265, col_dict[from_file][2]/265)
        context.move_to(right, len(gene_length)*65+20)
        context.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        #context.set_font_size(12)
        context.show_text(from_file)
        right += 100

surface.write_to_png("motifs.png")
surface.finish()
