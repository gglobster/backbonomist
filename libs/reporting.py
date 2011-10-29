from config import directories as dirs, genomes, segtype
from common import ensure_dir
from drawing import ContigDraw, PairwiseDraw
import numpy as np

def map_scaffolds(contig):
    """Generate annotated maps of backbone scaffolds."""
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    constructs_root = dirs['constructs_dir']+ctg_name+"/"
    maps_out_root = dirs['maps_dir']+ctg_name+"/"
    ensure_dir(maps_out_root)
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        cstrct_gbk = constructs_root+g_name+"_"+ctg_name+"_cstrct.gbk"
        map_file = maps_out_root+g_name+"_"+ctg_name+"_annot.pdf"
        print '\t', g_name
        # generate graphical map
        ContigDraw(g_name, cstrct_gbk, map_file)

def map_pairwise(contig):
    """Generate pairwise alignment maps."""
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    ref_ctg_file = dirs['ref_ctg_dir']+contig['file']
    constructs_root = dirs['constructs_dir']+ctg_name+"/"
    segments_root = dirs['segments_dir']+ctg_name+"/"
    maps_out_root = dirs['maps_dir']+ctg_name+"/"
    ensure_dir(maps_out_root)
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        cstrct_gbk = constructs_root+g_name+"_"+ctg_name+"_cstrct.gbk"
        segments_file = segments_root+g_name+"_"+ctg_name+"_segs.txt"
        map_file = maps_out_root+g_name+"_"+ctg_name+"_aln.pdf"
        print '\t', g_name
        # load segments TODO: add idp-based clumping
        segdata = np.loadtxt(segments_file, skiprows=1, dtype=segtype)
        # generate graphical map
        PairwiseDraw(g_name, cstrct_gbk, ref_ctg_file, segdata, map_file)