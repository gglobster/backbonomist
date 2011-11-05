from os import listdir
from config import directories as dirs, p_root_dir, genomes, segtype
from common import ensure_dir
from array_tetris import offset_q2r_coords
from drawing import ContigDraw, PairwiseDraw
from loaders import load_genbank
import numpy as np

def prep_maps(ref_ctg, run_id):
    """Set up generation of various maps."""
    # set inputs and outputs
    ref_ctg_n = ref_ctg['name']
    run_root = p_root_dir+run_id+"/"
    cst_root = run_root+dirs['constructs_dir']+ref_ctg_n+"/"
    segments_root = run_root+dirs['aln_seg_dir']+ref_ctg_n+"/"
    ctgs_root = dirs['gbk_contigs_dir']
    ctg_segs_root = segments_root+"contigs/"
    cst_segs_root = segments_root+"constructs/"
    maps_root = run_root+dirs['maps_dir']+ref_ctg_n+"/"
    ctg_aln_maps_dir = maps_root+"contig_alns/"
    cst_ann_maps_dir = maps_root+"constructs_annot/"
    cst_aln_maps_dir = maps_root+"constructs_aln/"
    ensure_dir([cst_root, ctg_segs_root, ctg_aln_maps_dir,
                cst_ann_maps_dir, cst_aln_maps_dir])
    print " ", ref_ctg_n
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        print "\t", g_name, "...",
        cstrct_gbk = cst_root+g_name+"_"+ref_ctg_n+"_cstrct.gbk"
        # maps of contigs aligned to reference
        map_ctg_alns(ref_ctg, genome, ctgs_root, ctg_segs_root,
                     ctg_aln_maps_dir)
        print "ctg_aln",
        # map of scaffold construct
        map_cst_annot(ref_ctg, genome, cstrct_gbk, cst_ann_maps_dir)
        print "cst_annot",
        # map of construct aligned to reference
        map_cst_aln(ref_ctg, genome, cstrct_gbk, cst_segs_root,
                    cst_aln_maps_dir)
        print "cst_aln"

def map_ctg_alns(ref_ctg, genome, ctgs_root, ctg_segs_root, ctg_aln_maps_dir):
    """Generate maps of contigs aligned to reference."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = ref_ctg['name']
    maps_root = ctg_aln_maps_dir+g_name+"/"
    segs_root = ctg_segs_root+g_name+"/"
    ctgs_dir = ctgs_root+g_name+"/"
    ref_ctg_file = dirs['ori_g_dir']+ref_ctg['file']
    # list genbank files in matches directory
    dir_contents = listdir(segs_root)
    for ctg_num in dir_contents:
        ctg_gbk = ctgs_dir+g_name+"_"+ctg_num+".gbk"
        seg_file = segs_root+ctg_num+"/"+ctg_num+"_"+ref_ctg_n+"_segs.txt"
        map_file = maps_root+ctg_num+"_ctg_aln2_"+ref_ctg_n+".pdf"
        # load segments TODO: add idp-based clumping
        segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
        # offset coordinates where desired
        g_offset = genome['offset']
        if g_offset[0] != 0 or g_offset[1] != 0:
            q_len = len(load_genbank(ctg_gbk).seq)
            segdata = offset_q2r_coords(segdata, q_len, g_offset)
        # determine whether to flip the query sequence (negative offset)
        if g_offset[1] < 0:
            q_invert = True
        else:
            q_invert = False
        # generate graphical map
        PairwiseDraw(ref_ctg_n, g_name, ctg_gbk, ref_ctg_file, segdata,
                     map_file, q_invert, g_offset)

def map_cst_annot(ref_ctg, genome, cstrct_gbk, cst_ann_maps_dir):
    """Generate map of annotated scaffold construct."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = ref_ctg['name']
    map_file = cst_ann_maps_dir+g_name+"_cst_m2_"+ref_ctg_n+".pdf"
    ContigDraw(g_name, cstrct_gbk, map_file)

def map_cst_aln(ref_ctg, genome, cstrct_gbk, segs_root, cst_aln_maps_dir):
    """Generate map of construct aligned to reference."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = ref_ctg['name'] 
    ref_ctg_file = dirs['ori_g_dir']+ref_ctg['file']
    seg_file = segs_root+g_name+"/"+g_name+"_"+ref_ctg_n+"_segs.txt"
    map_file = cst_aln_maps_dir+g_name+"_cst_aln2_"+ref_ctg_n+".pdf"
    # load segments TODO: add idp-based clumping
    segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
    # offset coordinates where desired
    g_offset = genome['offset']
    if g_offset[0] != 0 or g_offset[1] != 0:
        q_len = len(load_genbank(cstrct_gbk).seq)
        segdata = offset_q2r_coords(segdata, q_len, g_offset)
    # determine whether to flip the query sequence (negative offset)
    if g_offset[1] < 0:
        q_invert = True
    else:
        q_invert = False
    # generate graphical map
    PairwiseDraw(ref_ctg_n, g_name, cstrct_gbk, ref_ctg_file, segdata,
                 map_file, q_invert, g_offset)