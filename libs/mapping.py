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
    ref_gbk_root = run_root+dirs['ref_gbk_dir']
    ref_gbk = ref_gbk_root+ref_ctg_n+"_re-annot.gbk"
    cst_root = run_root+dirs['scaffolds_dir']+ref_ctg_n+"/"
    segments_root = run_root+dirs['aln_seg_dir']+ref_ctg_n+"/"
    ctg_segs_root = segments_root+"contigs/"
    cst_segs_root = segments_root+"constructs/"
    maps_root = run_root+dirs['maps_dir']+ref_ctg_n+"/"
    ctg_aln_maps_root = maps_root+"contig_alns/"
    cst_ann_maps_root = maps_root+"constructs_annot/"
    cst_aln_maps_root = maps_root+"constructs_aln/"
    ensure_dir([cst_root, ctg_segs_root, cst_segs_root, maps_root,
                ctg_aln_maps_root, cst_ann_maps_root, cst_aln_maps_root])
    print " ", ref_ctg_n, "...",
    # map of reference with segment details
    map_ref_segs(ref_ctg, run_id)
    print "ref_map"
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        print "\t", g_name, "...",
        scaff_gbk = cst_root+g_name+"_"+ref_ctg_n+"_scaffold.gbk"
        ctg_aln_maps_dir = ctg_aln_maps_root+g_name+"/"
        ensure_dir([ctg_aln_maps_dir])
        # maps of contigs aligned to reference
        print "ctg_aln",
        map_ctg_alns(ref_ctg, ref_gbk, genome, ctg_segs_root,
                     ctg_aln_maps_dir)
        # map of scaffold construct
        print "cst_annot",
        map_cst_annot(ref_ctg, genome, scaff_gbk, cst_ann_maps_root)
        # map of construct aligned to reference
        print "cst_aln"
        map_cst_aln(ref_ctg, ref_gbk, genome, scaff_gbk, cst_segs_root,
                    cst_aln_maps_root)

def map_ctg_alns(ref_ctg, ref_gbk, genome, ctg_segs_root, maps_root):
    """Generate maps of contigs aligned to reference."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = ref_ctg['name']
    segs_root = ctg_segs_root+g_name+"/"
    ctgs_dir = dirs['gbk_contigs_dir']+g_name+"/"
    # list genbank files in matches directory
    dir_contents = listdir(segs_root)
    for ctg_num in dir_contents:
        ctg_gbk = ctgs_dir+g_name+"_"+ctg_num+".gbk"
        seg_file = segs_root+ctg_num+"/"+ctg_num+"_"+ref_ctg_n+"_segs.txt"
        map_file = maps_root+g_name+"_"+ctg_num+"_vs_"+ref_ctg_n+".pdf"
        # start mapping
        try:
            # load segments TODO: add idp-based clumping
            segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
            # deactivate offsetting
            g_offset = (0,0)
            q_invert = False
            # generate graphical map
            PairwiseDraw(ref_ctg_n, g_name+"_"+ctg_num, ref_gbk, ctg_gbk,
                         segdata, map_file, q_invert, g_offset, 'dual',
                         'dual', 'm', 'fct', 'fct')
        except IOError:
            print "\nERROR: could not load segments data for", g_name, ctg_num
            print "\t\t\t",
        except StopIteration:
            print "\nERROR: could not make map for", g_name, ctg_num
            print "\t\t\t",

def map_cst_annot(ref_ctg, genome, scaff_gbk, maps_root):
    """Generate map of annotated scaffold construct."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = ref_ctg['name']
    map_file = maps_root+g_name+"_<"+ref_ctg_n+".pdf"
    # start mapping
    try: open(scaff_gbk, 'r')
    except IOError:
        print "WARNING: No scaffold construct to map"
    else:
        ContigDraw(g_name, scaff_gbk, map_file, 'm', 'fct')

def map_cst_aln(ref_ctg, ref_gbk, genome, scaff_gbk, segs_root, maps_root):
    """Generate map of construct aligned to reference."""
    # set inputs and outputs
    g_name = genome['name']
    ref_ctg_n = ref_ctg['name']
    seg_file = segs_root+g_name+"/"+g_name+"_"+ref_ctg_n+"_segs.txt"
    map_file = maps_root+g_name+"_vs_"+ref_ctg_n+".pdf"
    # start mapping
    try: open(scaff_gbk)
    except IOError:
        print "WARNING: No scaffold construct to map"
    else:
        try:
            # load segments TODO: add idp-based clumping
            segdata = np.loadtxt(seg_file, skiprows=1, dtype=segtype)
        except IOError:
            print "\nERROR: could not load segments data for", g_name
        except StopIteration:
            print "\nERROR: could not make map for", g_name, "construct"
        else:
            # offset coordinates where desired
            g_offset = genome['offset']
            if g_offset[0] != 0 or g_offset[1] != 0:
                q_len = len(load_genbank(scaff_gbk).seq)
                segdata = offset_q2r_coords(segdata, q_len, g_offset)
            # determine whether to flip the query sequence (negative offset)
            if g_offset[1] < 0:
                q_invert = True
            else:
                q_invert = False
            # generate graphical map
            PairwiseDraw(ref_ctg_n, g_name, ref_gbk, scaff_gbk, segdata,
                         map_file, q_invert, g_offset, 'dual', 'dual', 'm',
                         'fct', 'fct')

def map_ref_segs(ref_ctg, run_id):
    """Generate map of reference contig with segment details.

    This provides a comparison of the original reference and the
    re-annotated version.

    """
    # set inputs and outputs
    ref_ctg_n = ref_ctg['name']
    run_root = p_root_dir+run_id+"/"
    ori_file = dirs['ori_g_dir']+ref_ctg['file']
    ref_maps_root = run_root+dirs['ref_map_dir']
    ensure_dir([ref_maps_root])
    gbk_file = run_root+dirs['ref_gbk_dir']+ref_ctg_n+"_re-annot.gbk"
    map_file = ref_maps_root+ref_ctg_n+"_ref.pdf"
    # start mapping
    try:
        # make mock segment, full-length with 100% id
        record = load_genbank(gbk_file)
        length = len(record.seq)
        segdata = [[1, length, 1, length, 100]]
        # deactivate offsetting
        g_offset = (0,0)
        q_invert = False
        # generate graphical map
        PairwiseDraw(ref_ctg_n+"_ra", ref_ctg_n+"_ori", gbk_file, ori_file,
                     segdata, map_file, q_invert, g_offset, 'dual', 'dual',
                     'm', 'fct', 'product')
    except IOError:
        print "\nERROR: could not load segments data for", ref_ctg_n
        print "\t\t\t",
    except StopIteration:
        print "\nERROR: could not make map for", ref_ctg_n
        print "\t\t\t",
