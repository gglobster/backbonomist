import os, subprocess, re
from os import listdir
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from config import mauve_exec, directories as dirs, p_root_dir, genomes, \
    max_size, chop_mode
from common import ensure_dir
from parsing import mauver_load2_k0, parse_clustal_idstars
from array_tetris import chop_rows
from loaders import load_genbank, load_fasta
from writers import write_fasta

def align_clustal(file_name):
    """Make external call to ClustalW aligner."""
    cline = ClustalwCommandline(infile=file_name)
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse ClustalW errors
    return report

def align_muscle(infile_name, outfile_name, log_file):
    """Make external call to Muscle aligner."""
    cline = MuscleCommandline(input=infile_name, out=outfile_name, clw=True,
                              loga=log_file, quiet='y')
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse MUSCLE errors
    return report

def align_mauve(file_list, output):
    """Make external call to Mauve aligner."""
    input_files = ' '.join(file_list)
    cline = mauve_exec+" --output="+ output +" "+input_files
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse Mauve errors
    return report

def align_ctg2ref(contig, run_id):
    """Align contigs pairwise to the reference contig."""
    # set inputs and outputs
    ref_ctg_name = contig['name'] # reference contig
    run_root = p_root_dir+run_id+"/"
    ref_ctg_file = dirs['ori_g_dir']+contig['file']
    mauve_root = run_root+dirs['mauve_out_dir']+ref_ctg_name+"/contigs/"
    segments_root = run_root+dirs['aln_seg_dir']+ref_ctg_name+"/contigs/"
    q_ctgs_root = run_root+dirs['match_out_dir']+ref_ctg_name+"/"
    ensure_dir([segments_root])
    print " ", ref_ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        print "\t", g_name, "...",
        ctgs_fas_dir = q_ctgs_root+g_name+"/"
        mauve_dir = mauve_root+g_name+"/"
        aln_segs_root = segments_root+g_name+"/"
        ensure_dir([mauve_dir])
        # list genbank files in matches directory
        dir_contents = listdir(ctgs_fas_dir)
        for item in dir_contents:
            pattern = re.compile(r'.*_(\d*)\.fas$')
            match = pattern.match(item)
            if match:
                ctg_num = match.group(1)
                print ctg_num,
                # set inputs and outputs
                q_contig = ctgs_fas_dir+item
                file_list = (ref_ctg_file, q_contig)
                mauve_outfile = mauve_dir+ctg_num+".mauve"
                aln_segs_dir = aln_segs_root+ctg_num+"/"
                ensure_dir([aln_segs_dir])
                segfile = aln_segs_dir+ctg_num+"_"+ref_ctg_name+"_segs.txt"
                open(segfile, 'w').write('')
                # do Mauve alignment
                align_mauve(file_list, mauve_outfile)
                try:
                    # parse Mauve output (without initial clumping)
                    coords = mauver_load2_k0(mauve_outfile+".backbone", 0)
                    # chop segments that are too long
                    chop_array = chop_rows(coords, max_size, chop_mode)
                    # make detailed pairwise alignments of the segments
                    ref_rec = load_genbank(ref_ctg_file)
                    query_rec = load_fasta(q_contig)
                    iter_align(chop_array, ref_rec, query_rec, aln_segs_dir,
                               segfile)
                except IOError:
                    print "\nERROR: Mauve alignment failed"
                    print "\t\t\t",
        print ""

def align_cstrct2ref(contig, run_id):
    """Align constructs pairwise to the reference contig."""
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    run_root = p_root_dir+run_id+"/"
    ref_ctg_file = dirs['ori_g_dir']+contig['file']
    mauve_root = run_root+dirs['mauve_out_dir']+ctg_name+"/constructs/"
    segments_root = run_root+dirs['aln_seg_dir']+ctg_name+"/constructs/"
    scaff_root = run_root+dirs['scaffolds_dir']+ctg_name+"/"
    ensure_dir([segments_root])
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        scaff_gbk = scaff_root+g_name+"_"+ctg_name+"_scaffold.gbk"
        file_list = (ref_ctg_file, scaff_gbk)
        print "\t", g_name, "...",
        # set outputs
        mauve_dir = mauve_root+g_name+"/"
        aln_segs_dir = segments_root+g_name+"/"
        ensure_dir([mauve_dir, aln_segs_dir])
        mauve_outfile = mauve_dir+g_name+"_"+ctg_name+".mauve"
        segfile = aln_segs_dir+g_name+"_"+ctg_name+"_segs.txt"
        # abort if there is no scaffold construct
        try: open(scaff_gbk, 'r')
        except IOError:
            print "WARNING: No scaffold construct to align"
        else:
            # prep segments file
            open(segfile, 'w').write('')
            # purge any pre-existing sslist file
            sslist_file = scaff_gbk+".sslist"
            if os.path.isfile(sslist_file):
                try: os.remove(sslist_file)
                except Exception: raise
            # do Mauve alignment
            align_mauve(file_list, mauve_outfile)
            try:
                # parse Mauve output (without initial clumping)
                coords = mauver_load2_k0(mauve_outfile+".backbone", 0)
                print len(coords), '->',
                # chop segments that are too long
                chop_array = chop_rows(coords, max_size, chop_mode)
                print len(chop_array), 'segments <', max_size, 'bp',
                # make detailed pairwise alignments of the segments
                ref_rec = load_genbank(ref_ctg_file)
                query_rec = load_genbank(scaff_gbk)
                id = iter_align(chop_array, ref_rec, query_rec, aln_segs_dir, segfile)
                print "@", id, "% id. overall"
            except IOError:
                print "\nERROR: Mauve alignment failed"

def iter_align(coord_array, ref_rec, query_rec, aln_dir, segs_file):
    """Iterate through array of coordinates to make pairwise alignments."""
    # set up the root subdirectories
    seqs = aln_dir+"input_seqs/"
    alns = aln_dir+"output_alns/"
    ensure_dir([seqs, alns])
    aln_id = 0
    aln_len = 0
    # cycle through segments
    for segment_pair in coord_array:
        xa, xb, xc, xd = segment_pair
        # extract the corresponding sequence slices
        ref_seq = ref_rec[abs(xa):abs(xb)]
        query_seq = query_rec[abs(xc):abs(xd)]
        # reverse-complement sequences with negative sign
        if xa < 0 :
            ref_seq = ref_seq.reverse_complement()
        if xc < 0 :
            query_seq = query_seq.reverse_complement()
        # write sequences to file
        mscl_in = seqs+str(xa)+"_"+str(xb)+"_"+str(xc)+"_"+str(xd)+".fas"
        write_fasta(mscl_in, [ref_seq, query_seq])
        # skip segments that are too small to align
        if abs(abs(xa)-abs(xb)) < 10:
            idp = 0
        else:
            # set up outfiles
            mscl_out = alns+str(xa)+"_"+str(xb)+"_"+str(xc)+"_"+str(xd)+".aln"
            logfile = aln_dir+"muscle_log.txt"
            # perform alignment
            align_muscle(mscl_in, mscl_out, logfile)
            idntot = parse_clustal_idstars(mscl_out)
            idp = int((float(idntot)/len(query_seq))*100)
            aln_id += idntot
            aln_len += len(query_seq)
        # write details out to segments file
        line = "\t".join([str(xa), str(xb), str(xc), str(xd), str(idp)+"\n"])
        open(segs_file, 'a').write(line)
    overall_id = int((float(aln_id)/aln_len)*100)
    return overall_id
