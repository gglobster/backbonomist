import re
import numpy as np
from os import path, listdir
from loaders import load_genbank, load_multifasta
from writers import write_genbank, write_fasta
from string_ops import multisplit_finder
from common import ensure_dir
from config import separator, directories as dirs, p_root_dir, genomes, prox_D
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from parsing import mauver_load2_k0
from array_tetris import get_anchor_loc
from annotation import annot_ctgs

def unpack_genomes(genome):
    """Unpack genome files.

    Here, unpacking means extracting data and producing specific files to
    standardize how the information is made available to downstream analysis.
    Depending on the input file format, different unpacking methods are
    invoked. In all cases, this ensures that for each genome, there is a
    multifasta file of the contigs all together as well as a separate Genbank
    file for each contig.

    Supported input file formats are the following:
    - mfas: Basic whole genome sequence in multifasta file of contigs. This
    can be used to process a finished genome in a single Fasta file as well.
    - cgbk: All contigs concatenated in a single GenBank file (Genoscope,
    French WGS). This can be used to process a finished genome in a single
    GanBank file as well.
    # TODO: provide support for other possible input formats

    Unpacking 'cgbk' genomes involves an initial step to detect occurrences
    of the sequence separator and collect the start and stop coordinates of
    each contig. Each pair of coordinates can then be used to extract the
    contig sequence and create a SeqRecord for that contig, which SeqIO
    normally does when it unpacks multifasta files.

    """
    # set up inputs
    infile = genome['file'] #TODO: make GUI input loader (upstream)
    inpath = dirs['ori_g_dir']+infile
    g_name = genome['name']
    print " ", g_name, "...",
    # prep output destinations
    mfas_dir = dirs['mfas_contigs_dir']
    fas_dir = dirs['fas_contigs_dir']+g_name+"/"
    ensure_dir([mfas_dir, fas_dir])
    mfas_file = mfas_dir+g_name+"_contigs.fas"
    records = []
    # select unpacking method
    if genome['input'] is 'mfas':
        try: path.exists(inpath) is True
        except ValueError: raise Exception("Bad input file path")
        genome_recs = load_multifasta(inpath)
        # generate GenBank files
        counter = 0
        for rec in genome_recs:
            counter +=1
            ctg_num = str(counter)
            new_id = g_name+"_"+str(counter)  # workaround for long ids
            new_seq = rec.seq
            new_seq.alphabet = generic_dna
            new_rec = SeqRecord(seq=new_seq, id=new_id)
            records.append(new_rec)  # for multifasta output
            fas_file = fas_dir+new_id+".fas"
            write_fasta(fas_file, new_rec)
    elif genome['input'] is 'cgbk':
        # load in genome data
        genome_rec = load_genbank(inpath)
        g_string = genome_rec.seq
        # find split coordinates
        coord_pairs = multisplit_finder(g_string, separator)
        # split record
        counter = 0
        for (start, stop) in coord_pairs:
            counter +=1
            ctg_num = str(counter)
            new_record = genome_rec[start:stop]
            new_record.id = g_name+"_"+ctg_num
            records.append(new_record)  # for multifasta output
            fas_file = fas_dir+g_name+"_"+ctg_num+".fas"
            write_fasta(fas_file, new_record)
    else:
        xmsg = "Input file format "+genome['input']+" unspecified/unsupported"
        raise Exception(xmsg)
    print counter, "contigs"
    # write master file
    write_fasta(mfas_file, records)

def extract_seg(contig, run_id):
    """Extract reference segments using coordinates."""
    # file info
    ctg_name = contig['name']
    run_root = p_root_dir+run_id+"/"
    in_file = dirs['ori_g_dir']+contig['file']
    out_root = run_root+dirs['ref_seg_dir']+ctg_name+"/"
    ensure_dir([out_root])
    print " ", ctg_name, "...",
    # open record
    record = load_genbank(in_file)
    count = 0
    for ref in contig['refs']:
        # extract segment
        segment = record[ref['coords'][0]:ref['coords'][1]]
        segment.id = ctg_name+"_"+ref['type']
        # write to file
        out_file = out_root+ctg_name+"_"+ref['type']+".fas"
        write_fasta(out_file, segment)
        count +=1
    print count, "segments"

def build_scaffolds(contig, run_id):
    """Build a scaffold of contigs based on the reference.

    This takes contigs that gave positive hits when blasted with reference
    segments. The contigs were aligned against the complete reference in a
    previous step for mapping purposes. Now the output of that step is re-used
    determine their position. A caveat is that if there are natural local
    rearrangements in the sequence relative to the reference, they may not be
    resolved appropriately. The problem is somewhat moderated by the fact that
    this function takes the best (usually the largest) hit region as "anchor"
    to position the contig within the scaffold. But if the rearranged region
    takes up a significant portion of the contig length, the anchoring will
    probably not be called correctly. Visual inspection of the finalized
    maps should help diagnose any such problems. The order can be fixed
    manually using the Mauve Contig Mover, which is part of Mauve 2.

    In some cases it is clear from looking at the maps generated further on
    that some contigs are included based on spurious hits. It is possible to
    make them be ignored and left out of the scaffold by listing their ID
    number in the genome dictionaries in the config file then rerunning the
    pipeline from this step.

    """
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    run_root = p_root_dir+run_id+"/"
    ctgs_root = run_root+dirs['match_out_dir']+ctg_name+"/"
    mauve_root = run_root+dirs['mauve_out_dir']+ctg_name+"/contigs/"
    scaffolds_root = run_root+dirs['scaffolds_dir']+ctg_name+"/"
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        ctgs_dir = ctgs_root+"/"+g_name+"/"
        print "\t", g_name, "...",
        # set outputs
        mauve_dir = mauve_root+g_name+"/"
        scaffolds_dir = scaffolds_root+g_name+"/"
        ensure_dir([mauve_dir, scaffolds_dir])
        scaff_fas = scaffolds_dir+g_name+"_"+ctg_name+"_scaffold.fas"
        scaff_gbk = scaffolds_dir+g_name+"_"+ctg_name+"_scaffold.gbk"
        # list genbank files in matches directory
        dir_contents = listdir(ctgs_dir)
        anchors_array = np.zeros(1, dtype=[('ctg', 'i4'),
                                           ('start', 'i4'),
                                           ('end', 'i4'),
                                           ('orient', 'i2')])
        for item in dir_contents:
            pattern = re.compile(r'.*_(\d*)\.gbk$')
            match = pattern.match(item)
            if match:
                ctg_num = match.group(1)
                if int(ctg_num) not in genome['ignore']:
                    print ctg_num,
                    # set inputs
                    mauve_file = mauve_dir+ctg_num+".mauve"
                    bb_file = mauve_file+".backbone"
                    try:
                        # parse Mauve output
                        coords = mauver_load2_k0(bb_file, prox_D)
                        # determine which segment to use as anchor
                        anchor_seg = get_anchor_loc(coords)
                        anchors_array = np.insert(anchors_array, 0,
                                                  (ctg_num,
                                                   anchor_seg['start'],
                                                   anchor_seg['end'],
                                                   anchor_seg['orient']))
                    except IOError:
                        print "ERROR: Mauve alignment file unavailable"
                        print "\t\t",
        # abort if there is no valid contig to proceed with
        try:
            assert len(anchors_array) > 1 # always 1 from stub
        except AssertionError:
            print "WARNING: Contig list empty"
        else:
            # order contigs by anchor location
            anchors_array = np.sort(anchors_array, order='start')
            # load contig records from the genbank files in the matches directory
            ctg_list = []
            for ctg_anchor in anchors_array:
                ctg_num = ctg_anchor['ctg']
                if ctg_num > 0:
                    contig_gbk = ctgs_dir+g_name+"_"+str(ctg_num)+".gbk"
                    record = load_genbank(contig_gbk)
                    if ctg_anchor['orient'] == -1: # flip record
                        record.seq = record.seq.reverse_complement()
                    ctg_list.append(record)
                else: # workaround for having 0 value leftover from stub
                    pass # having it might come in handy in later dev
            # output scaffold files
            write_fasta(scaff_fas, ctg_list)
            scaff_record = SeqRecord('', id='temp')
            scaff_bumper = SeqRecord(separator, id='join')
            for record in ctg_list:
                feat_start = len(scaff_record.seq)
                scaff_record += record
                feat_stop = len(scaff_record.seq)
                scaff_record += scaff_bumper
                feat_loc = FeatureLocation(feat_start, feat_stop)
                feature = SeqFeature(location=feat_loc,
                                     type='contig',
                                     qualifiers={'id': record.id})
                scaff_record.features.append(feature)
            scaff_record.id = g_name+"_"+ctg_name
            write_genbank(scaff_gbk, scaff_record[:-100]) # rm last bumper
            print ""
