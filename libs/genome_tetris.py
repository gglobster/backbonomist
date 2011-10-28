import re
import numpy as np
from os import path, listdir
from loaders import load_genbank, load_multifasta
from writers import write_genbank, write_fasta
from string_ops import multisplit_finder
from common import ensure_dir
from config import separator, directories as dirs, genomes
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import generic_dna
from aligning import align_mauve
from parsing import mauver_load_k0
from array_tetris import get_anchor_loc

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
    gbk_dir = dirs['gbk_contigs_dir']+g_name+"/"
    ensure_dir(gbk_dir)
    fas_file = mfas_dir+g_name+"_contigs.fas"
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
            new_id = g_name+"_"+str(counter)  # workaround for long ids
            new_seq = rec.seq
            new_seq.alphabet = generic_dna
            new_rec = SeqRecord(seq=new_seq, id=new_id)
            records.append(new_rec)  # for multifasta output
            gbk_file = gbk_dir+new_id+".gbk"
            write_genbank(gbk_file, new_rec)
        print counter, "records"
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
            new_record = genome_rec[start:stop]
            new_record.id = g_name+"_"+str(counter)
            records.append(new_record)  # for multifasta output
            gbk_file = gbk_dir+g_name+"_"+str(counter)+".gbk"
            write_genbank(gbk_file, new_record)
        print counter, "records"
    else:
        xmsg = "Input file format "\
               +genome['input']\
               +" unspecified or unsupported"
        raise Exception(xmsg)
    # write master file
    write_fasta(fas_file, records)
    return genome['name']

def extract_seg(contig):
    """Extract reference segments using coordinates."""
    # file info
    ctg_name = contig['name']
    in_file = dirs['ref_ctg_dir']+contig['file']
    out_root = dirs['ref_seg_dir']+ctg_name+"/"
    ensure_dir(out_root)
    print " ", ctg_name, "...",
    # open record
    record = load_genbank(in_file)
    for ref in contig['refs']:
        # extract segment
        segment = record[ref['coords'][0]:ref['coords'][1]]
        segment.id = ctg_name+"_"+ref['type']
        # write to file
        out_file = out_root+ctg_name+"_"+ref['type']+".fas"
        write_fasta(out_file, segment)
        print ref['type'],
    print ""

def build_scaffolds(contig):
    """Build a scaffold of contigs based on the reference.

    This takes contigs that gave positive hits when blasted with reference
    segments. Now the contigs are aligned against the complete reference to
    determine their position. A caveat is that if there are natural local
    rearrangements in the sequence relative to the reference, they may not be
    resolved appropriately. The problem is somewhat moderated by the fact that
    this function takes the best (usually the largest) hit region as "anchor"
    to position the contig within the scaffold. But if the rearranged region
    takes up a significant portion of the contig length, the anchoring will
    probably not be called correctly. Visual inspection of the finalized
    maps should reveal any such problems. The order can be fixed manually
    using the Mauve Contig Mover, which is part of Mauve 2.

    """
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    ref_ctg_file = dirs['ref_ctg_dir']+contig['file']
    ctgs_root = dirs['match_out_dir']+ctg_name+"/"
    mauve_root = dirs['mauve_out_dir']+ctg_name+"/"
    scaffolds_root = dirs['scaffolds_dir']+ctg_name+"/"
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        ctgs_dir = ctgs_root+"/"+g_name+"/"
        print "\t", g_name,
        # set outputs
        mauve_dir = mauve_root+g_name+"/"
        scaffolds_dir = scaffolds_root+g_name+"/"
        ensure_dir(mauve_dir)
        ensure_dir(scaffolds_dir)
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
                print ctg_num,
                # set inputs
                file_list = (ref_ctg_file, ctgs_dir+item)
                # set Mauve output
                mauve_outfile = mauve_dir+ctg_num+".mauve"
                # do Mauve alignment
                align_mauve(file_list, mauve_outfile)
                # parse Mauve output
                coords = mauver_load_k0(mauve_outfile+".backbone", 2)
                # determine which segment to use as anchor
                anchor_seg = get_anchor_loc(coords)
                anchors_array = np.insert(anchors_array, 0,
                                          (ctg_num,
                                           anchor_seg['start'],
                                           anchor_seg['end'],
                                           anchor_seg['orient']))
        # order contigs by anchor location
        anchors_array = np.sort(anchors_array, order='start')
        # load contig records from the genbank files in the matches directory
        ctg_list = []
        for ctg_anchor in anchors_array:
            ctg_num = ctg_anchor['ctg']
            if ctg_num > 0:
                contig_gbk = ctgs_dir+g_name+"_"+str(ctg_num)+".gbk"
                record = load_genbank(contig_gbk)
                ctg_list.append(record)
            else: # workaround for having a 0 value leftover from stub
                pass # having it might come in handy in later dev, so it stays
        # output scaffold files (for genbank, need to concatenate records)
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
        write_genbank(scaff_gbk, scaff_record[:-100]) # leave out last bumper
        print ""
