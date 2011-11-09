from config import genomes, directories as dirs, p_root_dir, blast_prefs, \
    prot_db_name
import re, subprocess
from loaders import load_genbank, load_multifasta, load_fasta
from writers import write_genbank
from common import ensure_dir
from blasting import local_blastp_2file
from parsing import collect_cogs
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
from os import path

def train_prodigal(seq_file, training_file, mode):
    """Train Prodigal on the entire dataset."""
    cline = "prodigal "+mode+" -i "+seq_file+" -t "+training_file
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def run_prodigal(in_file, an_gbk, an_aa, trn_file, mode):
    """Annotate sequence records individually using Prodigal."""
    cline = "prodigal "+mode+" -i "+in_file\
                            +" -o "+an_gbk\
                            +" -a "+an_aa\
                            +" -t "+trn_file
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    return output

def annot_ctgs(genome, fas_file, ctg_num):
    """Predict ORFs on contigs (no functional annotation)."""
    # set inputs and outputs
    g_name = genome['name']
    g_file = dirs['ori_g_dir']+genome['file']
    ctg_annot_root = dirs['annot_ctg_dir']+g_name+"/"
    ctg_predict_dir = ctg_annot_root+"gbk_feat/"
    ctg_aa_dir = ctg_annot_root+"aa/"
    ctg_gbk_dir = dirs['gbk_contigs_dir']+g_name+"/"
    trn_dir = dirs['annot_trn_dir']
    training_file = dirs['annot_trn_dir']+g_name+"_annot.trn"
    predict_gbk = ctg_predict_dir+g_name+"_"+ctg_num+"_feats.gbk"
    predict_aa = ctg_aa_dir+g_name+"_"+ctg_num+"_aa.fas"
    full_gbk = ctg_gbk_dir+g_name+"_"+ctg_num+".gbk"
    ensure_dir([trn_dir, ctg_predict_dir, ctg_aa_dir, ctg_gbk_dir])
    # predictions
    if not path.exists(training_file):
        train_prodigal(g_file, training_file, "-q")
    if not path.exists(predict_aa):
        run_prodigal(fas_file, predict_gbk, predict_aa, training_file, "-q")
    # consolidate annotated genbank file
    record = load_fasta(fas_file)
    record.features = []
    aa_record = load_multifasta(predict_aa)
    counter = 1
    for aa_rec in aa_record:
        # get feature details from description line
        # necessary because prodigal output fails to load as gbk record
        defline = aa_rec.description
        pattern = re.compile('.+#\s(\d+)\s#\s(\d+)\s#\s(\S*1)\s#\sID.+')
        match = pattern.match(defline)
        start_pos = int(match.group(1))
        end_pos = int(match.group(2))
        strand_pos = int(match.group(3))
        feat_loc = FeatureLocation(start_pos, end_pos)
        l_tag = g_name+"_"+ctg_num+"_"+str(counter)
        # consolidation feature annotations
        quals = {'note': defline, 'locus_tag': l_tag, 'fct': 'no match',
                 'translation': aa_rec.seq}
        feature = SeqFeature(location=feat_loc,
                             strand=strand_pos,
                             id='cds_'+str(counter),
                             type='CDS',
                             qualifiers=quals)
        record.features.append(feature)
        counter +=1
    record.description = g_name+"_"+ctg_num
    record.name = g_name+"_"+ctg_num
    record.seq.alphabet = generic_dna
    write_genbank(full_gbk, record)

def annot_contigs(ref_contig, run_id):
    """Annotate contigs (predict ORFs and assign function)."""
    # locate the COG database
    prot_db = dirs['ref_dbs_dir']+prot_db_name
    # TODO: add other DB / pfams?
    # set inputs and outputs
    ref_name = ref_contig['name'] # reference contig
    run_root = p_root_dir+run_id+"/"
    fas_ctgs_root = run_root+dirs['match_out_dir']+ref_name+"/"
    ctg_cds_root = run_root+dirs['ctg_cds_dir']+ref_name+"/"
    ctg_prot_root = run_root+dirs['ctg_prot_dir']+ref_name+"/"
    gbk_ctgs_root = run_root+dirs['gbk_ctgs_dir']+ref_name+"/"
    annot_trn_root = dirs['annot_trn_dir']
    ensure_dir([ctg_cds_root, ctg_prot_root, gbk_ctgs_root, annot_trn_root])
    print " ", ref_name
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        fas_ctgs_dir = fas_ctgs_root+g_name+"/"
        # set output files
        training_file = annot_trn_root+g_name+"_annot.trn"
        g_file = dirs['ori_g_dir']+genome['file']

        # set output dirs




        ctg_gbk = scaff_root+g_name+"/"+g_name+"_"+ref_name+"_scaffold.gbk"
        print '\t', g_name, "...",


        blast_out_dir = ctg_annot_root+"blastp/"

        ensure_dir([gbk_out_dir, aa_out_dir, blast_out_dir])

        annot_gbk = gbk_out_dir+g_name+"_"+ref_name+"_annot.gbk"
        annot_aa = aa_out_dir+g_name+"_"+ref_name+"_aa.fas"
        blast_out = blast_out_dir+g_name+"_"+ref_name+".xml"
        fin_gbk_out = scaff_cons_root+g_name+"_"+ref_name+"_cstrct.gbk"
        # abort if there is no scaffold construct
        try: open(scaff_gbk, 'r')
        except IOError:
            print "WARNING: No scaffold construct to annotate"
        else:
            # gene prediction
            print "predict",
            if not path.exists(training_file):
                train_prodigal(g_file, training_file, "-q")
            if not path.exists(annot_aa):
                run_prodigal(scaff_gbk, annot_gbk, annot_aa, training_file,
                             "-q")
            # blast the amino acids against COG
            print "blastp",
            if not path.exists(blast_out):
                local_blastp_2file(annot_aa, prot_db, blast_out, blast_prefs)
            # collect best hits
            rec_cogs = collect_cogs(blast_out)
            # consolidate annotated genbank file
            record = load_genbank(scaff_gbk)
            ctg_feats = [feature for feature in record.features
                         if feature.type == 'contig'] # backup contig features
            record.features = [] # wipe previous features
            for contig in ctg_feats:
                record.features.append(contig) # restore contig features
            counter = 1
            gene_count = 0
            while counter < len(rec_cogs)/2:
                # TODO: improve this whole bit, e.g. add protein translation etc
                gene_count +=1
                # get feature details from description line
                # necessary because prodigal output fails to load as gbk record
                this_prot = 'Query_'+str(counter)
                annotation = rec_cogs[this_prot]
                defline = rec_cogs[this_prot+"_def"]
                pattern = re.compile('.+#\s(\d+)\s#\s(\d+)\s#\s(\S*1)\s#\sID.+')
                match = pattern.match(defline)
                start_pos = int(match.group(1))
                end_pos = int(match.group(2))
                strand_pos = int(match.group(3))
                feat_loc = FeatureLocation(start_pos, end_pos)
                # consolidation feature annotations
                quals = {'note': defline, 'fct': annotation}
                feature = SeqFeature(location=feat_loc,
                                     strand=strand_pos,
                                     id='cds_'+str(counter),
                                     type='CDS',
                                     qualifiers=quals)
                record.features.append(feature)
                counter +=1
            print "annot"
            record.description = g_name+"_"+ctg_name+"_scaffold"
            record.name = g_name+"_"+ctg_name
            record.dbxrefs = ["Project: "+ctg_name+"-like backbones"]
            record.seq.alphabet = generic_dna
            write_genbank(fin_gbk_out, record)
