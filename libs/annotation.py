from config import genomes, directories as dirs, blast_prefs, \
    prot_db_name, project_id, p_root_dir
import re, subprocess
from loaders import load_multifasta, load_fasta
from writers import write_genbank
from common import ensure_dir
from blasting import local_blastp_2file
from parsing import collect_cogs
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
from os import path, listdir
from shutil import copyfile

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

def annot_contigs(ref_contig, run_id):
    """Annotate contigs (predict ORFs and assign function)."""
    # locate the COG database
    prot_db = dirs['ref_dbs_dir']+prot_db_name
    # TODO: add other DB / pfams?
    # set inputs and outputs
    ref_name = ref_contig['name'] # reference contig
    run_root = p_root_dir+run_id+"/"
    fas_ctgs_root = run_root+dirs['match_out_dir']+ref_name+"/"
    ctg_cds_root = dirs['ctg_cds_dir']
    ctg_prot_root = dirs['ctg_prot_dir']
    ctg_blast_root = dirs['ctg_blast_dir']
    g_gbk_ctgs_root = dirs['gbk_contigs_dir']
    r_gbk_ctgs_root = run_root+dirs['run_gbk_ctgs_dir']+ref_name+"/"
    annot_trn_root = dirs['annot_trn_dir']
    print " ", ref_name
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        fas_ctgs_dir = fas_ctgs_root+g_name+"/"
        g_file = dirs['ori_g_dir']+genome['file']
        print '\t', g_name, "...",
        # set output files
        training_file = annot_trn_root+g_name+"_annot.trn"
        # set output dirs
        ctg_cds_dir = ctg_cds_root+g_name+"/"
        ctg_prot_dir = ctg_prot_root+g_name+"/"
        ctg_blast_dir = ctg_blast_root+g_name+"/"
        g_gbk_ctgs_dir = g_gbk_ctgs_root+g_name+"/"
        r_gbk_ctgs_dir = r_gbk_ctgs_root+g_name+"/"
        annot_trn_dir = annot_trn_root+g_name+"/"
        ensure_dir([ctg_cds_dir, ctg_prot_dir, ctg_blast_dir,
                    g_gbk_ctgs_dir, r_gbk_ctgs_dir, annot_trn_dir])
        # list fasta files in matches directory
        dir_contents = listdir(fas_ctgs_dir)
        for item in dir_contents:
            pattern = re.compile(r'.*_(\d*)\.fas$')
            match = pattern.match(item)
            if match:
                ctg_num = match.group(1)
                print ctg_num,
                # set inputs and outputs
                ctg_fas = fas_ctgs_dir+item
                g_ctg_gbk = g_gbk_ctgs_dir+g_name+"_"+ctg_num+".gbk"
                r_ctg_gbk = r_gbk_ctgs_dir+g_name+"_"+ctg_num+".gbk"
                annot_gbk = ctg_cds_root+g_name+"_"+ctg_num+"_cds.gbk"
                annot_aa = ctg_prot_root+g_name+"_"+ctg_num+"_aa.fas"
                blast_out = ctg_blast_dir+g_name+"_"+ctg_num+".xml"
                if not path.exists(r_ctg_gbk):
                    if not path.exists(g_ctg_gbk):
                        # gene prediction
                        if not path.exists(training_file):
                            train_prodigal(g_file, training_file, "-q")
                        if not path.exists(annot_aa):
                            run_prodigal(ctg_fas, annot_gbk, annot_aa, training_file,
                                         "-q")
                        # blast the amino acids against COG
                        if not path.exists(blast_out):
                            local_blastp_2file(annot_aa, prot_db, blast_out,
                                               blast_prefs)
                        # collect best hits
                        rec_cogs = collect_cogs(blast_out)
                         # consolidate annotated genbank file
                        record = load_fasta(ctg_fas)
                        record.features = []
                        aa_record = load_multifasta(annot_aa)
                        counter = 1
                        for aa_rec in aa_record:
                            this_prot = 'Query_'+str(counter)
                            annotation = rec_cogs[this_prot]
                            # get feature details from description line
                            # because prodigal output fails to load as valid genbank
                            defline = aa_rec.description
                            pattern = re.compile('.+#\s(\d+)\s#\s(\d+)\s#\s(\S*1)\s#\sID.+')
                            match = pattern.match(defline)
                            start_pos = int(match.group(1))
                            end_pos = int(match.group(2))
                            strand_pos = int(match.group(3))
                            feat_loc = FeatureLocation(start_pos, end_pos)
                            l_tag = g_name+"_"+ctg_num+"_"+str(counter)
                            # consolidation feature annotations
                            quals = {'note': defline, 'locus_tag': l_tag,
                                     'fct': annotation, 'translation': aa_rec.seq}
                            feature = SeqFeature(location=feat_loc,
                                                 strand=strand_pos,
                                                 id='cds_'+str(counter),
                                                 type='CDS',
                                                 qualifiers=quals)
                            record.features.append(feature)
                            counter +=1
                        record.description = g_name+"_"+ctg_num
                        record.name = g_name+"_"+ctg_num
                        record.dbxrefs = ["Project: "+project_id+"/"+ref_name
                                          +"-like backbones"]
                        record.seq.alphabet = generic_dna
                        write_genbank(g_ctg_gbk, record)
                        copyfile(g_ctg_gbk, r_ctg_gbk)
                    else:
                        copyfile(g_ctg_gbk, r_ctg_gbk)
        print ""