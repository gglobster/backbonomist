from config import genomes, directories as dirs, blast_prefs, prot_db_name
import re, subprocess
from loaders import load_multifasta, load_genbank
from writers import write_genbank
from common import ensure_dir
from blasting import local_blastp_2file
from parsing import collect_cogs
from Bio.Alphabet import generic_dna
from Bio.SeqFeature import SeqFeature, FeatureLocation
from drawing import ContigDraw

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

def annot_scaffolds(contig):
    """Annotate scaffolds."""
    # locate the COG database
    prot_db = dirs['ref_dbs_dir']+prot_db_name
    # TODO: add other DB / pfams?
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    scaff_root = dirs['scaffolds_dir']+ctg_name+"/"
    scaff_annot_root = dirs['scaff_annot_dir']+ctg_name+"/"
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs
        g_name = genome['name']
        g_file = dirs['ori_g_dir']+genome['file']
        scaff_gbk = scaff_root+g_name+"/"+g_name+"_"+ctg_name+"_scaffold.gbk"
        print '\t', g_name, "...",
        # set output dirs
        gbk_out_dir = scaff_annot_root+"predict/"
        aa_out_dir = scaff_annot_root+"aa/"
        blast_out_dir = scaff_annot_root+"blastp/"
        solid_out_dir = scaff_annot_root+"genbank/"
        ensure_dir(gbk_out_dir)
        ensure_dir(aa_out_dir)
        ensure_dir(blast_out_dir)
        ensure_dir(solid_out_dir)
        # set output files
        training_file = dirs['annot_trn_dir']+g_name+"_annot.trn"
        annot_gbk = gbk_out_dir+g_name+"_"+ctg_name+"_annot.gbk"
        annot_aa = aa_out_dir+g_name+"_"+ctg_name+"_aa.fas"
        blast_out = blast_out_dir+g_name+"_"+ctg_name+".xml"
        fin_gbk_out = solid_out_dir+g_name+"_"+ctg_name+"_fin.gbk"
        # gene prediction
        print "predict",
        train_prodigal(g_file, training_file, "-q")
        run_prodigal(scaff_gbk, annot_gbk, annot_aa, training_file, "-q")
        # blast the amino acids against COG
        print "blastp",
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
            # necessary because the prodigal output is not parser-friendly
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
