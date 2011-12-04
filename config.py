import numpy

# TODO: setup script to generate config file or through GUI

# Input data

project_id = 'BCSL'
project_date = '2011'
prot_db_name = 'Bacteria_prot'

from sets.broad import pXO1_positives as genomes
from sets.references import pXO1 as references

# Proximity thresholds for clumping
prox_D = 2000   # for ballpark estimation
prox_F = 100    # for fine alignment (set to 0 to skip clumping)

# Chopping size to limit length of detailed alignments
max_size = 3000
chop_mode = 'exact_size' # 'maxsize_bisect','maxsize_divisor','count_divisor'

# Identity percentage cutoffs and color coding
idpt = {95: '#444444',     # top similarity class (HexColor('#444444'))
        85: '#777777',     # upper middle class (HexColor('#777777'))
        70: '#BBBBBB',     # lower middle class (HexColor('#BBBBBB'))
        50: '#DDDDDD',     # low similarity class (HexColor('#DDDDDD'))
         0: '#FFFFFF'}     # lower than cutoff (HexColor('#FFFFFF'))

# Thresholds for binning contig sizes
ctg_thresholds = [100, 500, 1000]

# Function categories and legend (keys MUST be lowercase)

fct_flags = {'mge': ('transposase', 'transposon', 'intron',
                     'tyrosine recombinase', 'dna-invertase',
                     'reverse transcriptase'),
             'rep': ('something else', 'replication'),
             'syn': ('synthase', 'something else'),
             'tox': ('protective antigen', 'lethal factor',
                     'virulence factor'),
             'ger': ('spore germination protein', 'something else'),
             'tra': ('type iv', 'type ii/iv', 'topoisomerase',
                     'dna translocase ftsk', 'conjugal transfer',
                     'conjugation protein'),
             'ctl': ('transcriptional regulator', 'regulatory protein',
                     'regulator', 'transcriptional repressor'
                     'response regulator aspartate phosphatase'),
             'unk': ('uncharacterized', 'conserved domain protein'),
             'def': ('no match', 'hypothetical protein', 'pXO1-')} # default

fct_colors = {'mge': ('#66CC00', 'MGE'),
              'rep': ('#FF9900', 'Replication'),
              'syn': ('#FF00CC', 'Synthesis'),
              'tox': ('#FF0000', 'Pathogenesis'),
              'ger': ('#993333', 'Germination'),
              'tra': ('#6666FF', 'Transfer'),
              'ctl': ('#FFCC00', 'Control'),
              'unk': ('#666666', 'Uncharacterized'),
              'oth': ('#CCCCCC', 'Other'),
              'def': ('#FFFFFF', 'No match')}

# Concat contigs separator
separator = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"

# Project root directory
g_root_dir = 'data/'+project_id+'genomes/'
p_root_dir = g_root_dir+'runs/'

# run-independent directories
fixed_dirs = {
'ref_dbs_dir': 'data/ref_dbs/',
'ori_g_dir': g_root_dir+'original/',
'gbk_contigs_dir': g_root_dir+'genbank_ctgs/',
'fas_contigs_dir': g_root_dir+'fasta_ctgs/',
'mfas_contigs_dir': g_root_dir+'mfasta_ctgs/',
'blast_db_dir': g_root_dir+'blast_db/',
'annot_trn_dir': g_root_dir+'annotation/training/',
'ctg_cds_dir': g_root_dir+'annotation/genes/',
'ctg_prot_dir': g_root_dir+'annotation/proteins/',
'ctg_blast_dir': g_root_dir+'annotation/blastp/',
'ctg_stats': g_root_dir+'contig_stats/'
}

# run-dependent directories
run_dirs = {
'ref_pickles': 'references/',
'ref_seg_dir': 'references/segments/',
'ref_gbk_dir': 'references/genbank/',
'ref_fas_dir': 'references/fasta/',
'ref_map_dir': 'maps/references/',
'match_pickles': 'matching/',
'blast_out_dir': 'matching/blastn/',
'match_out_dir': 'matching/matches/',
'scaffolds_dir': 'scaffolds/',
'run_gbk_ctgs_dir': 'gbk_contigs/',
'mauve_out_dir': 'alignments/mauve_out/',
'aln_seg_dir': 'alignments/aln_segments/',
'maps_dir': 'maps/',
'reports': 'reports/'
}

# Blast parameters
blast_prefs = {'evalue': 0.01,
               'outfmt_pref': 6}
min_match = 500     # min size for a blast hit to be considered relevant
min_score = 1000    # min score for a blast hit to be considered relevant

# Blast results arrays datatypes
blast_dtypes = numpy.dtype([('query', 'S16'),
                           ('dbhit', 'S32'),
                           ('idp', 'float'),
                           ('mlen', 'uint8'),
                           ('mms', 'uint8'),
                           ('gaps', 'uint8'),
                           ('q_start', 'uint16'),
                           ('q_end', 'uint16'),
                           ('r_start', 'uint16'),
                           ('r_end', 'uint16'),
                           ('evalue', 'S5'),
                           ('bitscore', 'float')])

# Mauve executable
mauve_exec = 'progressiveMauve'

# Datatypes
mtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4')]
segtype = [('A','i4'), ('B','i4'), ('C', 'i4'), ('D','i4'), ('E','i4')]
