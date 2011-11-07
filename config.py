import numpy

# TODO: setup script to generate config file

# Input data

project_id = 'BCSL'
prot_db_name = 'Bacteria_prot'

genomes = [{'name': 'pXO1', 'input': 'cgbk', 'file': 'pXO1.gbk',
            'offset': (5000,0), 'ignore': (0, 0)},
#           {'name': 'pBc239', 'input': 'mfas', 'file': 'CP000228.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'pAH187_270', 'input': 'mfas', 'file': 'CP001179.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'pAH820_272', 'input': 'mfas', 'file': 'CP001285.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'p03B_179', 'input': 'mfas', 'file': 'CP001406.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'NZ_ACNE0', 'input': 'mfas', 'file': 'NZ_ACNE0.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'NZ_ACNF0', 'input': 'mfas', 'file': 'NZ_ACNF0.fas',
#            'offset': (0,0), 'ignore': (151, 0)},
#           {'name': 'NZ_ACNI0', 'input': 'mfas', 'file': 'NZ_ACNI0.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'NZ_ACNJ0', 'input': 'mfas', 'file': 'NZ_ACNJ0.fas',
#            'offset': (0,0), 'ignore': (204, 0)},
#           {'name': 'NZ_ACNK0', 'input': 'mfas', 'file': 'NZ_ACNK0.fas',
#            'offset': (15000,0), 'ignore': (0, 149)}, # 68?
#            {'name': 'VD_022a', 'input': 'mfas',
#             'file': 'G8177_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'VD_022c', 'input': 'mfas',
#             'file': 'G8177_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X1-1a', 'input': 'mfas',
#             'file': 'G13150_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X1-1c', 'input': 'mfas',
#             'file': 'G13150_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X1-2a', 'input': 'mfas',
#             'file': 'G13151_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X1-2c', 'input': 'mfas',
#             'file': 'G13151_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X2-1a', 'input': 'mfas',
#             'file': 'G13156_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X2-1c', 'input': 'mfas',
#             'file': 'G13156_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1O-2a', 'input': 'mfas',
#             'file': 'G13159_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1O-2c', 'input': 'mfas',
#             'file': 'G13159_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1O-3a', 'input': 'mfas',
#             'file': 'G13160_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1O-3c', 'input': 'mfas',
#             'file': 'G13160_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X1-3a', 'input': 'mfas',
#             'file': 'G13172_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG1X1-3c', 'input': 'mfas',
#             'file': 'G13172_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'AND1407a', 'input': 'mfas',
#             'file': 'G13175_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'AND1407c', 'input': 'mfas',
#             'file': 'G13175_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG2X1-3a', 'input': 'mfas',
#             'file': 'G13211_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
            {'name': 'BAG2X1-3c', 'input': 'mfas',
             'file': 'G13211_contigs.fasta',
             'offset': (0,0), 'ignore': (0, 0)}]

            # p03B_179 = p03BB102_179

references = [{'name': 'pXO1', 'file': 'pXO1.gbk',
               'refs': ({'coords': (1127,6561),'type': '3_10'},
                        {'coords': (8991,13756),'type': '12_17'},
                        {'coords': (19200,20894),'type': '20_20'},
                        {'coords': (21714,25117),'type': '23_25'},
                        {'coords': (28490,30979),'type': '30_34'},
                        {'coords': (38284,40203),'type': '45_46'},
                        {'coords': (52063,55796),'type': '64_65'},
                        {'coords': (55987,60586),'type': '66_69'},
                        {'coords': (61404,63518),'type': '71_74'},
                        {'coords': (64090,68866),'type': '75_79'},
                        {'coords': (69021,74102),'type': '80_84'},
                        {'coords': (74119,77871),'type': '85_90'},
                        {'coords': (79735,82867),'type': '93_97'},
                        {'coords': (83606,86696),'type': '98_101'},
                        {'coords': (87297,91464),'type': '103_107'},
                        {'coords': (91495,95163),'type': '108_108'},
                        {'coords': (95212,98847),'type': '109_114'},
                        {'coords': (100828,102495),'type': '117_118'},
                        {'coords': (155257,157413),'type': '182_184'},
                        {'coords': (174223,176584),'type': '207_212'})}]
                        # TODO: make coords loader

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
p_root_dir = 'data/'+project_id+'/'

directories = {
'ref_dbs_dir': 'data/ref_dbs/',
'ori_g_dir': p_root_dir+'genomes/original/',
'gbk_contigs_dir': p_root_dir+'genomes/genbank_ctgs/',
'mfas_contigs_dir': p_root_dir+'genomes/mfasta_ctgs/',
'blast_db_dir': p_root_dir+'genomes/blast_db/',
'annot_trn_dir': p_root_dir+'genomes/annot_trn/',
'annot_ctg_dir': p_root_dir+'genomes/annot_ctgs/',
# the following are run-dependent
'ref_seg_dir': 'ref_segments/',
'blast_out_dir': 'blast_out/',
'match_out_dir': 'matches/',
'scaffolds_dir': 'scaffolds/',
'mauve_out_dir': 'mauve_out/',
'scaff_annot_dir': 'annotation/',
'constructs_dir': 'constructs/',
'aln_seg_dir': 'aln_segments/',
'maps_dir': 'maps/'
}

# Blast parameters
blast_prefs = {'evalue': 0.01,
               'outfmt_pref': 6}
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

# Proximity thresholds for clumping
prox_D = 2000   # for ballpark estimation
prox_F = 100    # for fine alignment (set to 0 to skip clumping)

# Chopping size to limit length of detailed alignments
max_size = 3000
chop_mode = 'exact_size' # 'maxsize_bisect','maxsize_divisor','count_divisor'

# Identity percentage cutoffs and color coding
idpt = {95: '#444444',     # top similarity class (HexColor('#444444'))
        80: '#777777',     # upper middle class (HexColor('#777777'))
        70: '#BBBBBB',     # lower middle class (HexColor('#BBBBBB'))
        50: '#DDDDDD',     # low similarity class (HexColor('#DDDDDD'))
         0: '#FFFFFF'}     # lower than cutoff (HexColor('#FFFFFF'))

