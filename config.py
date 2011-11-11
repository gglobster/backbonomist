import numpy

# TODO: setup script to generate config file

# Input data

project_id = 'BCSL'
project_date = '2011'
prot_db_name = 'Bacteria_prot'

genomes = [{'name': 'pXO1', 'input': 'cgbk', 'file': 'pXO1.gbk',
            'offset': (0,0), 'ignore': (0, 0)},
           # known plasmids
#           {'name': 'pBCXO1', 'input': 'mfas', 'file': 'pBCXO1.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'pBc239', 'input': 'mfas', 'file': 'CP000228.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'pBc10987', 'input': 'mfas', 'file': 'AE017195.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'pAH187_270', 'input': 'mfas', 'file': 'CP001179.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'pAH820_272', 'input': 'mfas', 'file': 'CP001285.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'p03B_179', 'input': 'mfas', 'file': 'CP001406.fas',
#            'offset': (0,0), 'ignore': (0, 0)}, #p03BB102_179
#           {'name': 'pH30_258', 'input': 'mfas', 'file': 'CP001166.fas',
#            'offset': (0,0), 'ignore': (0, 0)}, #pH308197_258
#           # BI batch 1
#           {'name': 'IS075', 'input': 'mfas', 'file': 'IS075.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'Schrouff', 'input': 'mfas', 'file': 'Schrouff.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'TIAC129', 'input': 'mfas', 'file': 'TIAC129.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'VD022', 'input': 'mfas', 'file': 'VD022.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'VD142', 'input': 'mfas', 'file': 'VD142.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           {'name': 'VD022_454', 'input': 'mfas', 'file': 'VD022_454.fas',
#            'offset': (0,0), 'ignore': (0, 0)},
#           # BI batch 2
            {'name': 'VD_022a', 'input': 'mfas',
             'file': 'G8177_assembly.fasta',
             'offset': (0,0), 'ignore': (0, 0)},
            {'name': 'VD_022c', 'input': 'mfas',
             'file': 'G8177_contigs.fasta',
             'offset': (0,0), 'ignore': (0, 0)},
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
#             'offset': (0,61380), 'ignore': (0, 0)},
#            {'name': 'AND1407c', 'input': 'mfas',
#             'file': 'G13175_contigs.fasta',
#             'offset': (0,180750), 'ignore': (0, 0)},
#            {'name': 'BAG2X1-3a', 'input': 'mfas',
#             'file': 'G13211_assembly.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'BAG2X1-3c', 'input': 'mfas',
#             'file': 'G13211_contigs.fasta',
#             'offset': (0,0), 'ignore': (0, 0)},
#           # genbank genomes
#            {'name': 'NZ_AAEK0', 'input': 'mfas', 'file': 'NZ_AAEK0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_AAJM0', 'input': 'mfas', 'file': 'NZ_AAJM0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ABDA0', 'input': 'mfas', 'file': 'NZ_ABDA0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ABDL0', 'input': 'mfas', 'file': 'NZ_ABDL0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ABDM0', 'input': 'mfas', 'file': 'NZ_ABDL0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLS0', 'input': 'mfas', 'file': 'NZ_ACLS0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLT0', 'input': 'mfas', 'file': 'NZ_ACLT0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLU0', 'input': 'mfas', 'file': 'NZ_ACLU0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLV0', 'input': 'mfas', 'file': 'NZ_ACLV0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLW0', 'input': 'mfas', 'file': 'NZ_ACLW0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLX0', 'input': 'mfas', 'file': 'NZ_ACLX0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLY0', 'input': 'mfas', 'file': 'NZ_ACLY0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACLZ0', 'input': 'mfas', 'file': 'NZ_ACLZ0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMA0', 'input': 'mfas', 'file': 'NZ_ACMA0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMB0', 'input': 'mfas', 'file': 'NZ_ACMB0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMC0', 'input': 'mfas', 'file': 'NZ_ACMC0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMD0', 'input': 'mfas', 'file': 'NZ_ACMD0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACME0', 'input': 'mfas', 'file': 'NZ_ACME0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMF0', 'input': 'mfas', 'file': 'NZ_ACMF0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMG0', 'input': 'mfas', 'file': 'NZ_ACMG0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMH0', 'input': 'mfas', 'file': 'NZ_ACMH0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMI0', 'input': 'mfas', 'file': 'NZ_ACMI0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMJ0', 'input': 'mfas', 'file': 'NZ_ACMJ0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMK0', 'input': 'mfas', 'file': 'NZ_ACMK0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACML0', 'input': 'mfas', 'file': 'NZ_ACML0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMM0', 'input': 'mfas', 'file': 'NZ_ACMM0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMN0', 'input': 'mfas', 'file': 'NZ_ACMN0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMO0', 'input': 'mfas', 'file': 'NZ_ACMO0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMP0', 'input': 'mfas', 'file': 'NZ_ACMP0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMQ0', 'input': 'mfas', 'file': 'NZ_ACMQ0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMR0', 'input': 'mfas', 'file': 'NZ_ACMR0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMS0', 'input': 'mfas', 'file': 'NZ_ACMS0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMT0', 'input': 'mfas', 'file': 'NZ_ACMT0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMU0', 'input': 'mfas', 'file': 'NZ_ACMU0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMV0', 'input': 'mfas', 'file': 'NZ_ACMV0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMW0', 'input': 'mfas', 'file': 'NZ_ACMW0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMX0', 'input': 'mfas', 'file': 'NZ_ACMX0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMY0', 'input': 'mfas', 'file': 'NZ_ACMY0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACMZ0', 'input': 'mfas', 'file': 'NZ_ACMZ0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNA0', 'input': 'mfas', 'file': 'NZ_ACNA0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNB0', 'input': 'mfas', 'file': 'NZ_ACNB0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNC0', 'input': 'mfas', 'file': 'NZ_ACNC0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACND0', 'input': 'mfas', 'file': 'NZ_ACND0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
            {'name': 'NZ_ACNE0', 'input': 'mfas', 'file': 'NZ_ACNE0.fas',
             'offset': (0,129500), 'ignore': (0, 0)}
#            {'name': 'NZ_ACNF0', 'input': 'mfas', 'file': 'NZ_ACNF0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNG0', 'input': 'mfas', 'file': 'NZ_ACNG0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNH0', 'input': 'mfas', 'file': 'NZ_ACNH0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNI0', 'input': 'mfas', 'file': 'NZ_ACNI0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNJ0', 'input': 'mfas', 'file': 'NZ_ACNJ0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNK0', 'input': 'mfas', 'file': 'NZ_ACNK0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ACNL0', 'input': 'mfas', 'file': 'NZ_ACNL0.fas',
#             'offset': (0,0), 'ignore': (0, 0)},
#            {'name': 'NZ_ADFM0', 'input': 'mfas', 'file': 'NZ_ADFM0.fas',
#             'offset': (0,0), 'ignore': (0, 0)}
]

references = [{'name': 'pXO1', 'file': 'pXO1.gbk',
               'refs': ({'coords': (1127,6561),'name': '3_10'},
                        {'coords': (8991,13756),'name': '12_17'},
                        {'coords': (19200,20894),'name': '20_20'},
                        {'coords': (21714,25117),'name': '23_25'},
                        {'coords': (28490,30979),'name': '30_34'},
                        {'coords': (38284,40203),'name': '45_46'},
                        {'coords': (52063,55796),'name': '64_65'},
                        {'coords': (55987,60586),'name': '66_69'},
                        {'coords': (61404,63518),'name': '71_74'},
                        {'coords': (64090,68866),'name': '75_79'},
                        {'coords': (69021,74102),'name': '80_84'},
                        {'coords': (74119,77871),'name': '85_90'},
                        {'coords': (79735,82867),'name': '93_97'},
                        {'coords': (83606,86696),'name': '98_101'},
                        {'coords': (87297,91464),'name': '103_107'},
                        {'coords': (91495,95163),'name': '108_108'},
                        {'coords': (95212,98847),'name': '109_114'},
                        {'coords': (100828,102495),'name': '117_118'},
                        {'coords': (155257,157413),'name': '182_184'},
                        {'coords': (174223,176584),'name': '207_212'})}]
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
'fas_contigs_dir': p_root_dir+'genomes/fasta_ctgs/',
'mfas_contigs_dir': p_root_dir+'genomes/mfasta_ctgs/',
'blast_db_dir': p_root_dir+'genomes/blast_db/',
'annot_trn_dir': p_root_dir+'genomes/annotation/training/',
'ctg_cds_dir': p_root_dir+'genomes/annotation/genes/',
'ctg_prot_dir': p_root_dir+'genomes/annotation/proteins/',
'ctg_blast_dir': p_root_dir+'genomes/annotation/blastp/',
'ctg_stats': p_root_dir+'genomes/contig_stats/',
# the following are run-dependent
'ref_seg_dir': 'ref_segments/',
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
        85: '#777777',     # upper middle class (HexColor('#777777'))
        70: '#BBBBBB',     # lower middle class (HexColor('#BBBBBB'))
        50: '#DDDDDD',     # low similarity class (HexColor('#DDDDDD'))
         0: '#FFFFFF'}     # lower than cutoff (HexColor('#FFFFFF'))

# Thresholds for binning contig sizes
ctg_thresholds = [100, 500, 1000]
