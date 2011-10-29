from sys import argv, exit
from libs.common import ensure_dir
from libs.genome_tetris import unpack_genomes, extract_seg, build_scaffolds
from libs.blasting import make_genome_DB, basic_batch_blastn
from libs.parsing import glompX_blast_out
from libs.annotation import annot_scaffolds
from libs.reporting import map_scaffolds, map_pairwise
from libs.aligning import align2ref
from config import directories as dirs, backbones, genomes

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python main_script.py [step#]\n", \
          "For better results, run from within iPython."
    exit()

print "\n", \
      "##################################################\n", \
      "### Backbonomist v. 0.2                        ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) < 2:
    step = 0
else:
    step = int(argv[1])

if step is 0:
    print "\n###", step, ". Set up the work environment ###\n"
    for dir_name in dirs.keys():
        ensure_dir(dirs[dir_name])
    step +=1

if step is 1:
    print "\n###", step, ". Unpack genomes (separate contigs)###\n"
    for genome in genomes:
        unpack_genomes(genome)
    step +=1

if step is 2:
    print "\n###", step, ". Make Blast databases ###\n"
    for genome in genomes:
        make_genome_DB(genome)
    step +=1

if step is 3:
    print "\n###", step, ". Extract reference segments ###\n"
    for bb in backbones:
        extract_seg(bb)
    step +=1

if step is 4:
    print "\n###", step, ". Blast reference segments against genome DBs ###\n"
    for bb in backbones:
        basic_batch_blastn(bb)
    step +=1

if step is 5:
    print "\n###", step, ". Collect Blast results ###\n"
    for bb in backbones:
        glompX_blast_out(bb)
    step +=1

if step is 6:
    print "\n###", step, ". Build backbone scaffolds ###\n"
    for bb in backbones:
        build_scaffolds(bb)
    step +=1

if step is 7:
    print "\n###", step, ". Annotate backbone scaffolds ###\n"
    for bb in backbones:
        annot_scaffolds(bb)
    step +=1

if step is 8:
    print "\n###", step, ". Generate annotated scaffold maps ###\n"
    for bb in backbones:
        map_scaffolds(bb)
    step +=1

if step is 9:
    print "\n###", step, ". Align pairwise to reference ###\n"
    for bb in backbones:
        align2ref(bb)
    step +=1

if step is 10:
    print "\n###", step, ". Generate pairwise alignment maps ###\n"
    for bb in backbones:
        map_pairwise(bb)
    step +=1
    
if step > 10:
    print "\n### Nothing more to do! ###\n"