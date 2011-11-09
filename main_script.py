from sys import argv, exit
import time
from datetime import datetime
from libs.genome_tetris import unpack_genomes, extract_seg, build_scaffolds
from libs.blasting import make_genome_DB, basic_batch_blastn
from libs.parsing import glompX_blast_out
from libs.annotation import annot_contigs
from libs.mapping import prep_maps
from libs.aligning import align_cstrct2ref, align_ctg2ref
from config import references, genomes

print "\n", \
      "##################################################\n", \
      "### Backbonomist v. 0.2                        ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python main_script.py [run_id] [step#]\n", \
          "Note that these arguments are positional: order matters!\n"
    exit()

if len(argv) > 1:
    run_id = argv[1]
else:
    run_id = str(int(time.time())) # use timestamp as unique run identifier

if len(argv) > 2:
    step = int(argv[2])
else:
    step = 0

start_timestamp = str(datetime.now())

if step is 0:
    print "\n###", step, ". Set up logging & reporting ###\n"
    print "to do"
#    for dataset in datasets:
#        save_parameters(dataset, max_pairs, run_id, start_timestamp)
#        log_start_run(dataset, run_id, start_timestamp)
    step +=1

if step is 1:
    print "\n###", step, ". Unpack genomes (separate contigs) ###\n"
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
    for ref in references:
        extract_seg(ref, run_id)
    step +=1

if step is 4:
    print "\n###", step, ". Blast reference segments against genomes ###\n"
    for ref in references:
        basic_batch_blastn(ref, run_id)
    step +=1

if step is 5:
    print "\n###", step, ". Collect Blast results ###\n"
    for ref in references:
        glompX_blast_out(ref, run_id)
    step +=1

if step is 6:
    print "\n###", step, ". Annotate matching contigs ###\n"
    for ref in references:
        annot_contigs(ref, run_id)
    step +=1

if step is 7:
    print "\n###", step, ". Align contigs pairwise to reference ###\n"
    for ref in references:
        align_ctg2ref(ref, run_id)
    step +=1

if step is 8:
    print "\n###", step, ". Construct backbone-based scaffolds ###\n"
    for ref in references:
        build_scaffolds(ref, run_id)
    step +=1

if step is 9:
    print "\n###", step, ". Align constructs pairwise to reference ###\n"
    for ref in references:
        align_cstrct2ref(ref, run_id)
    step +=1

if step is 10:
    print "\n###", step, ". Generate maps ###\n"
    for ref in references:
        prep_maps(ref, run_id)
    step +=1

if step > 10:
    print "\n### Nothing more to do! ###\n"

