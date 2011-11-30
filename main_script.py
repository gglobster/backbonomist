from sys import argv, exit
import time
import cPickle as pickle
import os
from datetime import datetime
from libs.common import ensure_dir
from libs.genome_tetris import process_ref, unpack_genomes, build_scaffolds
from libs.blasting import make_genome_DB, basic_batch_blastn
from libs.parsing import glompX_blast_out
from libs.annotation import annot_genome_contigs
from libs.mapping import prep_maps
from libs.aligning import align_cstrct2ref, align_ctg2ref
from libs.reporting import save_datasumm, log_start_run, log_end_run, \
    init_reports, log_resume_run, matches_table
from config import references, genomes, run_dirs, fixed_dirs, p_root_dir

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

if len(argv) > 3:
    g_select = argv[3:]
else:
    g_select = None

start_timestamp = str(datetime.now())

# ensure existence of all directories
ensure_dir(fixed_dirs.values())
run_dirs_go = ["".join([p_root_dir, run_id, "/", rdir])
               for rdir in run_dirs.values()]
ensure_dir(run_dirs_go)

# check for pickled references
ref_pickles = p_root_dir+run_id+"/"+run_dirs['ref_pickles']+run_id+"_refs.p"
try: run_refs = pickle.load(open(ref_pickles, 'rb'))
except IOError:
    step = 0
    run_refs = []

# check for pickled matches
match_pickles = p_root_dir+run_id+"/"+run_dirs['match_pickles']+run_id\
                +"_matches.p"
try: run_matches = pickle.load(open(match_pickles, 'rb'))
except IOError:
    run_matches = []

if step is not 0:
    log_resume_run(run_id, start_timestamp, step)

elif step is 0:
    print "\n###", step, ". Set up logging & reporting ###\n"
    log_start_run(run_id, start_timestamp)
    save_datasumm(run_id, start_timestamp)
    init_reports(run_id, start_timestamp)
    step +=1

if step is 1:
    print "\n###", step, ". Prepare references ###\n"
    for ref in references:
        timestamp = str(datetime.now())
        ref_obj = process_ref(ref, run_id, timestamp)
        run_refs.append(ref_obj)
    if os.path.exists(ref_pickles):
        os.remove(ref_pickles)
    pickle.dump(run_refs, open(ref_pickles, 'wb'))
    step +=1

if step is 2:
    print "\n###", step, ". Prepare genomes ###\n"
    for genome in genomes:
        unpack_genomes(genome)
        make_genome_DB(genome)
    step +=1

if step is 3:
    print "\n###", step, ". Blast reference segments against genomes ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        basic_batch_blastn(ref, run_id, timestamp)
    step +=1

if step is 4:
    print "\n###", step, ". Collect Blast results ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        ref_hits, ctl_scores = glompX_blast_out(ref, run_id, timestamp)
        ref_matches = {'ref': ref, 'run': run_id, 'hits': ref_hits,
                       'ctl': ctl_scores}
        run_matches.append(ref_matches)
    if os.path.exists(match_pickles):
        os.remove(match_pickles)
    pickle.dump(run_matches, open(match_pickles, 'wb'))
    step +=1 

if step is 5:
    print "\n###", step, ". Annotate matching contigs ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        annot_genome_contigs(ref, run_id, timestamp)
    step +=1

if step is 6:
    print "\n###", step, ". Align contigs pairwise to reference ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        align_ctg2ref(ref, run_id, timestamp)
    step +=1

if step is 7:
    print "\n###", step, ". Construct backbone-based scaffolds ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        build_scaffolds(ref, run_id, timestamp)
    step +=1

if step is 8:
    print "\n###", step, ". Align constructs pairwise to reference ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        align_cstrct2ref(ref, run_id, timestamp)
    step +=1

if step is 9:
    print "\n###", step, ". Make match results table & graphs ###\n"
    for match_dict in run_matches:
        timestamp = str(datetime.now())
        matches_table(match_dict, timestamp)
    step +=1

if step is 10:
    print "\n###", step, ". Generate maps ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        prep_maps(ref, run_id, timestamp, g_select)
    step +=1

if step > 10:
    stop_timestamp = str(datetime.now())
    log_end_run(run_id, stop_timestamp)
    print "\n### Nothing more to do! ###\n"

