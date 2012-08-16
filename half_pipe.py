from sys import exit
import time
import cPickle as pickle
import os
from datetime import datetime
from libs.common import ensure_dir
from libs.genome_tetris import process_ref, unpack_genomes, add_refs_2g
from libs.blasting import make_genome_DB, basic_batch_blast
from libs.parsing import glompX_blast_out
from libs.reporting import save_datasumm, log_start_run, log_end_run, \
    init_reports, log_resume_run, matches_table

from run_config import *
from run_sets import references, genomes

print "\n", \
      "##################################################\n", \
      "### Backbonomist v. 0.3 HALF PIPE              ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python half_pipe.py [run_id] [blast_mode] [step#] [g_select]\n", \
          "Note that these arguments are positional: order matters!\n"
    exit()

if len(argv) > 1:
    run_id = argv[1]
else:
    run_id = str(int(time.time())) # use timestamp as unique run identifier

if len(argv) > 2:
    blast_mode = argv[2]
else:
    blast_mode = 'n' # nucleotide blast by default

if len(argv) > 3:
    step = int(argv[3])
else:
    step = 0

if len(argv) > 4:
    g_select = argv[4:]
else:
    g_select = None

start_timestamp = str(datetime.now())

# ensure existence of all directories
ensure_dir(fixed_dirs.values())
run_dirs_go = ["".join([r_root_dir, run_id, "/", rdir])
               for rdir in run_dirs.values()]
ensure_dir(run_dirs_go)

# check for pickles
pickle_root = r_root_dir+run_id+"/"+run_dirs['pickles']+run_id
ref_pickles = pickle_root+"_refs.p"
genome_pickles = pickle_root+"_genomes.p"
match_pickles = pickle_root+"_matches.p"
# references
try: run_refs = pickle.load(open(ref_pickles, 'rb'))
except IOError:
    step = 0
    run_refs = []
# genomes
try: run_gs = pickle.load(open(genome_pickles, 'rb'))
except IOError:
    step = 0
    run_gs = []
# matches
try: run_matches = pickle.load(open(match_pickles, 'rb'))
except IOError:
    run_matches = []

## pipeline

if step is not 0:
    log_resume_run(run_id, base_root, project_id, start_timestamp, step)

elif step is 0:
    print "\n###", step, ". Set up logging & reporting ###\n"
    log_start_run(run_id, base_root, project_id, run_dirs, start_timestamp)
    save_datasumm(run_id, blast_mode, r_root_dir, run_dirs, genomes,
                  references, project_id, project_date, start_timestamp)
    init_reports(run_id, fixed_dirs, ctg_thresholds, start_timestamp)
    step +=1

if step is 1:
    print "\n###", step, ". Prepare references ###\n"
    for ref in references:
        timestamp = str(datetime.now())
        ref_obj = process_ref(ref, ref_annot_flag, r_root_dir, fixed_dirs,
                              run_dirs, run_id, timestamp)
        run_refs.append(ref_obj)
    if os.path.exists(ref_pickles):
        os.remove(ref_pickles)
    pickle.dump(run_refs, open(ref_pickles, 'wb'))
    step +=1

if step is 2:
    print "\n###", step, ". Prepare genomes ###\n"
    for genome in genomes:
        unpack_genomes(genome, separator, fixed_dirs)
        make_genome_DB(genome)
    run_gs = add_refs_2g(genomes, references)
    if os.path.exists(genome_pickles):
        os.remove(genome_pickles)
    pickle.dump(run_gs, open(genome_pickles, 'wb'))
    step +=1

if step is 3:
    print "\n###", step, ". Blast reference segments against genomes ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        basic_batch_blast(run_gs, ref, blast_mode, run_id, timestamp)
    step +=1

if step is 4:
    print "\n###", step, ". Collect Blast results ###\n"
    for ref in run_refs:
        timestamp = str(datetime.now())
        ref_hits, ctl_scores = glompX_blast_out(run_gs, ref, blast_mode,
                                                run_id, timestamp)
        ref_matches = {'ref': ref, 'run': run_id, 'hits': ref_hits,
                       'ctl': ctl_scores}
        run_matches.append(ref_matches)
    if os.path.exists(match_pickles):
        os.remove(match_pickles)
    pickle.dump(run_matches, open(match_pickles, 'wb'))
    step +=1

if step is 5:
    print "\n###", step, ". Make match results table & graphs ###\n"
    for match_dict in run_matches:
        timestamp = str(datetime.now())
        matches_table(match_dict, r_root_dir, run_dirs, timestamp)
    step +=1

if step > 5:
    stop_timestamp = str(datetime.now())
    log_end_run(run_id, base_root, project_id, stop_timestamp)
    print "\n### Nothing more to do! ###\n"

