__author__ = 'GG'

from sys import argv, exit
import time
import cPickle as pickle
from datetime import datetime
from libs.common import ensure_dir
from config import run_dirs, fixed_dirs, p_root_dir

print "\n", \
      "##################################################\n", \
      "### Backbonomist v. 0.3                        ###\n", \
      "### Copyright 2011 Geraldine A. Van der Auwera ###\n", \
      "##################################################\n", \

if len(argv) > 1 and argv[1] == '-h':
    print "Basic usage: \n", \
          "$ python run_script.py [run_id] [step#]\n", \
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