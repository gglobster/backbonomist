from os import path
from config import directories as dirs, p_root_dir, project_id, \
    project_date, genomes, references, blast_prefs, max_size, ctg_thresholds
from common import ensure_dir
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def log_start_run(run_id, timestamp):
    """Record run initiated in the dataset log."""
    set_log = p_root_dir+"/"+project_id+"_log.html"
    param_file = run_id+"/"+dirs['reports']+run_id+"_dataset.txt"
    header = "<p><b>Log of processing runs for "+project_id+"</b></p><p><ul>"
    set_log_htm = ["<li><b>", run_id, "</b>&nbsp;- initiated ", timestamp,
                   " (<a href='", param_file, "'>dataset</a>)</li>"]
    linkline = "".join(set_log_htm)
    if not path.isfile(set_log): # first time running this dataset?
        open(set_log, 'w').write(header)
    open(set_log, 'a').write(linkline)

def log_end_run(run_id, timestamp):
    """Record run in the dataset log."""
    set_log = p_root_dir+"/"+project_id+"_log.html"
    run_report = run_id+"/"+run_id+"_report.html"
    set_log_htm = ["<li><b>", run_id, "</b>", "&nbsp;- completed ", timestamp,
                   " (<a href='", run_report, "'>report</a>)</li>"]
    linkline = "".join(set_log_htm)
    open(set_log, 'a').write(linkline)

def log_resume_run(run_id, timestamp, step):
    """Record run resumed in the dataset log."""
    set_log = p_root_dir+"/"+project_id+"_log.html"
    run_report = run_id+"/"+run_id+"_report.html"
    set_log_htm = ["<li><b>", run_id, "</b>", "&nbsp;- resumed ", timestamp,
                   " (from step "+str(step)+"</li>"]
    linkline = "".join(set_log_htm)
    open(set_log, 'a').write(linkline)

def save_datasumm(run_id, timestamp):
    """Save a summary of the dataset composition to file."""
    print " ", run_id,
    report_root = p_root_dir+run_id+"/"+dirs['reports']
    param_file = report_root+run_id+"_dataset.txt"
    ensure_dir([report_root])
    # genome data
    g_title = "\t".join(["## Genomes"])
    g_header = "\t".join(["## Name", "Offset", "Format", "File"])
    g_data = [g_title, g_header]
    for genome in genomes:
        g_data.append("\t".join([genome['name'], str(genome['offset'][1]),
                                 genome['input'], genome['file']]))
    g_block = "\n".join(g_data)
    # references data
    r_title = "\t".join(["## References"])
    r_header = "\t".join(["## Name", "Segments", "File"])
    r_data = [r_title, r_header]
    for ref in references:
        r_data.append("\t".join([ref['name'], str(len(ref['refs'])),
                                 ref['file']]))
        r_data.append("\t".join(["###", "Segments"]))
        r_data.append("\t".join(["###", "Name", "Start", "Stop"]))
        for seg in ref['refs']:
            r_data.append("\t".join(["", seg['name'],
                                     str(seg['coords'][0]),
                                     str(seg['coords'][1])]))
    r_block = "\n".join(r_data)
    # text block
    txt = ["# Project", project_id,
           "# Run ID", run_id,
           "# Date generated", project_date,
           "# Date processing initiated", timestamp,
           "# Dataset composition", g_block, r_block]
    # write to file
    open(param_file, 'w').write("\n".join(txt))
    print "dataset summary saved to file"

def init_reports(run_id, timestamp):
    """Record run initiated in the dataset log."""
    # set inputs and outputs
    ctg_stats_file = dirs['ctg_stats']+"contig_stats.txt"
    param_file = run_id+"/"+dirs['reports']+run_id+"_dataset.txt"
    ensure_dir([dirs['ctg_stats']])
    # initialize contig stats report
    kb1, kb2, kb3 = ctg_thresholds
    cs_header = "# Contig size distribution statistics\n"
    if not path.isfile(ctg_stats_file): # first time running?
        open(ctg_stats_file, 'w').write(cs_header)
    cs_txt = ["#\nRun ", run_id, " initiated ", timestamp,
              "\nGenome", "\t[<", str(kb1), "]",
              "\t[", str(kb1), "-", str(kb2), "]",
              "\t[", str(kb2), "-", str(kb3), "]",
              "\t[>", str(kb3), "]"]
    open(ctg_stats_file, 'a').write("".join(cs_txt))

def ctg_stats(g_name, records):
    """Report on size distribution of contigs in a genome."""
    i = 1000
    kb1, kb2, kb3 = ctg_thresholds
    # separate length categories
    small_ctgs = [len(rec)/i for rec in records if len(rec) <kb1*i]
    mid_ctgs = [len(rec)/i for rec in records if kb1*i <= len(rec) < kb2*i]
    big_ctgs = [len(rec)/i for rec in records if kb2*i <= len(rec) <= kb3*i]
    oversized = [len(rec)/i for rec in records if len(rec) >=kb3*i]
    # write to report file
    ctg_stats_file = dirs['ctg_stats']+"contig_stats.txt"
    txt = ["\n"+g_name, str(len(small_ctgs)), str(len(mid_ctgs)),
           str(len(big_ctgs)), str(len(oversized))]
    open(ctg_stats_file, 'a').write("\t".join(txt))
    # make figure
    fname = dirs['ctg_stats']+g_name+"_ctg_stats.png"
    ctg_cats = small_ctgs, mid_ctgs, big_ctgs
    plot_ctg_stats(ctg_cats, fname)

def plot_ctg_stats(ctg_cats, fname):
    """Plot distribution statistics of contigs per genome."""
    small_ctgs, mid_ctgs, big_ctgs = ctg_cats
    kb1, kb2, kb3 = ctg_thresholds
    # make bins
    factor = 10
    small_bins = np.arange(0, kb1, kb1/factor)
    mid_bins = np.arange(kb1, kb2, (kb2-kb1)/factor)
    big_bins = np.arange(kb2, kb3, (kb3-kb2)/factor)
    # plot histograms
    fig = plt.figure(figsize=(16,6))
    if not len(small_ctgs) == 0:
        small_x = fig.add_subplot(131)
        small_x.hist(small_ctgs, small_bins, rwidth=kb1/factor)
        small_x.set_title("length < "+str(kb1)+" kb")
        small_x.set_ylabel("number of contigs")
    if not len(mid_ctgs) == 0:
        mid_x = fig.add_subplot(132)
        mid_x.hist(mid_ctgs, mid_bins, rwidth=(kb2-kb1)/factor)
        mid_x.set_title(str(kb1)+" kb < length < "+str(kb2)+" kb")
        mid_x.set_ylabel("number of contigs")
    if not len(big_ctgs) == 0:
        big_x = fig.add_subplot(133)
        big_x.hist(big_ctgs, big_bins, rwidth=(kb3-kb2)/factor)
        big_x.set_title("length > "+str(kb2)+" kb")
        big_x.set_ylabel("number of contigs")
    plt.savefig(fname)
    plt.clf()

def matches_table(run_id):
    """Compile table of contig matches."""
