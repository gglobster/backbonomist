from os import path
from config import fixed_dirs, run_dirs, p_root_dir, project_id, \
    project_date, genomes, references, ctg_thresholds
from common import ensure_dir
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable, Size, AxesGrid


def log_start_run(run_id, timestamp):
    """Record run initiated in the dataset log."""
    set_log = p_root_dir+"/"+project_id+"_log.html"
    param_file = run_id+"/"+run_dirs['reports']+run_id+"_dataset.txt"
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
    set_log_htm = ["<li><b>", run_id, "</b>", "&nbsp;- resumed ", timestamp,
                   " (from step "+str(step)+")</li>"]
    linkline = "".join(set_log_htm)
    open(set_log, 'a').write(linkline)

def save_datasumm(run_id, timestamp):
    """Save a summary of the dataset composition to file."""
    print " ", run_id,
    report_root = p_root_dir+run_id+"/"+run_dirs['reports']
    param_file = report_root+run_id+"_dataset.txt"
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
        try:
            r_data.append("\t".join([ref['name'], str(len(ref['segs'])),
                                     ref['file']]))
            r_data.append("\t".join(["###", "Segments"]))
            r_data.append("\t".join(["###", "Name", "Start", "Stop"]))
            for seg in ref['segs']:
                r_data.append("\t".join(["", seg['name'],
                                         str(seg['coords'][0]),
                                         str(seg['coords'][1])]))
        except KeyError:
            r_data.append("\t".join([ref['name'], ref['file']]))
            r_data.append("\t".join(["###", "Segments generated from size",
                                     str(ref['chop_size'])]))
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
    ctg_stats_file = fixed_dirs['ctg_stats']+"contig_stats.txt"
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
    ctg_stats_file = fixed_dirs['ctg_stats']+"contig_stats.txt"
    txt = ["\n"+g_name, str(len(small_ctgs)), str(len(mid_ctgs)),
           str(len(big_ctgs)), str(len(oversized))]
    open(ctg_stats_file, 'a').write("\t".join(txt))
    # make figure
    fname = fixed_dirs['ctg_stats']+g_name+"_ctg_stats.png"
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

def matches_table(run_ref, run_id, ref_hits, ctl_scores):
    """Compile tables and graphs of contig matches."""
    # set imputs and outputs
    ref_n = run_ref.name
    run_root = p_root_dir+run_id+"/"
    report_root = run_root+run_dirs['reports']+ref_n+"/"
    hits_table_txt = report_root+run_id+"_"+ref_n+"_hits_table.txt"
    red_hits_fig = report_root+run_id+"_"+ref_n+"_sum_hits.png"
    mf_hits_fig = report_root+run_id+"_"+ref_n+"_all_hits.png"
    ensure_dir([report_root])
    rep_fhandle = open(hits_table_txt, 'w')
    rep_fhandle.write("# Matches to the reference segments of "+ref_n)
    segs = [seg['name'] for seg in run_ref.segs]
    red_g_list = []
    mf_g_list = []
    mf_ctgs = []
    g_names = []
    # traverse dict of results
    for g_name in sorted(ref_hits, reverse=True):
        genome_root = report_root+g_name+"/"
        ensure_dir([genome_root])
        g_table_fig = genome_root+g_name+"_hits_"+ref_n+".png"
        # collect scores
        g_list = []
        g_header = ["\n\n"+g_name] + segs
        g_lines = ["\t".join(g_header)]
        for contig in ref_hits[g_name]:
            ctg_hits = [ref_hits[g_name][contig][seg]
                        if seg in ref_hits[g_name][contig]
                        else 0 for seg in segs]
            g_list.append(ctg_hits)
            ctg_line = "\t".join([contig]+[str(score) for score in ctg_hits])
            g_lines.append(ctg_line)
        rep_fhandle.write("\n".join(g_lines))
        # normalize scores
        g_array = np.array(g_list)
        if not len(g_list):
            g_norm = np.zeros(len(ctl_scores))
            g_norm = np.reshape(g_norm, (1, len(ctl_scores)))
        else:
            try:
                g_norm = np.divide(g_array, ctl_scores)
            except ValueError:
                print "ERROR: problem processing hits for", g_name
                g_norm = np.zeros(len(ctl_scores))
                g_norm = np.reshape(g_norm, (1, len(ctl_scores)))
        # graph detailed scores per genome
        contigs = ref_hits[g_name].keys()
        hits_heatmap_multi(ref_n, segs, [g_name], [contigs], [g_norm],
                           g_table_fig)
        # prep for global graphs
        if len(g_array) > 1:
            g_reduce = np.amax(g_array, 0)
            g_reduce = np.divide(g_reduce, ctl_scores)
        else:
            g_reduce = np.reshape(g_norm, (len(ctl_scores),))
        red_g_list.append(g_reduce)
        mf_g_list.append(g_norm)
        mf_ctgs.append(contigs)
        g_names.append(g_name)
    # graph detailed scores for all genomes
    hits_heatmap_multi(ref_n, segs, g_names, mf_ctgs, mf_g_list, mf_hits_fig)
    # graph summarized scores for all genomes
    red_g_array = np.array(red_g_list)
    g_names.reverse()
    hits_heatmap_multi(ref_n, segs, [1], [g_names], [red_g_array],
                       red_hits_fig)
    # done
    rep_fhandle.close()

def hits_heatmap_multi(ref_n, segs, g_names, contigs, scores, imgfile):
    """Combine matches heatmaps in subplots."""
    g_names.reverse()
    contigs.reverse()
    scores.reverse()
    ctg_count = 0
    for item in contigs:
        ctg_count += len(item)
    fig = plt.figure(figsize=(len(segs)/2+1,ctg_count/2+3))
    grid = AxesGrid(fig, 111, nrows_ncols=(len(g_names), 1), axes_pad=0.4,
                    cbar_location="top", cbar_mode="single", cbar_size=0.1)
    for i in range(len(g_names)):
        hmap = grid[i].pcolor(scores[i], cmap='hot', vmin=0, vmax=1)
        grid[i].xaxis.set_major_locator(MaxNLocator(len(segs)))
        grid[i].yaxis.set_major_locator(MaxNLocator(len(scores[i])))
        grid.cbar_axes[i].colorbar(hmap)
        grid.cbar_axes[i].set_xticks([0, 0.5, 1])
        grid.cbar_axes[i].set_xticklabels(['Low', 'Medium', 'High'])
        grid[i].set_xticklabels(segs, size='small')
        grid[i].set_xlabel(ref_n+" reference segments", size='small')
        grid[i].set_title(g_names[i], size='small')
        grid[i].set_yticklabels("")
        labels = grid[i].get_xticklabels()
        for label in labels:
            label.set_rotation(30)
        y_index = 0.5
        for contig in contigs[i]:
            grid[i].text(-0.2, y_index, contig, size='small',
                         horizontalalignment='right',
                         verticalalignment='center',)
            y_index +=1
    plt.savefig(imgfile)
    plt.clf()

def hits_heatmap(title, segs, contigs, scores, imgfile):
    """Graph matches as a heatmap.

    Deprecated in favor of the multi-plot version.

    """
    fig = plt.figure()
    h = [Size.Fixed(1.5), Size.Fixed(15)]
    v = [Size.Fixed(2.5), Size.Fixed(5.)]
    ax = plt.subplot(111)
    hmap = ax.pcolor(scores, cmap='hot', vmin=0, vmax=1)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="5%", pad=1)
    cbar = plt.colorbar(hmap, cax=cax, ticks=[0, 0.5, 1],
                    orientation='horizontal')
    cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])
    ax.xaxis.set_major_locator(MaxNLocator(len(segs)))
    ax.yaxis.set_major_locator(MaxNLocator(len(contigs)))
    ax.set_xticklabels(segs, size='small')
    ax.set_yticklabels("")
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30)
    y_index = 0.5
    for contig in contigs:
        ax.text(-0.2, y_index, contig, size='small',
                horizontalalignment='right',
                verticalalignment='center',)
        y_index +=1
    ax.set_xlabel("Reference segments", size='small')
    ax.set_title(title)
    plt.savefig(imgfile)
    plt.clf()
