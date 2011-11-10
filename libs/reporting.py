from os import path
from config import directories as dirs, p_root_dir, project_id, \
    project_date, genomes, references, blast_prefs, max_size
from common import ensure_dir

def log_start_run(run_id, timestamp):
    """Record run initiated in the dataset log."""
    set_log = p_root_dir+"/"+project_id+"_log.html"
    param_file = run_id+"/"+dirs['reports']+run_id+"_parameters.txt"
    header = "<p><b>Log of processing runs for "+project_id+"</b></p><p><ul>"
    set_log_htm = ["<li><b>", run_id, "</b>&nbsp;- initiated ", timestamp,
                   " (<a href='", param_file, "'>parameters</a>)</li>"]
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

def save_datasumm(run_id, timestamp):
    """Save a summary of the dataset composition to file."""
    print " ", run_id,
    report_root = p_root_dir+run_id+"/"+dirs['reports']
    param_file = report_root+run_id+"_parameters.txt"
    ensure_dir([report_root])
    # genome data
    g_title = "\t".join(["##", "Genomes"])
    g_header = "\t".join(["##", "Name", "Offset", "Format", "File"])
    g_data = [g_title, g_header]
    for genome in genomes:
        g_data.append("\t".join(["", genome['name'], str(genome['offset'][1]),
                                 genome['input'], genome['file']]))
    g_block = "\n".join(g_data)
    # references data
    r_title = "\t".join(["##", "References"])
    r_header = "\t".join(["##", "Name", "Segments", "File"])
    r_data = [r_title, r_header]
    for ref in references:
        r_data.append("\t".join(["", ref['name'], str(len(ref['refs'])),
                                 ref['file']]))
        r_data.append("\t".join(["", "###", "Segments"]))
        r_data.append("\t".join(["", "###", "Name", "Start", "Stop"]))
        for seg in ref['refs']:
            r_data.append("\t".join(["", "", seg['name'],
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
    