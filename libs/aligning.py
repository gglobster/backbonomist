import subprocess
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from config import mauve_exec

def align_clustal(file_name):
    """Make external call to ClustalW aligner."""
    cline = ClustalwCommandline("clustalw", infile=file_name)
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse ClustalW errors
    return report

def align_muscle(infile_name, outfile_name, log_file):
    """Make external call to Muscle aligner."""
    cline = MuscleCommandline(input=infile_name, out=outfile_name, clw=True,
                              loga=log_file, quiet='y')
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse MUSCLE errors
    return report

def align_mauve(file_list, output):
    """Make external call to Mauve aligner."""
    input_files = ' '.join(file_list)
    cline = mauve_exec+" --output="+ output +" "+input_files
    child = subprocess.Popen(str(cline), stdout=subprocess.PIPE, shell=True)
    output, error = child.communicate()
    report = {'output': output, 'error': error}
    # TODO: should set up something to parse Mauve errors
    return report