from config import directories as dirs, genomes
from common import ensure_dir
from drawing import ContigDraw

def map_scaffolds(contig):
    """Generate annotated maps of backbone scaffolds."""
    # set inputs and outputs
    ctg_name = contig['name'] # reference contig
    scaff_annot_root = dirs['scaff_annot_dir']+ctg_name+"/genbank/"
    maps_out_root = dirs['maps_dir']+ctg_name+"/"
    ensure_dir(maps_out_root)
    print " ", ctg_name
    # cycle through genomes
    for genome in genomes:
        # set inputs and outputs
        g_name = genome['name']
        fin_gbk_out = scaff_annot_root+g_name+"_"+ctg_name+"_fin.gbk"
        map_file = maps_out_root+g_name+"_"+ctg_name+".pdf"
        print '\t', g_name
        # generate graphical map
        ContigDraw(g_name, fin_gbk_out, map_file)