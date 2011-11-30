from array_tetris import coord_chop
from Bio.SeqFeature import SeqFeature, FeatureLocation
from writers import write_fasta

class Reference(object):
    """Persistent reference object."""

    def __init__(self, name, file, input, seg_mode, fas_out, gbk_out,
                 seg_out, logfile):
        self.name = name
        self.file = file
        self.input = input
        self.seg_mode = seg_mode
        self.segs = []
        self.fas = fas_out
        self.gbk = gbk_out
        self.segs_dir = seg_out
        self.log = logfile

    def get_segs_from_list(self, list):
        self.segs = list

    def get_segs_from_chop(self, length, size):
        pair_list = coord_chop(length, size, 'exact_size')
        counter = 0
        for a,b in pair_list:
            counter +=1
            seg = {'coords': (a, b), 'name': str(counter), 'note': str(a)+'_'+str(b)}
            self.segs.append(seg)

    def extract_segs_seqs(self, record, out_dir):
        count = 0
        for seg in self.segs:
            # unpack segment coords
            seg_start, seg_stop = seg['coords'][0], seg['coords'][1]
            # extract segment sequence
            segment = record[seg_start:seg_stop]
            segment.id = self.name+"_"+seg['name']
            # write to individual file
            out_file = out_dir+self.name+"_"+seg['name']+".fas"
            write_fasta(out_file, segment)
            # record segment feature
            feat_loc = FeatureLocation(seg_start, seg_stop)
            feature = SeqFeature(location=feat_loc,
                                 type='ref_seg',
                                 qualifiers={'id': seg['name']})
            record.features.append(feature)
            count +=1
        return record
