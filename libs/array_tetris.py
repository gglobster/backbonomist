import numpy as np
from config import mtype

def extract_nonzero(array_o, array_d):
    """Extract rows that are non-zero into an array of same shape."""
    for row in array_o:
        a, b, c, d = row
        if a == 0 or c == 0:
            pass
        else:
            array_d = np.append(array_d, row)
    return array_d

def clump_rows(array, dist):
    """Collapse rows of coordinates that are close and in same orientation."""
    # sort array to ensure optimal clumping
    array = np.sort(array, order=mtype[0][0])
    # perform clumping
    r_ix = 1
    new_len = len(array)
    while r_ix < new_len:
        xa1, xb1, xc1, xd1 = array[r_ix-1]
        xa2, xb2, xc2, xd2 = array[r_ix]
        # clump if test conditions are fulfilled
        if float(xa1)/float(xc1) > 0 : # same sign
            #print xa2, xb1, xc2, xd1
            #print abs(xa2)-abs(xb1), abs(xc2)-abs(xd1)
            if abs(abs(xa2)-abs(xb1)) < dist \
            and abs(abs(xc2)-abs(xd1)) < dist:
                #print "case A yes", r_ix
                xa,xb = xa1,xb2
                xc,xd = xc1,xd2
                new_row = xa, xb, xc, xd
                # delete the original rows
                #print "before del", len(array)
                array = np.delete(array, (r_ix-1, r_ix), 0)
                #print "after del", len(array)
                # insert the updated row
                array = np.insert(array, r_ix-1, new_row)
                #print "after append", len(array)
            else :
                #print "case A no", r_ix
                r_ix += 1
        elif float(xa1)/float(xc1) < 0 : # different sign
            # (!!! may be broken for N>2) TODO: fix this!
            if abs(abs(xa1)-abs(xb2)) < dist \
            and abs(abs(xc2)-abs(xd1)) < dist:
                #print "case B yes", r_ix
                xa,xb = xa2,xb1
                xc,xd = xc1,xd2
                new_row = xa,xb,xc,xd
                # delete the original rows
                #print "before del", len(array)
                array = np.delete(array, (r_ix-1, r_ix), 0)
                #print "after del", len(array)
                # insert the updated row
                array = np.insert(array, r_ix-1, new_row)
                #print "after append", len(array)
            else :
                #print "case B no", r_ix
                r_ix += 1
        else : # this should not happen
            print "WTF? This shouldn't happen."
            r_ix += 1
        new_len = len(array)
    return array

def get_anchor_loc(quad_array):
    """Determine which segment is largest and where it hits."""
    # set up an array stub to receive size and location data
    sl_array = np.ones(1, dtype=[('size', 'i4'),
                                 ('start', 'i4'),
                                 ('end', 'i4'),
                                 ('orient', 'i2')])
    # collect data
    for row in quad_array:
        a, b, c, d = row
        size = -abs(abs(a)-abs(b)) # size negative'd so biggest shows up first
        loc_start = abs(a)
        loc_end = abs(b)
        orient = abs(a)/a
        data = size, loc_start, loc_end, orient
        sl_array = np.insert(sl_array, 0, data)
    # evaluate which segment is largest
    sl_array = np.sort(sl_array, order='size')
    anchor = sl_array[0]
    return anchor

