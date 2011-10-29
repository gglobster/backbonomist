from Bio import GenBank
from reportlab.pdfgen import canvas
from reportlab.lib.units import cm
from reportlab.lib.colors import black, white, HexColor
from loaders import load_genbank
from config import fct_flags, fct_colors

### Hardcoded presets ###
## General proportions ##
u 		= 0.0125		# conversion factor for sequence length-related values (0.0125, increase for short sequences)
hmar 	= 1*cm			# horizontal margin to canvas (1*cm)
vmar 	= 2*cm			# vertical margin to canvas (2*cm)
pNsize 	= 4*cm			# width to set aside for plasmid names (2.5*cm)
di 		= 0.18*cm		# half-length of interruption/frameshift tick marks (0.18*cm)
doL 	= 0.7*cm		# distance of ORF labels from their ORF (0.7*cm)
minL 	= 700			# minimum ORF size (in base pairs) for 'full arrows' (700, decrease for short seq)
w 		= 0.3*cm     	# half-width of the tail (0.3*cm)
h 		= 0.175*cm		# distance from the side tips of the head to the neck (0.175*cm)
dBL	 	= 3*cm			# distance between plasmid baselines (4*cm)
da 		= 0.4*cm		# distance of alignment tick marks from the corresponding plasmid baseline (0.4*cm)
ck_hsp 	= 170		    # horizontal spacing between columns (170)
ck_vsp 	= 17		    # vertical spacing between items (17)

## Typeface, fonts, sizes ##
rFont 	= "Helvetica"			# regular type ("Helvetica")
bFont 	= "Helvetica-Bold"		# bold type ("Helvetica-Bold")
LfSize 	= 14					# large font size (14)
NfSize 	= 12					# normal font size (12)
SfSize 	= 10					# small font size (10)

## Special feature settings ##
SFX 	= 'off'			# allow drawing special sequence features beside CDS/ORF ('on' or 'off')
osym 	= '*'			# symbol for origin of sequence (*)
snp		= '*'			# symbol for SNP (*)

## Scale of sequence, in legend ##
scX 	= 0*u			# starting point on the X axis(0)
incrT 	= 1000*u		# increment size in basepairs (1000, signifying 1 kb)
incrN 	= 5				# increment number (5, or 1 for short sequences)
dip 	= -2*cm			# starting point on Y axis (-2*cm)
dop 	= 0.3*cm		# size of vertical tick

## Canvasser : Initialize canvas ##
def Canvasser(hCan,vCan,transX,transY,outfile) :
    canvasN = canvas.Canvas(outfile, pagesize=(hCan,vCan))
    canvasN.translate(transX,transY)
    canvasN.setStrokeColor(black)
    canvasN.setFillColor(white)
    canvasN.setLineWidth(1)
    canvasN.setLineJoin(1)
    canvasN.setLineCap(0)
    return canvasN

## BaseDraw : draws contig baselines and features ##
def BaseDraw(canvas, pName, cLen, feat_list) :
    # draw plasmid baseline
    BaseL(cLen, canvas)
    # label the baseline with plasmid name and size
    LabeL(pName, cLen, canvas)
    # filter and draw annotation features
    ORFcnt = 0
    Y_pos = 0
    for feature in feat_list :
        if feature.type == 'CDS' or feature.type == 'cds':
            ORFcnt += 1
            Y_pos = ORFeus(canvas, feature, ORFcnt, Y_pos)
        if feature.type == 'contig':
            TigTicker(canvas, feature)

## BaseL : Draws phage baselines ##
def BaseL(cLen, canvas) :
    canvas.setLineWidth(3)
    y0 = 0
    x0 = 0              # all sequences are aligned on the left
    x1 = cLen*u         # convert sequence length to workspace dimensions
    pBL = canvas.beginPath()
    pBL.moveTo(x0,y0)
    pBL.lineTo(x1,y0)
    canvas.drawPath(pBL, stroke=1, fill=0)
    pBL.close()

## LabeL : Labels baselines with plasmid name and size ##
def LabeL(cName, cLen, canvas) :
    canvas.setFillColor(black)
    y0 = 0
    x0 = -pNsize                # retreat into the left margin to write out name and size
    y1 = y0 + ck_vsp/10         # bump name up a bit from BL level
    y2 = y0 - ck_vsp            # bump size down a bit from BL level
    pLenStr = str(float(cLen)/1000) # converts number to suitable form
    canvas.setFont(bFont,LfSize)
    canvas.drawString(x0,y1,cName)
    canvas.setFont(rFont,NfSize)
    canvas.drawString(x0,y2,pLenStr+' kb')

## ORFeus : Draws plasmid ORFs ##
def ORFeus(canvas, feature, count, Y_pos) :
    canvas.setLineWidth(1)
    # extract relevant properties of the CDS feature
    flag = feature.qualifiers.get('flag') # not presently in use
    shape = feature.qualifiers.get('shape')
    annotation = feature.qualifiers.get('fct')
    # evaluate what strand the ORF is on (determines direction of arrow)
    cstrand = feature.strand
    if cstrand == None:
        cstrand = 1
    # take start and end points
    location 	= feature.location
    Zs 			= location.nofuzzy_start
    Ze 			= location.nofuzzy_end
    featL 		= Ze - Zs
    # calculate X axis coordinates (expr of cstrand has changed!)
    if cstrand == -1 :	# reverse orientation
        xs,xe 	= Ze*u,Zs*u		# start and end
        xn 		= xe+minL*u		# neck of arrow
    else :				# forward orientation
        xs,xe 	= Zs*u,Ze*u		# start and end
        xn 		= xe-minL*u		# neck of arrow
    midLZ 	= ((Zs+Ze)/2)*u	# middle of ORF for optional label
    # set Y axis coordinates
    y0 = 0
    yt,yb,ynt,ynb = y0+w,y0-w,y0+h,y0-h
    # initialize path
    pORF = canvas.beginPath()
    # draw square-shaped ORFS
    if shape == ['square'] :
        pORF.moveTo(xs,ynt)
        pORF.lineTo(xe,ynt)
        pORF.lineTo(xe,ynb)
        pORF.lineTo(xs,ynb)
        pORF.lineTo(xs,ynt)
    # draw triangle-shaped ORFS
    elif featL <= minL :
        pORF.moveTo(xs,yt)
        pORF.lineTo(xe,y0)
        pORF.lineTo(xs,yb)
        pORF.lineTo(xs,yt)
    # draw arrow-shaped ORFS
    elif featL > minL :
        pORF.moveTo(xs,ynt)
        pORF.lineTo(xn,ynt)
        pORF.lineTo(xn,yt)
        pORF.lineTo(xe,y0)
        pORF.lineTo(xn,yb)
        pORF.lineTo(xn,ynb)
        pORF.lineTo(xs,ynb)
        pORF.lineTo(xs,ynt)
    # evaluate function category and set fill color
    fct_key = annot_color(fct_flags, annotation)
    color_hex = HexColor(fct_colors[fct_key][0])
    canvas.setFillColor(color_hex)
    # finalize object path
    canvas.drawPath(pORF, stroke=1, fill=1)
    pORF.close()
    canvas.setFillColor(black)
    canvas.setFont(rFont, SfSize)
    # write CDS numbers
    canvas.drawCentredString(midLZ,y0-doL,str(count))
    # write annotation line
    if annotation[0] == 'no match':
        pass
    else:
        y_annot = y0-doL-dop*12-(Y_pos*ck_vsp)
        canvas.drawString(-pNsize, y_annot, str(count)+'. '+annotation[0])
        canvas.setFont(rFont, NfSize)
        Y_pos +=1
    return Y_pos

def TigTicker(canvas, feature):
    """Draw contig separators."""
    # get contig name
    name = feature.qualifiers.get('id')[0]
    # take start and end points
    location 	= feature.location
    Zs 			= location.nofuzzy_start
    Ze 			= location.nofuzzy_end
    xs,xe 	= Zs*u,Ze*u		# start and end
    # set Y axis coordinates
    y0 = h+dop
    # draw
    canvas.setLineWidth(2)
    ttl = canvas.beginPath()
    ttl.moveTo(xs,y0)
    ttl.lineTo(xs,y0+h*2.5)
    ttl.lineTo(xs+dop,y0+h*2.5)
    canvas.drawPath(ttl, stroke=1, fill=0)
    canvas.drawString(xs+dop*2,y0+h*2,name)
    ttl.close()

def SeqScale(canvas, scX, incrT, incrN, dip, dop) :
    """Draws the sequence scale bar."""
    canvas.setLineWidth(1.2)
    canvas.setFillColor(black)
    incrCNT = 0							# initialize count of increments
    psc = canvas.beginPath()
    psc.moveTo(scX,dip-dop)				# start at beginning (duh!)
    psc.lineTo(scX+incrT*incrN,dip-dop)	# draw the scale baseline
    # draw ticks until the max number of increments is reached
    while incrCNT <= incrN :
        psc.moveTo(scX+incrT*incrCNT,dip-dop)
        psc.lineTo(scX+incrT*incrCNT,dip)
        incrCNT += 1
    canvas.drawPath(psc, stroke=1, fill=0)
    psc.close()
    # write out scale extremities values (needs hand-fix if not using kbs)
    canvas.setFont(rFont,NfSize)
    canvas.drawRightString(scX,dip+dop,'0')
    canvas.drawString(scX+incrT*incrN,dip+dop,str(incrN)+' kb')

def ContigDraw(cName, in_file, out_file):
    """Draw sequence map of a single contig to file."""
    # load contig record
    seq_record = load_genbank(in_file)
    ctg_length = len(seq_record.seq)
    features = seq_record.features
    cds = [feature for feature in features
           if feature.type == 'CDS' or feature.type == 'cds']
    annot_cds = [1 for feature in cds
                 if feature.qualifiers.get('fct')[0] != 'no match']
    annot_cnt = sum(annot_cds)
    # calculate main canvas dimensions
    if ctg_length < 25000:
        hCan = 32*cm
    else:
        hCan = hmar*2 + pNsize + ctg_length*u
    vCan = dBL + vmar*2 + annot_cnt*ck_vsp
    transX = hmar + pNsize
    transY = dBL + vmar + annot_cnt*ck_vsp
    # set up main canvas
    canvas_main = Canvasser(hCan, vCan, transX, transY, out_file)
    # draw contig baseline and features
    BaseDraw(canvas_main, cName, ctg_length, features)
    # draw scale
    SeqScale(canvas_main, scX, incrT, incrN, dip, dop )
    # write to file and finalize the figure
    canvas_main.showPage()
    canvas_main.save()

def annot_color(fct_flags, annotation):
    """Look up the color to use based on annotation keywords."""
    annot_line = annotation[0].lower()
    fct_key = 'oth'
    for key in fct_flags:
        i = 0
        while i < len(fct_flags[key]):
            if annot_line.find(fct_flags[key][i]) > -1:
                fct_key = key
                break
            i +=1
    return fct_key

def PairwiseDraw(g_name, cstrct_gbk, ref_ctg_file, segdata, map_file):
    """Draw pairwise alignment map with similarity shading."""
    print "map"