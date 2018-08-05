import ROOT as rt
from larcv import larcv
import sys
sys.argv.append("-b")


superafile = "/home/jmills/workdir/uresnet/larflow/larcv/larcv/app/UBImageMod/test/larcv.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "croppedsegment.root" )
io.initialize()

# -------------------------------------
# UBSplitDetector

split_cfg="""Verbosity:0
InputProducer: \"wire\"
OutputBBox2DProducer: \"detsplit\"
CropInModule: true
OutputCroppedProducer: \"detsplit\"
BBoxPixelHeight: 512
BBoxPixelWidth: 832
CoveredZWidth: 310
FillCroppedYImageCompletely: true
DebugImage: false
MaxImages: 10
RandomizeCrops: true
MaxRandomAttempts: 50
MinFracPixelsInCrop: 0.0001
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,split_cfg
fcfg.close()
split_pset = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )

# -------------------------------------
# UBSegmentCropDetector

lfcrop_cfg="""Verbosity:0
InputBBoxProducer: \"detsplit\"
InputADCProducer: \"wire\"
InputLabelsProducer: \"segment\"
InputCroppedADCProducer: \"detsplit\"
OutputCroppedADCProducer: \"adc\"
OutputLabelsProducer: \"labels\"
OutputWeightsProducer: \"weights\"
OutputCroppedMetaProducer: \"flowmeta\"
OutputFilename: \"croppedsegment_lf.root\"
CheckFlow: true
MakeCheckImage: true

MaxImages: 10
LimitOverlap: false
MaxOverlapFraction: 0.2

BackgroundBias: 0
TrackBias:0
TrackEndBias: 0
ShowerBias: 0
"""

lfcfg = open("ubsegmentcrop.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "ubsegmentcrop.cfg", "UBCropSegment" )

# -------------------------------------
# ALGOS

split_algo = larcv.UBSplitDetector()
split_algo.configure(split_pset)
split_algo.initialize()

lfcrop_algo = larcv.UBCropSegment()
lfcrop_algo.configure(lfpset)
lfcrop_algo.initialize()

# -------------------------------------

nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 10
#Number of entries total: nentries * MaxImages * 3 (planes)
for n in range(0,nentries):
    io.read_entry(n)
    split_algo.process( io )
    lfcrop_algo.process( io );

lfcrop_algo.finalize()
io.finalize()


print "FIN"
