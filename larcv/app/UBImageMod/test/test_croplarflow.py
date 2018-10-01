import ROOT as rt
from larcv import larcv
import sys
sys.argv.append("-b")


#superafile = "/media/hdd2/taritree/larflow/xfer/larcv_5477923_0.root"
superafile = "../../../../../testdata/larcv_5482426_95.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "baka.root" )
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
MaxImages: -1
RandomizeCrops: false
MaxRandomAttempts: 1000
MinFracPixelsInCrop: 0.0
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,split_cfg
fcfg.close()
split_pset = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )

# -------------------------------------
# UBLArFlowCropDetector

lfcrop_cfg="""Verbosity:0
InputBBoxProducer: \"detsplit\"
InputCroppedADCProducer: \"detsplit\"
InputADCProducer: \"wire\"
InputVisiProducer: \"pixvisi\"
InputFlowProducer: \"pixflow\"
OutputCroppedADCProducer:  \"adc\"
OutputCroppedVisiProducer: \"visi\"
OutputCroppedFlowProducer: \"flow\"
OutputCroppedMetaProducer: \"flowmeta\"
OutputFilename: \"baka_lf.root\"
SaveOutput: false
CheckFlow:  true
MakeCheckImage: true
DoMaxPool: false
RowDownsampleFactor: 2
ColDownsampleFactor: 2
MaxImages: -1
LimitOverlap: false
RequireMinGoodPixels: false
MaxOverlapFraction: 0.2
IsMC: true
UseVectorizedCode: true
"""

lfcfg = open("ublarflowcrop.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "ublarflowcrop.cfg", "UBCropLArFlow" )

# -------------------------------------
# ALGOS

split_algo = larcv.UBSplitDetector()
split_algo.configure(split_pset)
split_algo.initialize()

lfcrop_algo = larcv.UBCropLArFlow()
lfcrop_algo.configure(lfpset)
lfcrop_algo.initialize()

# -------------------------------------

nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 1

for n in range(0,nentries):
    io.read_entry(n)
    split_algo.process( io )
    lfcrop_algo.process( io );
    print "finished event"
    break
    
lfcrop_algo.finalize()
io.finalize()


print "FIN"
