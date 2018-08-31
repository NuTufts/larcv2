import ROOT as rt
from larcv import larcv
import sys
sys.argv.append("-b")


superafile = "larcv_062218.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "baka_cropinfill.root" )
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
DebugImage: true
MaxImages: 50
RandomizeCrops: true
MaxRandomAttempts: 50
MinFracPixelsInCrop: 0
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,split_cfg
fcfg.close()
split_pset = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )

# -------------------------------------
# UBLArFlowCropDetector

lfcrop_cfg="""Verbosity:0
InputBBoxProducer: \"detsplit\"
InputWireProducer: \"wire\"
InputLabelsProducer: \"Labels\"
InputADCProducer: \"ADC\"
OutputCroppedWireProducer: \"wire\"
OutputCroppedLabelsProducer: \"Labels\"
OutputCroppedADCProducer: \"ADC\"
OutputCroppedWeightsProducer: \"Weights\"
OutputCroppedMetaProducer: \"meta\"
OutputFilename: \"baka_cropinfill.root\"
CheckFlow: false
MakeCheckImage: false
DoMaxPool: false
RowDownsampleFactor: 2
ColDownsampleFactor: 2
MaxImages: 10
LimitOverlap: false
MaxOverlapFraction: 0.2
"""

lfcfg = open("ublarflowcrop.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "ublarflowcrop.cfg", "UBCropInfill" )

# -------------------------------------
# ALGOS

split_algo = larcv.UBSplitDetector()
split_algo.configure(split_pset)
split_algo.initialize()

lfcrop_algo = larcv.UBCropInfill()
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

lfcrop_algo.finalize()
io.finalize()


print "FIN"
