import ROOT as rt
from larcv import larcv
import sys
sys.argv.append("-b")


superafile = "output_infill.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "infill.root" )
io.initialize()

# -------------------------------------
# UBLArFlowCropDetector

lfcrop_cfg="""Verbosity:0
InputWireProducer: \"wire\"
InputNoFillProducer: \"nofill\"
InputFillProducer: \"fill\"
InputWeightsProducer: \"weights\"
OutputOverlayProducer: \"overlay\"
OutputInfillProducer: \"Infill_outputImagesProcessFactory\"
OutputFilename: \"infill.root\"
"""

lfcfg = open("finalinfill.cfg",'w')
print >>lfcfg,lfcrop_cfg
lfcfg.close()
lfpset = larcv.CreatePSetFromFile( "finalinfill.cfg", "InfilloutputImages" )

# -------------------------------------
# ALGOS

lfcrop_algo = larcv.InfilloutputImages()
lfcrop_algo.configure(lfpset)
lfcrop_algo.initialize()

# -------------------------------------

nentries = io.get_n_entries()
print "Num Entries: ",nentries

for n in range(0,nentries):
    io.read_entry(n)
    split_algo.process( io )
    lfcrop_algo.process( io );

lfcrop_algo.finalize()
io.finalize()


print "FIN"
