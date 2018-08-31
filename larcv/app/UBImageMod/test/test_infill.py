import ROOT as rt
from larcv import larcv


superafile = "larcv_062218.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "baka_infill.root" )
io.initialize()

# -------------------------------------
# UBSplitDetector

scfg="""Verbosity: 0
InputProducer: \"wire\"
LabelsInput: \"Labels\"
ADCInput: \"ADC\"
OutputBBox2DProducer: \"detsplit\"
CropInModule: true
OutputCroppedProducer: \"cwire\"
OutputLabelsProducer: \"cLabels\"
OutputADCProducer: \"cADC\"
OutputWeightsProducer: \"cWeights\"
BBoxPixelHeight: 512
BBoxPixelWidth: 832
CoveredZWidth: 310
FillCroppedYImageCompletely: true
DebugImage: false
MaxImages: 10
RandomizeCrops: true
MaxRandomAttempts: 50
MinFracPixelsInCrop: 0
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,scfg
fcfg.close()


cfg = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector_Infill" )
algo = larcv.UBSplitDetector_Infill()
algo.configure(cfg)
algo.initialize()

# -------------------------------------

cfgcrop="""
"""


nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 1

for n in range(nentries):
    io.read_entry(n)
    algo.process( io )

    detsplit = io.get_data( "image2d", "detsplit" )

    h_v = {}
    for i in range(0,detsplit.image2d_array().size()):
        h_v[i] = larcv.as_th2d( detsplit.image2d_array()[i], "test%d"%(i) )
    print h_v

    c = rt.TCanvas("c","c",1500,400)
    c.Divide(3,1)
    for i in range(len(h_v)/3):
        for j in range(3):
            c.cd(j+1)
            h_v[3*i+j].Draw("COLZ")

        c.Update()
        raw_input()
    io.save_entry()

io.finalize()

print "FIN"
