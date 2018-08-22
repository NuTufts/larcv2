import ROOT as rt
from larcv import larcv


#superafile = "/media/hdd2/taritree/larflow/xfer/larcv_5477923_0.root"
superafile = "../../../../..//testdata/larcv_5482426_95.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "baka.root" )
io.initialize()

# -------------------------------------
# UBSplitDetector

scfg="""Verbosity: 0
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
MaxRandomAttempts: 50
MinFracPixelsInCrop: 0.0001
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,scfg
fcfg.close()


cfg = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )
algo = larcv.UBSplitDetector()
algo.configure(cfg)
algo.initialize()

# -------------------------------------

cfgcrop="""
"""


nentries = io.get_n_entries()
print "Num Entries: ",nentries
nentries = 3

for n in range(nentries):
    io.read_entry(n)
    algo.process( io )

    detsplit = io.get_data( "image2d", "detsplit" )

    if False:
        # visualize
        h_v = {}
        for i in range(0,detsplit.image2d_array().size()):
            h_v[i] = larcv.as_th2d( detsplit.image2d_array()[i], "test%d"%(i) )
            h_v[i].GetZaxis().SetRangeUser(-10,100)
        #print h_v

        c = rt.TCanvas("c","c",1500,400)
        c.Divide(3,1)
        for i in range(len(h_v)/3):
            for j in range(3):
                c.cd(j+1)
                h_v[3*i+j].Draw("COLZ")

            c.Update()
            raw_input()
    print "save entry"
    io.save_entry()

algo.printElapsedTime()
io.finalize()

io2 = larcv.IOManager(larcv.IOManager.kREAD)
io2.add_in_file( "baka.root" )
io2.initialize()

for i in xrange(io2.get_n_entries()):
    io2.read_entry(i)
    evdetsplit = io2.get_data("image2d","detsplit")
    print "entry ",i,": desplit images = ",evdetsplit.as_vector().size()

print "FIN"
