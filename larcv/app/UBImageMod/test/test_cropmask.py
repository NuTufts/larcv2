import ROOT as rt
from ROOT import TFile, TTree
from array import array
from larcv import larcv
import sys
sys.argv.append("-b")


superafile = "larcv.root"

io = larcv.IOManager(larcv.IOManager.kBOTH)
io.add_in_file( superafile )
io.set_out_file( "croppedmask.root" )
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
MaxRandomAttempts: 60
MinFracPixelsInCrop: 0.0001
"""

fcfg = open("ubsplit.cfg",'w')
print >>fcfg,split_cfg
fcfg.close()
split_pset = larcv.CreatePSetFromFile( "ubsplit.cfg", "UBSplitDetector" )

# -------------------------------------
# UBMaskCropDetector

maskcrop_cfg="""Verbosity:2
InputBBoxProducer: \"detsplit\"
InputADCProducer: \"wire\"
InputMasksProducer: \"cluster\"
InputCroppedADCProducer: \"detsplit\"
OutputCroppedADCProducer: \"adc\"
OutputMasksProducer: \"masks\"
OutputWeightsProducer: \"weights\"
OutputCroppedMetaProducer: \"flowmeta\"
OutputFilename: \"croppedmask_lf.root\"
CheckFlow: true
MakeCheckImage: true

MaxImages: 10
LimitOverlap: false
MaxOverlapFraction: 0.2


"""

maskcfg = open("ubmaskcrop.cfg",'w')
print >>maskcfg,maskcrop_cfg
maskcfg.close()
maskpset = larcv.CreatePSetFromFile( "ubmaskcrop.cfg", "UBCropMask" )

# -------------------------------------
# ALGOS

split_algo = larcv.UBSplitDetector()
split_algo.configure(split_pset)
split_algo.initialize()

maskcrop_algo = larcv.UBCropMask()
maskcrop_algo.configure(maskpset)
maskcrop_algo.initialize()

# -------------------------------------
# TTREE for Mask Class Count
ClassCounts = TTree( 'ClassCounts_Tree', 'Tree to hold Counts for Each Class')
count_cosmic = array('I', [0])
count_electron = array('I', [0])
count_neutron = array('I', [0])
count_proton = array('I', [0])
count_neutrino = array('I', [0])
count_other = array('I', [0])

ClassCounts.Branch('count_cosmics', count_cosmic, 'count_cosmics/I')
ClassCounts.Branch('count_electrons', count_electron,  'count_electron/I')
ClassCounts.Branch('count_neutrons', count_neutron,  'count_neutron/I' )
ClassCounts.Branch('count_protons', count_proton,  'count_proton/I' )
ClassCounts.Branch('count_neutrinos', count_neutrino,  'count_neutrino/I' )
ClassCounts.Branch('count_others', count_other,  'count_other/I' )

# -------------------------------------
# Loop
nentries = io.get_n_entries()
print "Num Entries: ",nentries
# nentries = 10
#Number of entries total: nentries * MaxImages * 3 (planes)
for n in range(0,nentries):
    print "Entry:   ", n
    io.read_entry(n)
    split_algo.process( io )
    maskcrop_algo.process( io );
    count_cosmic[0] += maskcrop_algo.get_count_cosmic()
    # print "cosmic:   ", maskcrop_algo.get_count_cosmic()
    count_electron[0] += maskcrop_algo.get_count_electron()
    # print "electron:   ", maskcrop_algo.get_count_electron()
    count_proton[0] += maskcrop_algo.get_count_proton()
    # print "proton:   ", maskcrop_algo.get_count_proton()
    count_neutron[0] += maskcrop_algo.get_count_neutron()
    # print "neutron:   ", maskcrop_algo.get_count_neutron()
    count_neutrino[0] += maskcrop_algo.get_count_neutrino()
    # print "neutrino:   ", maskcrop_algo.get_count_neutrino()
    count_other[0] += maskcrop_algo.get_count_other()
    # print "other:   ", maskcrop_algo.get_count_other()

ClassCounts.Fill()
ClassCounts.Write()

maskcrop_algo.finalize()
io.finalize()


print "FIN"
