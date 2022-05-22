import sys
import time
import numpy as np
from ROOT import TFile, TTree
from ROOT import gROOT, gSystem

t0 = time.time()

dbname = "pulseA"
if len(sys.argv)>1:
    dbname = sys.argv[1]

pulse = np.load("npy/"+dbname+".npy")
pulse = pulse.astype("float64") # float16 is not supported in pyroot
print("Pulse Data: ")
print(" data type ",type(pulse))
print(" data dtype ",pulse.dtype)
print(" data size ",pulse.size)
print(" data shape ",pulse.shape)
print(" data ndim ",pulse.ndim)

print("load pulse done\n")

position = np.load("npy/"+dbname+"pos.npy")

print("Positions Data: ")
print(" data type ",type(position))
print(" data dtype ",position.dtype)
print(" data size ",position.size)
print(" data shape ",position.shape)
print(" data ndim ",position.ndim)

print("load positions done\n")

spulsepoints = pulse.shape[1]*(pulse.shape[2]-1)

fout = TFile("pulsedb/"+dbname+".root","RECREATE")
tree = TTree("tree","pulse shape db")
seg    = np.empty((1), dtype="int")
pos    = np.empty((3), dtype="double")
core   = np.empty((pulse.shape[1]), dtype="double")
spulse = np.empty((spulsepoints), dtype="double")
tree.Branch("seg", seg, "seg/I")
tree.Branch("pos", pos, "pos[3]/D")
tree.Branch("core", core, "core[%d]/D" % pulse.shape[1])
tree.Branch("spulse", spulse, "spulse[%d]/D" % spulsepoints)

print("convert npy/%s.npy to pulsedb/%s.root\n" % (dbname,dbname))
spulsetmp=np.empty((spulsepoints), dtype="double")

for ipoint in range(pulse.shape[0]):
    if ipoint%1000==0 :
        print("\r finish %d / %d points..." % (ipoint,pulse.shape[0]), end = ' ', flush=True)
    
    pos[0:3]=position[ipoint].copy()
    core[0:pulse.shape[1]]=pulse.transpose(0,2,1)[ipoint][0][0:pulse.shape[1]]
    for i in range(pulse.shape[2]-1):
        spulsetmp[i*pulse.shape[1]:(i+1)*pulse.shape[1]]=pulse.transpose(0,2,1)[ipoint][i+1][0:pulse.shape[1]]

    # segment of max pulse
    seg[0] = int(np.where(spulsetmp==np.amax(spulsetmp))[0][0]/pulse.shape[1])+1

    spulse[0:spulsepoints] = spulsetmp
    tree.Fill()

print("\r finish %d / %d points...\n" % (ipoint,pulse.shape[0]))
        
fout.cd()
tree.Write()
fout.Close()
    
print("============ Elapsed time: %.1f seconds =============" % (time.time()-t0))
