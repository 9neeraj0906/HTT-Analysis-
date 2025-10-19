import ROOT 
import math
import time

#import data
dataPath = "/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root"

# Load Tree
f = ROOT.TFile.Open(dataPath)
tree = f.Get("Events")

# create histogram
# Syntax: ROOT.TH1F//1D histogram of Floats//("", "Graph title, xaxis, yaxis", bins, min, max)
hMt = ROOT.TH1F("hMtData", "Transverse mass Data; MT(Mu,MET) [GeV]; Events",100, 0.0, 150.0 )


# get the number of entries
numEntries = tree.GetEntries()
print("Total events:", numEntries)

start = time.time()

for i in range(numEntries):
	tree.GetEntry(i)

	#condition for atleast onemuon
	if tree.nMuon <= 0:
		continue
	#Initialise muon
	muPt = tree.Muon_pt[0]
	muEta = tree.Muon_eta[0]
	muPhi = tree.Muon_phi[0]
	if tree.nTau <=0:
		continue

	tauPt = tree.Tau_pt[0]


	# Muon should pass either of the isomu22 or mu19tau20
	passMuon = bool(tree.Muon_tightId[0])
	passIsoMu22 = passMuon and abs(muEta) <  2.1 and muPt >23
	passMu19Tau20 = (tree.nTau > 0 and 20<muPt<23 and tauPt > 25)


	triggerPass  = passIsoMu22 or passMu19Tau20

	if not triggerPass:
		continue



#Transverse mass
	metPt = tree.MET_pt
	metPhi = tree.MET_phi
	mt = math.sqrt(2 *muPt* metPt *(1-math.cos(muPhi - metPhi)))

# Fill histogram
	hMt.Fill(mt)
# To print the number of events processed for every 100000 events
	if i%100000 == 0 and i> 0:
		elapsed = time.time() - start
		rate = i/elapsed
		print("processed %d / %d (%.1f%%) at %.1f ev/s"%
		(i, numEntries, 100.0 * i /  numEntries, rate))


totalTime = time.time() - start
print("Loop finished in %.2f  seconds" % totalTime)


# Draw and save
c = ROOT.TCanvas("canvas_mt" , "Transeverse mass plot" , 800, 600)
hMt.Draw()
c.SaveAs("transeverseMassData.pdf")
