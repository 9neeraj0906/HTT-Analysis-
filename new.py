import ROOT
import math
import time



def deltaR(eta1, phi1, eta2, phi2):
    dphi = phi1 - phi2
    deta = eta1 - eta2
    if dphi > math.pi:
        dphi -= 2*math.pi
    elif dphi < -math.pi:
        dphi += 2*math.pi
    return math.sqrt(abs(dphi) ** 2 + abs(deta) **2)



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
    # condition for atleast one muon
    if tree.nMuon <= 0:
        continue

    #Check if the muon selected is from one of the HLT paths or not
    HLT_paths = [
        tree.HLT_IsoMu22[i],
        tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20[i],
        tree.HLT_IsoTkMu22[i],
        tree.HLT_IsoTkMu22_eta2p1[i],
        tree.HLT_IsoMu22_eta2p1[i],
        tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20_singleL1[i],
    ]
    if not any(HLT_paths):
        continue

    # Furter Cuts w.r.t to the path fired
    HLTIsoMu22Fired = tree.HLT_IsoMu22_eta2p1[i]
    HLTIsoMu19Tau20Fired = tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20[i] or tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20_singleL1[i]

    #We have nMuon for every different entry, hence we have nMuon number of muonic variables/properties. For each nMuon there can be nTrigObj
    # Loop over all reconstructed muons
    for j in range(tree.nMuon):
        muPt  = tree.Muon_pt[j]
        muEta = tree.Muon_eta[j]
        muPhi = tree.Muon_phi[j]

        # Loop over all trigger objects for each muon
        for k in range(tree.nTrigObj):
            trigEta = tree.TrigObj_eta[k]
            trigPhi = tree.TrigObj_phi[k]
            trigId  = tree.TrigObj_id[k]
            filterBits = tree.TrigObj_filterBits[k]

            # Only keep muon trigger objects
            if abs(trigId) != 13:
                continue

            # ΔR matching between reco muon and trig muon
            if deltaR(muEta, muPhi, trigEta, trigPhi) >= 0.5:
                continue

            # Case 1: IsoMu22 fired
            if HLTIsoMu22Fired:
                isoFilterBit = 1  # TrkIsoVVL
                if filterBits & (1 << isoFilterBit):
                    if muPt > 23 and abs(muEta) < 2.1:
                        pass  # matched to IsoMu22

            # Case 2: Mu19Tau20 fired
            elif HLTIsoMu19Tau20Fired:
                    lastFilterBitMuon = 8  # IsoTkMu (last filter for muon)
                    if filterBits & (1 << lastFilterBitMuon):
                        if 20 < muPt < 23 and abs(muEta) < 2.1:
                            # Now look for matching tau trigger object
                            for t in range(tree.nTrigObj):
                                if abs(tree.TrigObj_id[t]) != 15:
                                    continue  # skip non-tau objects

                                tauEta = tree.TrigObj_eta[t]
                                tauPhi = tree.TrigObj_phi[t]
                                tauFilterBits = tree.TrigObj_filterBits[t]
                                # ΔR match between muon and tau trig objs
                                if deltaR(muEta, muPhi, tauEta, tauPhi) >= 0.5:
                                    continue
                                # Tau should pass its last filter ( bit 8)
                                if tauFilterBits & (1 << 8):  # adjust bit number if needed
                                    pass  #  Both muon and tau matched successfully
