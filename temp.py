cd cernimport uproot
import awkward as ak
import numpy as np
import time

def deltaR(eta1, phi1, eta2, phi2):
    dphi = phi1 - phi2
    dphi = (dphi + np.pi) % (2 * np.pi) - np.pi  # wrap to [-pi, pi]
    deta = eta1 - eta2
    return np.hypot(deta, dphi)

# Open ROOT file
dataPath = "/home/neersj-shetty96/cernproj/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root"
f = uproot.open(dataPath)
tree = f["Events"]

# Load only relevant branches
branches = [
    "nMuon", "Muon_pt", "Muon_eta", "Muon_phi",
    "nTrigObj", "TrigObj_pt", "TrigObj_eta", "TrigObj_phi",
    "TrigObj_id", "TrigObj_filterBits",
    "MET_pt", "MET_phi",
    "HLT_IsoMu22", "HLT_IsoMu22_eta2p1",
    "HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
    "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
    "HLT_IsoTkMu22", "HLT_IsoTkMu22_eta2p1"
]
events = tree.arrays(branches, library="ak")

numEvents = len(events["nMuon"])
print("Loaded branches, number of events:", numEvents)

start = time.time()
mt_list = []

for i in range(numEvents):
    # progress print every 100,000 events
    if (i + 1) % 100000 == 0 or i == numEvents - 1:
        elapsed = time.time() - start
        rate = (i + 1) / elapsed
        print("Processed %d / %d events (%.1f%%) at %.1f ev/s" %
              (i + 1, numEvents, 100.0 * (i + 1) / numEvents, rate))

    nMuon = events["nMuon"][i]
    if nMuon == 0:
        continue

    # HLT paths fired
    HLTIsoMu22Fired = events["HLT_IsoMu22_eta2p1"][i]
    HLTIsoMu19Tau20Fired = (
        events["HLT_IsoMu19_eta2p1_LooseIsoPFTau20"][i] or
        events["HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1"][i]
    )

    if not (HLTIsoMu22Fired or HLTIsoMu19Tau20Fired):
        continue

    mu_pts = events["Muon_pt"][i]
    mu_etas = events["Muon_eta"][i]
    mu_phis = events["Muon_phi"][i]

    nTrigObj = events["nTrigObj"][i]
    trig_pts = events["TrigObj_pt"][i]
    trig_etas = events["TrigObj_eta"][i]
    trig_phis = events["TrigObj_phi"][i]
    trig_ids = events["TrigObj_id"][i]
    trig_filterBits = events["TrigObj_filterBits"][i]

    metPt = events["MET_pt"][i]
    metPhi = events["MET_phi"][i]

    for j in range(nMuon):
        muPt = mu_pts[j]
        muEta = mu_etas[j]
        muPhi = mu_phis[j]

        matched_muon = False
        matched_tau = False

        for k in range(nTrigObj):
            if abs(trig_ids[k]) != 13:
                continue
            if deltaR(muEta, muPhi, trig_etas[k], trig_phis[k]) >= 0.5:
                continue

            # IsoMu22 case
            if HLTIsoMu22Fired:
                if (trig_filterBits[k] & (1 << 1)) and (muPt > 23) and (abs(muEta) < 2.1):
                    matched_muon = True
                    break

            # Mu19Tau20 case
            elif HLTIsoMu19Tau20Fired:
                if (trig_filterBits[k] & (1 << 8)) and (20 < muPt < 23) and (abs(muEta) < 2.1):
                    matched_muon = True
                    # now check for tau match
                    for t in range(nTrigObj):
                        if abs(trig_ids[t]) != 15:
                            continue
                        if deltaR(muEta, muPhi, trig_etas[t], trig_phis[t]) >= 0.5:
                            continue
                        if trig_filterBits[t] & (1 << 8):
                            matched_tau = True
                            break
                    if matched_tau:
                        break

        if matched_muon:
            mt = np.sqrt(2 * muPt * metPt * (1 - np.cos(muPhi - metPhi)))
            mt_list.append(mt)

totalTime = time.time() - start
plot(mt_list)
print("Loop finished in %.2f seconds" % totalTime)
print("Number of muons passing all cuts:", len(mt_list))
