import ROOT

# --- Open your file ---
f = ROOT.TFile.Open("/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root")
tree = f.Get("Events")

# --- Keep track of tau filterBits values ---
seen_tau_bits = set()

# --- Loop through a few events ---
n = min(tree.GetEntries(), 5000)  # check up to 5000 events
for i in range(n):
    tree.GetEntry(i)

    # Check if Mu19Tau20 path fired
    if not getattr(tree, "HLT_IsoMu19_eta2p1_LooseIsoPFTau20", False):
        continue

    # Loop over all trigger objects in this event
    for j in range(tree.nTrigObj):
        trigId = int(tree.TrigObj_id[j])
        filterBits = int(tree.TrigObj_filterBits[j])

        # Tau trigger objects have id = 15
        if abs(trigId) == 15:
            seen_tau_bits.add(filterBits)

# --- Show unique filterBits values for taus ---
print("Unique filterBits values seen for tau trigger objects:")
for v in sorted(seen_tau_bits):
    print("%5d -> binary: %s" % (v, bin(v)))
