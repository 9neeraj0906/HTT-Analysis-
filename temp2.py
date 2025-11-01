import ROOT
import math
import time


def deltaR(eta1, phi1, eta2, phi2):
    dphi = phi1 - phi2
    if dphi > math.pi:
        dphi -= 2*math.pi
    elif dphi < -math.pi:
        dphi += 2*math.pi
    deta = eta1 - eta2
    return math.hypot(deta, dphi)

# Analysis function
def Analysis(file_path, hist_name="hMt", is_mc = False):
    f = ROOT.TFile.Open(file_path)
    tree = f.Get("Events")
    #Makes global directory the owner of the histogram created rather f
    ROOT.gDirectory.cd("PyROOT:/")
    # Histogram
    hMt = ROOT.TH1F(hist_name, "Transverse mass Data; MT(Mu,MET) [GeV]; Events", 100, 0.0, 150.0)

    numEvents = tree.GetEntries()
    print "Total events:", numEvents
    start = time.time()

    # Event loop
    for i in range(numEvents):
        tree.GetEntry(i)

        if tree.nMuon <= 0:
            continue

        MetFilters = [
        tree.Flag_goodVertices,
        tree.Flag_globalSuperTightHalo2016Filter,
        tree.Flag_HBHENoiseFilter,
        tree.Flag_HBHENoiseIsoFilter,
        tree.Flag_EcalDeadCellTriggerPrimitiveFilter,
        tree.Flag_BadPFMuonFilter,
        tree.Flag_eeBadScFilter,
        tree.Flag_ecalBadCalibFilter
        ]


        # Check if *all* filters are True; if not, skip the event
        if not all(MetFilters):
            continue

        # B-Jet veto
        deepcsv_loose = 0.2217   # 2016
        deepcsv_medium = 0.6321  # 2016

        veto_event = False
        for j in range(tree.nJet):
            jetPt = tree.Jet_pt[j]
            jetEta = tree.Jet_eta[j]
            btag = tree.Jet_btagDeepB[j]

            if jetPt > 25 and abs(jetEta) < 2.4:
                if btag > deepcsv_medium:  # Medium WP veto
                    veto_event = True
                    break  # no need to check more jets
                elif btag > deepcsv_loose:  # Loose WP veto
                    veto_event = True
                    break

        if veto_event:
            continue
        # Third lwpton veto
        vetoThirdLepton = False
        tau_candidates = []
        for p in range(tree.nTau):
            tauPt = tree.Tau_pt[p]
            tauEta = tree.Tau_eta[p]
            tauPhi = tree.Tau_phi[p]
            if tauPt > 20 and abs(tauEta) < 2.3:
                tau_candidates.append((tauEta, tauPhi))

        for e in range(tree.nElectron):
            electronPt = tree.Electron_pt[e]
            electronEta = tree.Electron_eta[e]
            electronPhi = tree.Electron_phi[e]
            electronDxy = tree.Electron_dxy[e]
            electronDz = tree.Electron_dz[e]

            if electronPt < 10 or abs(electronEta) > 2.5 or abs(electronDz) > 0.2 or abs(electronDxy) > 0.045:
                continue
            electronPass =(
                tree.Electron_mvaFall17V2noIso_WP90[e]
                and tree.Electron_convVeto[e]
            )
            missingInnerHits = tree.Electron_lostHits[e]
            if missingInnerHits > 1:
                continue
            if tree.Electron_pfRelIso03_all[e] >= 0.3*electronPt:
                continue

            if not electronPass:
                continue
            separated = all(deltaR(electronEta, electronPhi, tauEta, tauPhi) > 0.5
                            for tauEta, tauPhi in tau_candidates)
            if separated:
                vetoThirdLepton = True
                break

        if vetoThirdLepton:
            continue

        #Muon vriteria for all final states
        for j in range(tree.nMuon):
            muPt = tree.Muon_pt[j]
            muEta = tree.Muon_eta[j]
            muPhi = tree.Muon_phi[j]
            if muPt <= 10.0 or abs(muEta) >= 2.4:
                continue

            # Cut 2: Vertex Displacement (dZ < 0.2 cm and dxy < 0.045 cm)
            if abs(tree.Muon_dz[j]) >= 0.2 or abs(tree.Muon_dxy[j]) >= 0.045:
                 continue

            # Cut 3: Identification (Pass Medium ID)
            if not tree.Muon_mediumId[j]:
                 continue

            # Cut 4: Isolation (I_mu < 0.3 * P_T^mu), assuming pfRelIso04_all is I_mu / P_T^mu
            if tree.Muon_pfRelIso04_all[j] >= 0.3:
                 continue
        # Primaary vertex
        if tree.PV_npvsGood < 1:
            continue  # reject event

        # Check if the muon selected is from one of the HLT paths
        HLT_paths = [
            tree.HLT_IsoMu22,
            tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20,
            tree.HLT_IsoTkMu22,
            tree.HLT_IsoTkMu22_eta2p1,
            tree.HLT_IsoMu22_eta2p1,
            tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1,
        ]
        if not any(HLT_paths):
            continue

        HLTIsoMu22Fired = tree.HLT_IsoMu22_eta2p1
        HLTIsoMu19Tau20Fired = tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20 or tree.HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1

        if not (HLTIsoMu22Fired or HLTIsoMu19Tau20Fired):
            continue

        # Loop over offline muons
        for j in range(tree.nMuon):
            muPt = tree.Muon_pt[j]
            muEta = tree.Muon_eta[j]
            muPhi = tree.Muon_phi[j]


            if tree.Muon_dxy[j] >= 0.05 or \
               tree.Muon_dz[j] >= 0.2 or \
               tree.Muon_mediumId[j] == False or \
               tree.Muon_pfRelIso04_all[j] >= 0.15:
                continue

            muMatched = False

            # Loop over all trigger objects
            for k in range(tree.nTrigObj):
                trigId = tree.TrigObj_id[k]
                trigEta = tree.TrigObj_eta[k]
                trigPhi = tree.TrigObj_phi[k]
                filterBits = tree.TrigObj_filterBits[k]

                # Only muon trigger objects
                if abs(trigId) != 13:
                    continue

                # delta r match
                if deltaR(muEta, muPhi, trigEta, trigPhi) >= 0.5:
                    continue

                # Check trigger path
                if HLTIsoMu22Fired:
                    if (filterBits & (1 << 1)) and muPt > 23 and abs(muEta) < 2.1:
                        muMatched = True
                        break
                elif HLTIsoMu19Tau20Fired:
                    if (filterBits & (1 << 8)) and 20 < muPt < 23 and abs(muEta) < 2.1:
                        muMatched = True
                        break

            if not muMatched:
                continue

            # If Mu19Tau20 trigger, check taus
            tauMatched = True  # default for single-muon trigger
            if HLTIsoMu19Tau20Fired:
                tauMatched = False
                for p in range(tree.nTau):
                    tauPt = tree.Tau_pt[p]
                    tauEta = tree.Tau_eta[p]
                    tauPhi = tree.Tau_phi[p]

                    if tauPt < 20 or abs(tauEta) > 2.3:
                        continue
                    if tree.Tau_idDeepTau2017v2p1VSjet[p] <  16:
                        continue
                    if tree.Tau_idDeepTau2017v2p1VSe[p] < 1:
                        continue
                    if tree.Tau_idDeepTau2017v2p1VSmu[p] < 8:
                        continue
                    if tree.Tau_decayMode[p] != 0 and tree.Tau_decayMode[p] != 1 and tree.Tau_decayMode[p] != 10:
                        continue
                    # Loop over trigger objects again
                    for t in range(tree.nTrigObj):
                        trigIdTau = tree.TrigObj_id[t]
                        trigEtaTau = tree.TrigObj_eta[t]
                        trigPhiTau = tree.TrigObj_phi[t]
                        filterBitsTau = tree.TrigObj_filterBits[t]

                        if abs(trigIdTau) != 15:
                            continue


                        if deltaR(tauEta, tauPhi, trigEtaTau, trigPhiTau) >= 0.5:
                            continue

                        if filterBitsTau & (1 << 8):
                            tauMatched = True
                            break
                    if tauMatched:
                        break

                if not tauMatched:
                    continue

            # Compute MT
            metPt = tree.MET_pt
            metPhi = tree.MET_phi
            mt = math.sqrt(2 * muPt * metPt * (1 - math.cos(muPhi - metPhi)))
            event_weight = 1.0 # Default weight for data
            if is_mc:
                # If genWeight exists (it's an MC file), use its value
                event_weight = tree.genWeight
            hMt.Fill(mt, event_weight)

        # Progress
        if i % 100000 == 0 and i > 0:
            elapsed = time.time() - start
            rate = float(i)/elapsed
            print "Processed %d/%d (%.1f%%) at %.1f ev/s" % (i, numEvents, 100.0*i/numEvents, rate)

    totalTime = time.time() - start
    print "Loop finished in %.2f seconds" % totalTime
    f.Close()
    return hMt


dataPath = "/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root"
signalMCPath = "/data/smc/8544982A-CEC6-6B4D-B7F8-9B5AF213B725.root"
backgroundMCPath = "/data/bmc/0718C107-8960-6B44-B96A-C60D53D52A95.root"

hMtData = Analysis(dataPath, hist_name="hMtData", is_mc = False)
hMtSignal = Analysis(signalMCPath, hist_name="hMtSignal", is_mc = True)
hMtBackground = Analysis(backgroundMCPath, hist_name="hMtBackground", is_mc = True)


# Draw
c = ROOT.TCanvas("canvas_mt", "Transverse mass plot", 800, 600)
hMtData.SetLineColor(ROOT.kBlack)
hMtSignal.SetLineColor(ROOT.kRed)
hMtBackground.SetLineColor(ROOT.kBlue)

hMtData.Draw("HIST")
hMtSignal.Draw("HIST SAME")
hMtBackground.Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.65,0.7,0.88,0.88)
leg.AddEntry(hMtData,"Data","l")
leg.AddEntry(hMtSignal,"Signal MC","l")
leg.AddEntry(hMtBackground,"Background MC","l")
leg.Draw()

c.SaveAs("/pythonn/results/transverseMassData_simple1.pdf")
