# -*- coding: utf-8 -*-
from __future__ import print_function

import ROOT
import uproot
import awkward as ak
import numpy as np
import time

# ================================================================
# USE THESE 5 FILES FIRST
# ================================================================
FILES = [
    "/data/bmc/D20B2931-4C27-FF4E-B640-3404DFD65150.root",
    "/data/bmc/AE5A627E-C8D7-974F-9698-C136D2BFFF1D.root",
    "/data/bmc/5544E66D-B271-1649-A85B-C08F01375C57.root",
    "/data/bmc/D6094FF2-1AAE-784B-9991-AA8441B9FB45.root",
    "/data/bmc/E215B4B9-11EA-D84F-A6D7-9FE57515B2CE.root"
]

# ================================================================
# Delta R
# ================================================================
def dR(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = np.arctan2(np.sin(phi1 - phi2), np.cos(phi1 - phi2))
    return np.sqrt(deta * deta + dphi * dphi)

# ================================================================
# Generator-level Z decay classifier.
# We start from every Z and walk down its daughter chain until we
# reach a lepton/neutrino. Tau is stopped at the tau itself, so a
# Z->tautau event is not mistaken for Z->mumu from tau decay.
# ================================================================
def classify_z_decay(pdg, mother):
    pdg = np.asarray(pdg)
    mother = np.asarray(mother)

    z_indices = np.where(np.abs(pdg) == 23)[0]

    for z in z_indices:
        found = []
        stack = []

        children = np.where(mother == z)[0]
        for child in children:
            stack.append(int(child))

        visited = set()

        while len(stack) > 0:
            p = stack.pop()

            if p in visited:
                continue
            visited.add(p)

            apdg = abs(int(pdg[p]))

            # Stop at the primary Z decay particle.
            if apdg in (11, 12, 13, 14, 15, 16):
                found.append(apdg)
                continue

            daughters = np.where(mother == p)[0]
            for child in daughters:
                stack.append(int(child))

        # Need two primary decay products of the same species.
        if found.count(15) >= 2:
            return "Z->tautau"
        if found.count(13) >= 2:
            return "Z->mumu"
        if found.count(11) >= 2:
            return "Z->ee"
        if (found.count(12) + found.count(14) + found.count(16)) >= 2:
            return "Z->nunu"

    return "other/unknown"

# ================================================================
# Generator tau matching.
# ================================================================
def has_gen_tau_match(gen_pt, gen_eta, gen_phi, reco_eta, reco_phi):
    gen_pt = np.asarray(gen_pt)
    gen_eta = np.asarray(gen_eta)
    gen_phi = np.asarray(gen_phi)

    mask = (np.abs(np.asarray(gen_pt)) >= 0.0) & (np.abs(np.asarray(gen_eta)) < 10.0)
    gen_eta = gen_eta[mask]
    gen_phi = gen_phi[mask]

    if len(gen_eta) == 0:
        return False

    dr = dR(gen_eta, gen_phi, reco_eta, reco_phi)
    return bool(np.min(dr) < 0.5)

# ================================================================
# Process one file.
# ================================================================
def process_file(file_path):
    file_start = time.time()

    print("\n============================================================")
    print("Processing:")
    print(file_path)
    print("============================================================")

    tree = uproot.open(file_path)["Events"]

    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
        "Muon_dz", "Muon_dxy", "Muon_mediumId", "Muon_pfRelIso04_all",
        "MET_pt", "MET_phi",
        "Jet_pt", "Jet_eta", "Jet_btagDeepB",
        "Electron_pt", "Electron_eta", "Electron_phi",
        "Electron_dxy", "Electron_dz",
        "Electron_mvaFall17V2noIso_WP90", "Electron_convVeto",
        "Electron_lostHits", "Electron_pfRelIso03_all",
        "Tau_pt", "Tau_eta", "Tau_phi", "Tau_charge", "Tau_dz",
        "Tau_idDeepTau2017v2p1VSe", "Tau_idDeepTau2017v2p1VSmu",
        "Tau_idDeepTau2017v2p1VSjet", "Tau_decayMode",
        "Tau_rawDeepTau2017v2p1VSjet",
        "PV_npvsGood",
        "Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter",
        "Flag_eeBadScFilter", "Flag_ecalBadCalibFilter",
        "HLT_IsoMu22_eta2p1",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
        "TrigObj_id", "TrigObj_eta", "TrigObj_phi", "TrigObj_filterBits"
    ]

    gen_branches = [
        "GenPart_pdgId", "GenPart_genPartIdxMother",
        "GenPart_pt", "GenPart_eta", "GenPart_phi"
    ]

    available = set(tree.keys())
    missing = [x for x in branches + gen_branches if x not in available]
    if len(missing) > 0:
        print("ERROR: missing branches:")
        for x in missing:
            print("  " + x)
        return None

    # Read reco arrays and generator arrays from the ORIGINAL event order.
    br = tree.arrays(branches, namedecode="utf-8")
    gen = tree.arrays(gen_branches, namedecode="utf-8")

    n_original = len(br["MET_pt"])
    event_index = np.arange(n_original, dtype=np.int64)

    print("Initial events:", n_original)

    def event_filter(mask):
        nonlocal_event = None
        # Python 2 has no nonlocal; this helper directly modifies the
        # dictionaries and event_index through explicit return values below.
        return mask

    # ------------------------------------------------------------
    # 1. MET filters
    # ------------------------------------------------------------
    met_mask = (
        br["Flag_goodVertices"] & br["Flag_globalSuperTightHalo2016Filter"] &
        br["Flag_HBHENoiseFilter"] & br["Flag_HBHENoiseIsoFilter"] &
        br["Flag_EcalDeadCellTriggerPrimitiveFilter"] &
        br["Flag_BadPFMuonFilter"] & br["Flag_eeBadScFilter"] &
        br["Flag_ecalBadCalibFilter"]
    )
    for key in br:
        br[key] = br[key][met_mask]
    event_index = event_index[met_mask]
    print("After MET filters:", len(br["MET_pt"]))

    # ------------------------------------------------------------
    # 2. b-jet veto
    # ------------------------------------------------------------
    jets_good = (br["Jet_pt"] > 20.0) & (np.abs(br["Jet_eta"]) < 2.4)
    bjet_veto = (br["Jet_btagDeepB"][jets_good] > 0.6321).any()
    keep = ~bjet_veto
    for key in br:
        br[key] = br[key][keep]
    event_index = event_index[keep]
    print("After b-jet veto:", len(br["MET_pt"]))

    # ------------------------------------------------------------
    # 3. Current electron veto used in your analysis
    # ------------------------------------------------------------
    tau_veto = (br["Tau_pt"] > 20.0) & (np.abs(br["Tau_eta"]) < 2.3)

    ele_good = (
        (br["Electron_pt"] >= 10.0) &
        (np.abs(br["Electron_eta"]) <= 2.5) &
        (np.abs(br["Electron_dz"]) <= 0.2) &
        (np.abs(br["Electron_dxy"]) <= 0.045) &
        (br["Electron_mvaFall17V2noIso_WP90"] == True) &
        (br["Electron_convVeto"] == True) &
        (br["Electron_lostHits"] <= 1) &
        (br["Electron_pfRelIso03_all"] < 0.15)
    )

    ele_sel = ak.JaggedArray.zip(
        eta=br["Electron_eta"][ele_good],
        phi=br["Electron_phi"][ele_good]
    )
    tau_sel = ak.JaggedArray.zip(
        eta=br["Tau_eta"][tau_veto],
        phi=br["Tau_phi"][tau_veto]
    )

    pairs_et = ele_sel.cross(tau_sel, nested=True)
    e_et = pairs_et.i0
    t_et = pairs_et.i1
    ele_separated = (dR(e_et["eta"], e_et["phi"], t_et["eta"], t_et["phi"]) > 0.5).all()
    veto_ele = ele_separated.any()

    keep = ~veto_ele
    for key in br:
        br[key] = br[key][keep]
    event_index = event_index[keep]
    print("After electron veto:", len(br["MET_pt"]))

    # ------------------------------------------------------------
    # 4. Muon selection
    # ------------------------------------------------------------
    mu_strict = (
        (br["Muon_pt"] > 10.0) &
        (np.abs(br["Muon_eta"]) < 2.4) &
        (np.abs(br["Muon_dz"]) < 0.2) &
        (np.abs(br["Muon_dxy"]) < 0.045) &
        (br["Muon_mediumId"] == True) &
        (br["Muon_pfRelIso04_all"] < 0.15) &
        (br["PV_npvsGood"] >= 1)
    )

    # ------------------------------------------------------------
    # 5. Trigger decision
    # ------------------------------------------------------------
    hlt_mu22 = br["HLT_IsoMu22_eta2p1"]
    hlt_mu19 = (
        br["HLT_IsoMu19_eta2p1_LooseIsoPFTau20"] |
        br["HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1"]
    )
    event_hlt_ok = hlt_mu22 | hlt_mu19

    # ------------------------------------------------------------
    # 6. Muon trigger matching
    # ------------------------------------------------------------
    mu = ak.JaggedArray.zip(
        pt=br["Muon_pt"],
        eta=br["Muon_eta"],
        phi=br["Muon_phi"],
        charge=br["Muon_charge"]
    )

    trig_mu_mask = np.abs(br["TrigObj_id"]) == 13
    trig_mu = ak.JaggedArray.zip(
        eta=br["TrigObj_eta"][trig_mu_mask],
        phi=br["TrigObj_phi"][trig_mu_mask],
        filterBits=br["TrigObj_filterBits"][trig_mu_mask]
    )

    pairs_mu = mu.cross(trig_mu, nested=True)
    m_mu = pairs_mu.i0
    t_mu = pairs_mu.i1

    dr_mu = dR(m_mu["eta"], m_mu["phi"], t_mu["eta"], t_mu["phi"])
    filter_mu22 = (t_mu["filterBits"] & (1 << 1)) != 0
    filter_mu19 = (t_mu["filterBits"] & (1 << 8)) != 0

    mu22_kin = (m_mu["pt"] > 23.0) & (np.abs(m_mu["eta"]) < 2.1)
    mu19_kin = (m_mu["pt"] > 20.0) & (m_mu["pt"] < 23.0) & (np.abs(m_mu["eta"]) < 2.1)

    match_mu22 = dr_mu < 0.5
    match_mu22 = match_mu22 & filter_mu22 & mu22_kin & hlt_mu22

    match_mu19 = dr_mu < 0.5
    match_mu19 = match_mu19 & filter_mu19 & mu19_kin & (~hlt_mu22) & hlt_mu19

    mu22_matched = match_mu22.any()
    mu19_matched = match_mu19.any()
    mu_matched = mu22_matched | mu19_matched

    # ------------------------------------------------------------
    # 7. Tau base selection
    # ------------------------------------------------------------
    tau_base = (
        (np.abs(br["Tau_eta"]) < 2.3) &
        (np.abs(br["Tau_dz"]) < 0.2) &
        (br["Tau_idDeepTau2017v2p1VSe"] >= 1) &
        (br["Tau_idDeepTau2017v2p1VSmu"] >= 8) &
        (br["Tau_idDeepTau2017v2p1VSjet"] >= 16) &
        ((br["Tau_decayMode"] == 0) |
         (br["Tau_decayMode"] == 1) |
         (br["Tau_decayMode"] == 10))
    )

    tau = ak.JaggedArray.zip(
        pt=br["Tau_pt"][tau_base],
        eta=br["Tau_eta"][tau_base],
        phi=br["Tau_phi"][tau_base],
        charge=br["Tau_charge"][tau_base],
        djet=br["Tau_rawDeepTau2017v2p1VSjet"][tau_base]
    )

    # ------------------------------------------------------------
    # 8. Tau trigger matching
    # ------------------------------------------------------------
    trig_tau_mask = np.abs(br["TrigObj_id"]) == 15
    trig_tau = ak.JaggedArray.zip(
        eta=br["TrigObj_eta"][trig_tau_mask],
        phi=br["TrigObj_phi"][trig_tau_mask],
        filterBits=br["TrigObj_filterBits"][trig_tau_mask]
    )

    pairs_tau = tau.cross(trig_tau, nested=True)
    ta = pairs_tau.i0
    tt = pairs_tau.i1
    dr_tau = dR(ta["eta"], ta["phi"], tt["eta"], tt["phi"])
    filter_tau20 = (tt["filterBits"] & (1 << 8)) != 0
    tau_trig_matched = ((dr_tau < 0.5) & filter_tau20).any()

    # ------------------------------------------------------------
    # 9. Attach validity flags to muons and taus
    # ------------------------------------------------------------
    mu_selected = ak.JaggedArray.zip(
        pt=br["Muon_pt"],
        eta=br["Muon_eta"],
        phi=br["Muon_phi"],
        charge=br["Muon_charge"],
        valid=(mu_strict & mu_matched & event_hlt_ok),
        mu22=mu22_matched,
        mu19tau=mu19_matched
    )

    tau_selected = ak.JaggedArray.zip(
        pt=tau["pt"],
        eta=tau["eta"],
        phi=tau["phi"],
        charge=tau["charge"],
        djet=tau["djet"],
        trigmatched=tau_trig_matched
    )

    # ------------------------------------------------------------
    # 10. mu x tau pair selection
    # ------------------------------------------------------------
    pairs = mu_selected.cross(tau_selected)
    pair_mu = pairs.i0
    pair_tau = pairs.i1

    pair_dr_ok = dR(
        pair_mu["eta"], pair_mu["phi"],
        pair_tau["eta"], pair_tau["phi"]
    ) > 0.5

    pair_os = (pair_mu["charge"] * pair_tau["charge"]) < 0

    tau_pt_ok = (
        (pair_mu["mu19tau"] & (pair_tau["pt"] > 32.0)) |
        (pair_mu["mu22"] & (pair_tau["pt"] > 30.0))
    )

    tau_trigger_ok = (
        pair_mu["mu22"] |
        (pair_mu["mu19tau"] & pair_tau["trigmatched"])
    )

    valid_pair = (
        pair_mu["valid"] &
        pair_dr_ok &
        pair_os &
        tau_pt_ok &
        tau_trigger_ok
    )

    valid_mu = pair_mu[valid_pair]
    valid_tau = pair_tau[valid_pair]

    best_idx = valid_tau["djet"].argmax()
    selected_mu_pt = valid_mu["pt"][best_idx]
    selected_mu_eta = valid_mu["eta"][best_idx]
    selected_mu_phi = valid_mu["phi"][best_idx]
    selected_tau_pt = valid_tau["pt"][best_idx]
    selected_tau_eta = valid_tau["eta"][best_idx]
    selected_tau_phi = valid_tau["phi"][best_idx]

    has_pair = selected_tau_pt.counts > 0

    selected_mu_pt = selected_mu_pt[has_pair]
    selected_mu_eta = selected_mu_eta[has_pair]
    selected_mu_phi = selected_mu_phi[has_pair]
    selected_tau_pt = selected_tau_pt[has_pair]
    selected_tau_eta = selected_tau_eta[has_pair]
    selected_tau_phi = selected_tau_phi[has_pair]

    final_original_indices = event_index[has_pair]

    # ------------------------------------------------------------
    # 11. Visible mass
    # ------------------------------------------------------------
    deta = selected_mu_eta - selected_tau_eta
    dphi = np.arctan2(
        np.sin(selected_mu_phi - selected_tau_phi),
        np.cos(selected_mu_phi - selected_tau_phi)
    )

    mvis = np.sqrt(
        2.0 * selected_mu_pt * selected_tau_pt *
        (np.cosh(deta) - np.cos(dphi))
    )

    mvis_flat = np.asarray(mvis.flatten())
    final_indices = np.asarray(final_original_indices)

    print("Final mu-tau events:", len(mvis_flat))
    print("Final Mvis 60--80 GeV:", int(np.sum((mvis_flat >= 60.0) & (mvis_flat < 80.0))))

    # ------------------------------------------------------------
    # 12. Generator composition
    # ------------------------------------------------------------
    gen_pdg = gen["GenPart_pdgId"]
    gen_mother = gen["GenPart_genPartIdxMother"]
    gen_pt = gen["GenPart_pt"]
    gen_eta = gen["GenPart_eta"]
    gen_phi = gen["GenPart_phi"]

    counts = {
        "Z->tautau": 0,
        "Z->mumu": 0,
        "Z->ee": 0,
        "Z->nunu": 0,
        "other/unknown": 0
    }

    counts_60_80 = {
        "Z->tautau": 0,
        "Z->mumu": 0,
        "Z->ee": 0,
        "Z->nunu": 0,
        "other/unknown": 0
    }

    tau_match_total = 0
    tau_match_60_80 = 0

    for j in range(len(final_indices)):
        idx = int(final_indices[j])
        category = classify_z_decay(gen_pdg[idx], gen_mother[idx])
        counts[category] += 1

        in_60_80 = (mvis_flat[j] >= 60.0 and mvis_flat[j] < 80.0)
        if in_60_80:
            counts_60_80[category] += 1

        matched = has_gen_tau_match(
            gen_pt[idx],
            gen_eta[idx],
            gen_phi[idx],
            selected_tau_eta.flatten()[j],
            selected_tau_phi.flatten()[j]
        )

        if matched:
            tau_match_total += 1
            if in_60_80:
                tau_match_60_80 += 1

    print("\nGENERATOR COMPOSITION")
    print("---------------------")
    total = len(final_indices)
    for key in ["Z->tautau", "Z->mumu", "Z->ee", "Z->nunu", "other/unknown"]:
        pct = 100.0 * counts[key] / total if total else 0.0
        print("%-15s %6d  %7.2f %%" % (key, counts[key], pct))

    n6080 = int(np.sum((mvis_flat >= 60.0) & (mvis_flat < 80.0)))
    print("\nGENERATOR COMPOSITION: Mvis 60--80 GeV")
    print("--------------------------------------")
    for key in ["Z->tautau", "Z->mumu", "Z->ee", "Z->nunu", "other/unknown"]:
        pct = 100.0 * counts_60_80[key] / n6080 if n6080 else 0.0
        print("%-15s %6d  %7.2f %%" % (key, counts_60_80[key], pct))

    tau_pct = 100.0 * tau_match_total / total if total else 0.0
    tau_pct_6080 = 100.0 * tau_match_60_80 / n6080 if n6080 else 0.0

    print("\nRECONSTRUCTED TAU MATCH")
    print("-----------------------")
    print("All final events: %d / %d = %.2f %%" % (tau_match_total, total, tau_pct))
    print("Mvis 60--80 GeV: %d / %d = %.2f %%" % (tau_match_60_80, n6080, tau_pct_6080))

    print("\nFile completed in %.1fs" % (time.time() - file_start))

    return counts, counts_60_80, tau_match_total, tau_match_60_80

# ================================================================
# Main
# ================================================================
if __name__ == "__main__":
    total_counts = {
        "Z->tautau": 0,
        "Z->mumu": 0,
        "Z->ee": 0,
        "Z->nunu": 0,
        "other/unknown": 0
    }

    total_counts_60_80 = {
        "Z->tautau": 0,
        "Z->mumu": 0,
        "Z->ee": 0,
        "Z->nunu": 0,
        "other/unknown": 0
    }

    total_tau_match = 0
    total_tau_match_60_80 = 0

    start = time.time()

    for i in range(len(FILES)):
        print("\n################ FILE %d / %d ################" % (i + 1, len(FILES)))
        result = process_file(FILES[i])

        if result is None:
            continue

        counts, counts_60_80, tau_match, tau_match_60_80 = result

        for key in total_counts:
            total_counts[key] += counts[key]
            total_counts_60_80[key] += counts_60_80[key]

        total_tau_match += tau_match
        total_tau_match_60_80 += tau_match_60_80

    total = sum(total_counts.values())
    total_6080 = sum(total_counts_60_80.values())

    print("\n")
    print("============================================================")
    print("COMBINED GENERATOR COMPOSITION")
    print("============================================================")
    for key in ["Z->tautau", "Z->mumu", "Z->ee", "Z->nunu", "other/unknown"]:
        pct = 100.0 * total_counts[key] / total if total else 0.0
        print("%-15s %7d  %7.2f %%" % (key, total_counts[key], pct))

    print("\n============================================================")
    print("COMBINED Mvis 60--80 GeV")
    print("============================================================")
    print("Total events:", total_6080)
    for key in ["Z->tautau", "Z->mumu", "Z->ee", "Z->nunu", "other/unknown"]:
        pct = 100.0 * total_counts_60_80[key] / total_6080 if total_6080 else 0.0
        print("%-15s %7d  %7.2f %%" % (key, total_counts_60_80[key], pct))

    print("\n============================================================")
    print("RECONSTRUCTED TAU MATCH")
    print("============================================================")
    tau_pct = 100.0 * total_tau_match / total if total else 0.0
    tau_pct_6080 = 100.0 * total_tau_match_60_80 / total_6080 if total_6080 else 0.0
    print("All final events: %d / %d = %.2f %%" % (total_tau_match, total, tau_pct))
    print("Mvis 60--80 GeV: %d / %d = %.2f %%" % (total_tau_match_60_80, total_6080, tau_pct_6080))

    print("\nTotal runtime: %.1fs" % (time.time() - start))
    print("============================================================")
