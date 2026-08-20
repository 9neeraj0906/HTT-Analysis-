# -*- coding: utf-8 -*-
import ROOT
import uproot
import awkward as ak
import numpy as np
import math
import time

#=====================================================================
# VECTORIZED VERSION -- FIXED
#
# What was wrong before:
#   1. Event-level cuts (met filter, b-jet veto, 3rd-lepton veto) collapsed
#      muon info using ak.any(), which is fine -- those really are event-level.
#   2. But the FINAL muon selection used ak.firsts(Muon_pt[muon_strict]),
#      i.e. "the first muon in the event that passes strict ID", completely
#      ignoring which muon (if any) was actually trigger-matched. The loop
#      version, by contrast, fills once for EVERY muon in the event that
#      passes strict ID *and* is itself the one matched to a trigger object
#      *and* survives the tau-matching requirement.
#   3. Because of that collapse, events with >=2 qualifying muons only ever
#      got ONE fill in the vectorized version (vs. possibly N fills in the
#      loop), and even single-muon events could get filled with the WRONG
#      muon's pt/phi if muon ordering in the branch didn't happen to put the
#      matched muon first.
#
# Fix: keep everything at muon-level (event, muon) shape all the way through
# to the fill step. Only flatten at the very end, once we know exactly which
# (event, muon) pairs are meant to produce a histogram entry -- this lets an
# event contribute zero, one, or multiple entries, exactly like the loop.
#
# NORMALIZATION FIX (added after the background scale investigation):
#   genWeight for this sample is NOT simply +-1 -- it's a large, varying
#   MadGraph-style generator weight. The correct normalization divides by
#   the SUM of genWeight across the ENTIRE physical sample (all files),
#   not the raw per-file event count, and NOT any single file's own local
#   genEventSumw (each file's Runs-tree genEventSumw is local to that file
#   only -- confirmed by direct inspection: it matches that file's own
#   sum(genWeight) closely). xsec is also no longer hardcoded, since a
#   signal sample and a background sample must never share one constant.
# =====================================================================

def get_total_gen_event_sumw(file_list):
    """
    Sums genEventSumw across every file in file_list. This is the CMS-
    standard normalization denominator for MC samples with non-unit
    genWeight values (e.g. MadGraph/aMC@NLO), and it MUST be the total
    across the whole physical sample -- not any single file's local
    value -- since each file's Runs tree only covers events generated
    for that file. Using a per-file local value as its own denominator
    would make every file normalize itself as if it were the entire
    sample, inflating the summed yield by roughly (number of files)x.

    Call this ONCE per sample (signal, or each background process)
    before looping over files, then pass the result into every
    Analysis_up(..., genEventSumw=total) call for that sample.
    """
    total = 0.0
    for fp in file_list:
        runs = uproot.open(fp)["Runs"]
        total += runs.array("genEventSumw").sum()
    return total


def Analysis_up(filePath, hist_name="hMt_vec", is_mc=False, xsec=6077.22, luminosity=35920.0, genEventSumw=None):
    """
    uproot3 + awkward0 (Python 2) implementation. Same logic as the
    for-loop version, verified cut-by-cut, ported to the legacy
    JaggedArray API instead of modern ak.Array/ak.zip/ak.cartesian.

    Key API differences from the modern-awkward version:
      - tree.arrays(branches, namedecode="utf-8") returns a plain dict,
        not an ak.Array with attribute access -- so we index with [] and
        stash results in a dict called `br`.
      - JaggedArray.cross(other, nested=True) replaces ak.cartesian; the
        two sides come back as .i0 / .i1 (not named fields).
      - JaggedArray.any() / .all() reduce the innermost (last) axis by
        default -- equivalent to ak.any(..., axis=-1).
      - Plain numpy per-event arrays broadcast automatically against
        JaggedArrays with ordinary operators (+, *, &, etc.) -- no
        ak.broadcast_arrays equivalent needed.
      - .all() / .any() on an empty sublist is vacuously True / False,
        same as Python's built-in all([]) / any([]), matching the loop.

    For is_mc=True, `genEventSumw` MUST be supplied -- the total sum of
    genWeight across the ENTIRE physical sample (all files), obtained
    via get_total_gen_event_sumw(). This is required because genWeight
    is not simply +-1 for this generator; the correct normalization
    divides by the sum of generator weights, not the raw event count,
    and that sum must be the dataset-wide total, not any single file's
    local contribution.
    """
    tree = uproot.open(filePath)["Events"]

    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi",
        "Muon_dz", "Muon_dxy", "Muon_mediumId",
        "Muon_pfRelIso04_all",

        "MET_pt", "MET_phi",

        "Jet_pt", "Jet_eta", "Jet_btagDeepB",

        "Electron_pt", "Electron_eta", "Electron_phi",
        "Electron_dxy", "Electron_dz",
        "Electron_mvaFall17V2noIso_WP90", "Electron_convVeto",
        "Electron_lostHits", "Electron_pfRelIso03_all",

        "Tau_pt", "Tau_eta", "Tau_phi",
        "Tau_idDeepTau2017v2p1VSe",
        "Tau_idDeepTau2017v2p1VSmu",
        "Tau_idDeepTau2017v2p1VSjet",
        "Tau_decayMode",

        "PV_npvsGood",

        "Flag_goodVertices",
        "Flag_globalSuperTightHalo2016Filter",
        "Flag_HBHENoiseFilter",
        "Flag_HBHENoiseIsoFilter",
        "Flag_EcalDeadCellTriggerPrimitiveFilter",
        "Flag_BadPFMuonFilter",
        "Flag_eeBadScFilter",
        "Flag_ecalBadCalibFilter",

        "HLT_IsoMu22_eta2p1",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
        "HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",

        "TrigObj_id", "TrigObj_eta", "TrigObj_phi",
        "TrigObj_filterBits",
    ]

    if is_mc:
        branches.append("genWeight")

    br = tree.arrays(branches, namedecode="utf-8")

    # numEvents is still tracked for informational printing, but is NO
    # LONGER used to normalize MC weights (see genEventSumw parameter).
    numEvents = len(br["PV_npvsGood"])
    print "Total events:", numEvents

    if is_mc and genEventSumw is None:
        raise ValueError(
            "Analysis_up called with is_mc=True but no genEventSumw was "
            "supplied. Compute it once via get_total_gen_event_sumw(file_list) "
            "for the WHOLE sample (all files), then pass it in explicitly -- "
            "using a per-file local value here would silently reproduce the "
            "normalization bug this parameter exists to prevent."
        )

    def dR(eta1, phi1, eta2, phi2):
        deta = eta1 - eta2
        dphi = np.arctan2(np.sin(phi1 - phi2), np.cos(phi1 - phi2))
        return np.sqrt(deta**2 + dphi**2)

    def event_filter(mask):
        """Apply an event-level boolean mask to every branch in br."""
        for key in br:
            br[key] = br[key][mask]

    # =================================================================
    # 1. MET Filters (event-level -- safe to filter directly)
    # =================================================================
    met_mask = (
        br["Flag_goodVertices"] &
        br["Flag_globalSuperTightHalo2016Filter"] &
        br["Flag_HBHENoiseFilter"] &
        br["Flag_HBHENoiseIsoFilter"] &
        br["Flag_EcalDeadCellTriggerPrimitiveFilter"] &
        br["Flag_BadPFMuonFilter"] &
        br["Flag_eeBadScFilter"] &
        br["Flag_ecalBadCalibFilter"]
    )
    event_filter(met_mask)

    # =================================================================
    # 2. B-Jet Veto (event-level)
    # =================================================================
    deepcsv_medium = 0.6321
    jets_good = (br["Jet_pt"] > 25) & (abs(br["Jet_eta"]) < 2.4)
    bjetVeto = (br["Jet_btagDeepB"][jets_good] > deepcsv_medium).any()
    event_filter(~bjetVeto)

    # =================================================================
    # 3. Third Lepton Veto (Electrons) (event-level)
    # =================================================================
    tau_good = (br["Tau_pt"] > 20) & (abs(br["Tau_eta"]) < 2.3)

    ele_good = (
        (br["Electron_pt"] >= 10) &
        (abs(br["Electron_eta"]) <= 2.5) &
        (abs(br["Electron_dz"]) <= 0.2) &
        (abs(br["Electron_dxy"]) <= 0.045) &
        (br["Electron_mvaFall17V2noIso_WP90"] == True) &
        (br["Electron_convVeto"] == True) &
        (br["Electron_lostHits"] <= 1) &
        (br["Electron_pfRelIso03_all"] < 0.3 * br["Electron_pt"])
    )

    ele_sel = ak.JaggedArray.zip(eta=br["Electron_eta"][ele_good],
                                       phi=br["Electron_phi"][ele_good])
    tau_sel = ak.JaggedArray.zip(eta=br["Tau_eta"][tau_good],
                                       phi=br["Tau_phi"][tau_good])

    pairs_et = ele_sel.cross(tau_sel, nested=True)
    e_et, t_et = pairs_et.i0, pairs_et.i1
    dr_et = dR(e_et["eta"], e_et["phi"], t_et["eta"], t_et["phi"])

    ele_separated = (dr_et > 0.5).all()   # reduces tau axis -> per electron
    veto_ele      = ele_separated.any()   # reduces electron axis -> per event
    event_filter(~veto_ele)

    # From here on we do NOT do any more event-level filtering that
    # discards which muon qualified. Everything stays at (event, muon) shape.

    # =================================================================
    # 4. Muon strict selection mask (preselection cuts are a subset of
    #    these, so this single mask reproduces both loop stages combined)
    # =================================================================
    def get_muon_strict(b):
        return (
            (b["Muon_pt"] > 10) &
            (abs(b["Muon_eta"]) < 2.4) &
            (abs(b["Muon_dz"]) < 0.2) &
            (abs(b["Muon_dxy"]) < 0.045) &
            (b["Muon_mediumId"] == True) &
            (b["Muon_pfRelIso04_all"] < 0.3) &
            (b["PV_npvsGood"] >= 1) &
            (abs(b["Muon_dxy"]) < 0.05) &      # redundant given |dxy|<0.045 above, kept for fidelity
            (b["Muon_pfRelIso04_all"] < 0.15)
        )

    muon_strict = get_muon_strict(br)   # shape: (event, muon)

    # =================================================================
    # 5. HLT paths (event-level flags, but NOT used to drop events yet --
    #    only to gate the muon-level fill mask, same as the loop does)
    # =================================================================
    # IMPORTANT: the loop's HLTIsoMu22Fired is ONLY HLT_IsoMu22_eta2p1 --
    # not the OR of all four IsoMu22-like paths. The broader 6-path
    # "any HLT fired" pre-check earlier in the loop is a strict superset
    # of (HLTIsoMu22Fired or HLTIsoMu19Tau20Fired) and therefore never
    # rejects anything this narrower check wouldn't already reject, so it
    # doesn't need its own representation here.
    HLTIsoMu22Fired = br["HLT_IsoMu22_eta2p1"]
    HLTIsoMu19Tau20Fired = (
        br["HLT_IsoMu19_eta2p1_LooseIsoPFTau20"] |
        br["HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1"]
    )
    event_hlt_ok = HLTIsoMu22Fired | HLTIsoMu19Tau20Fired   # flat numpy, shape: (event,)

    # =================================================================
    # 6. Muon trigger matching -- computed for ALL muons (not just the
    #    strict-passing ones), then ANDed with muon_strict at the end.
    #    This keeps per-muon identity intact.
    # =================================================================
    mu = ak.JaggedArray.zip(pt=br["Muon_pt"], eta=br["Muon_eta"], phi=br["Muon_phi"])

    trig_mu_mask = abs(br["TrigObj_id"]) == 13
    trig_mu = ak.JaggedArray.zip(
        eta=br["TrigObj_eta"][trig_mu_mask],
        phi=br["TrigObj_phi"][trig_mu_mask],
        filterBits=br["TrigObj_filterBits"][trig_mu_mask],
    )

    pairs_mu = mu.cross(trig_mu, nested=True)
    m_mu, t_mu = pairs_mu.i0, pairs_mu.i1
    dr_mu = dR(m_mu["eta"], m_mu["phi"], t_mu["eta"], t_mu["phi"])

    dr_ok          = dr_mu < 0.5
    filter_mu22    = (t_mu["filterBits"] & (1 << 1)) != 0
    filter_mu19tau = (t_mu["filterBits"] & (1 << 8)) != 0

    mu22_kin = (m_mu["pt"] > 23) & (abs(m_mu["eta"]) < 2.1)
    mu19_kin = (
        (m_mu["pt"] > 20) &
        (m_mu["pt"] < 23) &
        (abs(m_mu["eta"]) < 2.1)
    )

    # loop uses `if HLTIsoMu22Fired: ... elif HLTIsoMu19Tau20Fired: ...`
    # so the mu19tau branch is only ever reached when mu22 did NOT fire --
    # it is not an independent OR condition. Flat per-event booleans
    # broadcast automatically into the (event, muon, trigobj) structure.
    match_mu22    = dr_ok & filter_mu22    & mu22_kin & HLTIsoMu22Fired
    match_mu19tau = dr_ok & filter_mu19tau & mu19_kin & (~HLTIsoMu22Fired & HLTIsoMu19Tau20Fired)

    # per-muon boolean, shape (event, muon) -- reduces the trigobj axis
    muMatched = (match_mu22 | match_mu19tau).any()

    # =================================================================
    # 7. Tau trigger matching (event-level -- doesn't depend on which
    #    muon we're looking at, exactly like in the loop)
    # =================================================================
    tau_trigsel = (
        (br["Tau_pt"] >= 20) &
        (abs(br["Tau_eta"]) <= 2.1) &
        (br["Tau_idDeepTau2017v2p1VSjet"] >= 16) &
        (br["Tau_idDeepTau2017v2p1VSe"]   >= 1)  &
        (br["Tau_idDeepTau2017v2p1VSmu"]  >= 8)  &
        (
            (br["Tau_decayMode"] == 0) |
            (br["Tau_decayMode"] == 1) |
            (br["Tau_decayMode"] == 10)
        )
    )

    tau_cands = ak.JaggedArray.zip(eta=br["Tau_eta"][tau_trigsel],
                                         phi=br["Tau_phi"][tau_trigsel])

    trig_tau_mask = abs(br["TrigObj_id"]) == 15
    trig_tau = ak.JaggedArray.zip(
        eta=br["TrigObj_eta"][trig_tau_mask],
        phi=br["TrigObj_phi"][trig_tau_mask],
        filterBits=br["TrigObj_filterBits"][trig_tau_mask],
    )

    pairs_tau = tau_cands.cross(trig_tau, nested=True)
    ta, tt = pairs_tau.i0, pairs_tau.i1
    dr_tau = dR(ta["eta"], ta["phi"], tt["eta"], tt["phi"])

    filter_tau20     = (tt["filterBits"] & (1 << 8)) != 0
    tau_trig_matched = ((dr_tau < 0.5) & filter_tau20).any()  # per tau candidate
    has_matched_tau  = tau_trig_matched.any()                  # per event

    # only REQUIRED when the mu19+tau20 path is the one that fired
    event_tau_ok = (~HLTIsoMu19Tau20Fired) | has_matched_tau   # flat numpy, shape: (event,)

    # =================================================================
    # 8. Final per-(event, muon) fill mask -- combines everything the
    #    loop checks before calling hMt.Fill(), muon by muon. Flat
    #    per-event arrays broadcast automatically against the jagged
    #    muon-level mask.
    # =================================================================
    final_muon_mask = muon_strict & muMatched & event_hlt_ok & event_tau_ok

    # =================================================================
    # 9. Build MT per surviving (event, muon). MET_pt/MET_phi are flat
    #    per-event arrays and broadcast automatically across however
    #    many muons pass in that event.
    # =================================================================
    sel_mu_pt  = br["Muon_pt"][final_muon_mask]
    sel_mu_phi = br["Muon_phi"][final_muon_mask]

    mt = np.sqrt(
        2.0 * sel_mu_pt * br["MET_pt"] *
        (1.0 - np.cos(sel_mu_phi - br["MET_phi"]))
    )

    if is_mc:
        # xsec and luminosity are now function parameters -- pass the
        # correct cross-section for whichever sample this is (signal vs.
        # each background process). Normalize by genEventSumw (the TOTAL
        # sum of genWeight across the whole sample, all files combined,
        # not this file's local value) rather than raw event count, since
        # genWeight is not simply +-1 for this generator.
        event_weight = br["genWeight"] * xsec * luminosity / float(genEventSumw)
    else:
        # use the CURRENT (post-filtering) event count, not the original
        # numEvents -- by this point br has already been through the
        # MET-filter / b-jet-veto / 3rd-lepton-veto event_filter() calls,
        # so its length no longer matches the pre-filter total.
        event_weight = np.ones(len(br["MET_pt"]))

    # broadcast the flat per-event weight into the same jagged shape as mt,
    # then flatten both together so mt[i] and weight[i] always line up
    weight_jagged = sel_mu_pt.zeros_like() + event_weight

    mt_flat     = mt.flatten()
    weight_flat = weight_jagged.flatten()

    # =================================================================
    # 10. Fill ROOT histogram -- one Fill() call per qualifying muon,
    #     matching the loop's no-break behavior exactly.
    # =================================================================
    ROOT.gDirectory.cd("PyROOT:/")
    hMt = ROOT.TH1F(hist_name, "Transverse Mass; MT(Mu,MET) [GeV]; Events", 100, 0.0, 150.0)

    for i in range(len(mt_flat)):
        hMt.Fill(mt_flat[i], weight_flat[i])

    print "Histogram filled with %d entries" % len(mt_flat)
    return hMt
