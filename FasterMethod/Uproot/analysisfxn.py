# -*- coding: utf-8 -*-

import ROOT
import uproot
import awkward as ak
import numpy as np


# ============================================================
# CMS 2016 PILEUP WEIGHTING
# ============================================================

PU_DIR = (
    "/code/CMSSW_10_6_30/src/PhysicsTools/NanoAODTools/"
    "python/postprocessing/data/pileup/"
)

PU_MC_FILE = PU_DIR + "pileup_profile_Summer16.root"
PU_DATA_FILE = PU_DIR + "PileupData_GoldenJSON_Full2016.root"


# Load the official CMS 2016 PU distributions
pu_mc_file = ROOT.TFile.Open(PU_MC_FILE)
pu_data_file = ROOT.TFile.Open(PU_DATA_FILE)

pu_mc_hist = pu_mc_file.Get("pu_mc")
pu_data_hist = pu_data_file.Get("pileup")

if not pu_mc_hist:
    raise RuntimeError(
        "Could not load pu_mc from %s" % PU_MC_FILE
    )

if not pu_data_hist:
    raise RuntimeError(
        "Could not load pileup from %s" % PU_DATA_FILE
    )


# Normalize MC and data distributions
pu_mc_integral = pu_mc_hist.Integral()
pu_data_integral = pu_data_hist.Integral()


# Precalculate PU weights:
#
#       data PU PDF
# w = ----------------
#       MC PU PDF
#
pu_weights = {}

for ibin in range(1, pu_mc_hist.GetNbinsX() + 1):

    mc_value = pu_mc_hist.GetBinContent(ibin)
    data_value = pu_data_hist.GetBinContent(ibin)

    # Follow CMS convention of not allowing negative
    # probability contents.
    if mc_value < 0:
        mc_value = 0.0

    if data_value < 0:
        data_value = 0.0

    if pu_mc_integral > 0:
        mc_pdf = mc_value / pu_mc_integral
    else:
        mc_pdf = 0.0

    if pu_data_integral > 0:
        data_pdf = data_value / pu_data_integral
    else:
        data_pdf = 0.0

    if mc_pdf > 0:
        pu_weights[ibin] = data_pdf / mc_pdf
    else:
        pu_weights[ibin] = 1.0


def get_pu_weight(ntrue):

    """
    Get the 2016 CMS pileup weight corresponding to
    Pileup_nTrueInt.
    """

    if ntrue is None:
        return 1.0

    try:
        ntrue = float(ntrue)
    except:
        return 1.0

    ibin = pu_mc_hist.FindBin(ntrue)

    # Outside the PU histogram range
    if ibin < 1 or ibin > pu_mc_hist.GetNbinsX():
        return 1.0

    return pu_weights.get(ibin, 1.0)


# ============================================================
# TOTAL GEN EVENT SUMW
# ============================================================

def get_total_gen_event_sumw(file_list):

    total = 0.0

    for fp in file_list:
        total += uproot.open(fp)["Runs"].array(
            "genEventSumw"
        ).sum()

    return total


# ============================================================
# ANALYSIS
# ============================================================

def Analysis_up(filePath, hist_name="hMvis", is_mc=False,
                xsec=6077.22, luminosity=35920.0,
                genEventSumw=None):

    tree = uproot.open(filePath)["Events"]

    branches = [
        "Muon_pt", "Muon_eta", "Muon_phi", "Muon_charge",
        "Muon_dz", "Muon_dxy", "Muon_mediumId",
        "Muon_pfRelIso04_all",

        "MET_pt", "MET_phi",

        "Jet_pt", "Jet_eta", "Jet_btagDeepB",

        "Electron_pt", "Electron_eta", "Electron_phi",
        "Electron_dxy", "Electron_dz",
        "Electron_mvaFall17V2noIso_WP90",
        "Electron_convVeto", "Electron_lostHits",
        "Electron_pfRelIso03_all",

        "Tau_pt", "Tau_eta", "Tau_phi", "Tau_charge",
        "Tau_dz",
        "Tau_idDeepTau2017v2p1VSe",
        "Tau_idDeepTau2017v2p1VSmu",
        "Tau_idDeepTau2017v2p1VSjet",
        "Tau_decayMode",
        "Tau_rawDeepTau2017v2p1VSjet",

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
        "TrigObj_filterBits"
    ]

    # --------------------------------------------------------
    # MC-only branches
    # --------------------------------------------------------

    if is_mc:
        branches.append("genWeight")
        branches.append("Pileup_nTrueInt")

    br = tree.arrays(
        branches,
        namedecode="utf-8"
    )

    if is_mc and genEventSumw is None:
        raise ValueError(
            "genEventSumw must be supplied for MC"
        )

    # ========================================================
    # DELTA R
    # ========================================================

    def dR(eta1, phi1, eta2, phi2):

        deta = eta1 - eta2

        dphi = np.arctan2(
            np.sin(phi1 - phi2),
            np.cos(phi1 - phi2)
        )

        return np.sqrt(deta**2 + dphi**2)

    # ========================================================
    # EVENT FILTER
    # ========================================================

    def event_filter(mask):

        for key in br:
            br[key] = br[key][mask]

    # ========================================================
    # MET FILTERS
    # ========================================================

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

    # ========================================================
    # b-JET VETO
    # ========================================================

    jets_good = (
        (br["Jet_pt"] > 20) &
        (abs(br["Jet_eta"]) < 2.4)
    )

    bjet_veto = (
        br["Jet_btagDeepB"][jets_good] > 0.6321
    ).any()

    event_filter(~bjet_veto)

    # ========================================================
    # ELECTRON VETO
    # ========================================================

    tau_veto = (
        (br["Tau_pt"] > 20) &
        (abs(br["Tau_eta"]) < 2.3)
    )

    ele_good = (
        (br["Electron_pt"] >= 10) &
        (abs(br["Electron_eta"]) <= 2.5) &
        (abs(br["Electron_dz"]) <= 0.2) &
        (abs(br["Electron_dxy"]) <= 0.045) &
        (br["Electron_mvaFall17V2noIso_WP90"] == True) &
        (br["Electron_convVeto"] == True) &
        (br["Electron_lostHits"] <= 1) &
        (br["Electron_pfRelIso03_all"] < 0.15)
    )

    electrons = ak.JaggedArray.zip(
        eta=br["Electron_eta"][ele_good],
        phi=br["Electron_phi"][ele_good]
    )

    taus = ak.JaggedArray.zip(
        eta=br["Tau_eta"][tau_veto],
        phi=br["Tau_phi"][tau_veto]
    )

    pairs = electrons.cross(
        taus,
        nested=True
    )

    e = pairs.i0
    t = pairs.i1

    separated = (
        dR(
            e["eta"], e["phi"],
            t["eta"], t["phi"]
        ) > 0.5
    ).all()

    event_filter(~separated.any())

    # ========================================================
    # MUON SELECTION
    # ========================================================

    muon_strict = (
        (br["Muon_pt"] > 10) &
        (abs(br["Muon_eta"]) < 2.4) &
        (abs(br["Muon_dz"]) < 0.2) &
        (abs(br["Muon_dxy"]) < 0.045) &
        (br["Muon_mediumId"] == True) &
        (br["Muon_pfRelIso04_all"] < 0.15) &
        (br["PV_npvsGood"] >= 1)
    )

    # ========================================================
    # TRIGGER
    # ========================================================
    # USE ONLY MU-TAU TRIGGER
    hlt_mu19tau = br["HLT_IsoMu19_eta2p1_LooseIsoPFTau20"]

    mu = ak.JaggedArray.zip(
        pt=br["Muon_pt"],
        eta=br["Muon_eta"],
        phi=br["Muon_phi"],
        charge=br["Muon_charge"]
    )

    trig_mu_mask = abs(br["TrigObj_id"]) == 13

    trig_mu = ak.JaggedArray.zip(
        eta=br["TrigObj_eta"][trig_mu_mask],
        phi=br["TrigObj_phi"][trig_mu_mask],
        filterBits=br["TrigObj_filterBits"][trig_mu_mask]
    )

    pairs = mu.cross(
        trig_mu,
        nested=True
    )

    m = pairs.i0
    tm = pairs.i1

    dr = dR(
        m["eta"], m["phi"],
        tm["eta"], tm["phi"]
    ) < 0.5

    f19 = (
        (tm["filterBits"] & (1 << 8)) != 0
    )

    mu19_kin = (
        (m["pt"] > 20) &
        (m["pt"] < 23) &
        (abs(m["eta"]) < 2.1)
    )

    match19 = (
        dr &
        f19 &
        mu19_kin &
        hlt_mu19tau
    )

    mu19 = match19.any()
    # ========================================================
    # TAU SELECTION
    # ========================================================

    tau_base = (
        (abs(br["Tau_eta"]) < 2.3) &
        (abs(br["Tau_dz"]) < 0.2) &
        (br["Tau_idDeepTau2017v2p1VSe"] >= 1) &
        (br["Tau_idDeepTau2017v2p1VSmu"] >= 8) &
        (br["Tau_idDeepTau2017v2p1VSjet"] >= 16) &
        (
            (br["Tau_decayMode"] == 0) |
            (br["Tau_decayMode"] == 1) |
            (br["Tau_decayMode"] == 10)
        )
    )

    tau = ak.JaggedArray.zip(
        pt=br["Tau_pt"][tau_base],
        eta=br["Tau_eta"][tau_base],
        #phi=br["Tau_phi"][tau_base],
        charge=br["Tau_charge"][tau_base],
        djet=br["Tau_rawDeepTau2017v2p1VSjet"][tau_base]
        # iMPLEMENT ALL THE DEFINED VARIAABLES IN TAU_BASE SPECIALLLY ANTI ELECTRON AND ANTI MUON. CHECK THE CHARGE AAND PHI
    )

    # ========================================================
    # TAU TRIGGER MATCHING
    # ========================================================

    trig_tau_mask = abs(br["TrigObj_id"]) == 15

    trig_tau = ak.JaggedArray.zip(
        eta=br["TrigObj_eta"][trig_tau_mask],
        phi=br["TrigObj_phi"][trig_tau_mask],
        filterBits=br["TrigObj_filterBits"][trig_tau_mask]
    )

    pairs = tau.cross(
        trig_tau,
        nested=True
    )

    ta = pairs.i0
    tt = pairs.i1

    tau_match = (
        (
            dR(
                ta["eta"], ta["phi"],
                tt["eta"], tt["phi"]
            ) < 0.5
        ) &
        (
            (tt["filterBits"] & (1 << 8)) != 0
        )
    ).any()

    # ========================================================
    # OBJECTS
    # ========================================================

    mu = ak.JaggedArray.zip(
        pt=br["Muon_pt"],
        eta=br["Muon_eta"],
        phi=br["Muon_phi"],
        charge=br["Muon_charge"],

        valid=(
            muon_strict &
            (mu22 | mu19) &
            (hlt_mu22 | hlt_mu19tau)
        ),

        mu22=mu22,
        mu19=mu19
    )

    tau = ak.JaggedArray.zip(
        pt=tau["pt"],
        eta=tau["eta"],
        phi=tau["phi"],
        charge=tau["charge"],
        djet=tau["djet"],
        trigmatched=tau_match
    )

    # ========================================================
    # MU-TAU PAIRS
    # ========================================================

    pairs = mu.cross(tau)

    pmu = pairs.i0
    ptau = pairs.i1

    pair_dr = (
        dR(
            pmu["eta"], pmu["phi"],
            ptau["eta"], ptau["phi"]
        ) > 0.5
    )

    opposite_sign = (
        pmu["charge"] * ptau["charge"] < 0
    )

    tau_pt = (
        (pmu["mu19"] & (ptau["pt"] > 32)) |
        (pmu["mu22"] & (ptau["pt"] > 30))
    )

    tau_trigger = (
        pmu["mu22"] |
        (
            pmu["mu19"] &
            ptau["trigmatched"]
        )
    )

    valid = (
        pmu["valid"] &
        pair_dr &
        opposite_sign &
        tau_pt &
        tau_trigger
    )

    pmu = pmu[valid]
    ptau = ptau[valid]

    # ========================================================
    # SELECT BEST PAIR
    # ========================================================

    best = ptau["djet"].argmax()

    mu_pt = pmu["pt"][best]
    mu_eta = pmu["eta"][best]
    mu_phi = pmu["phi"][best]

    tau_pt = ptau["pt"][best]
    tau_eta = ptau["eta"][best]
    tau_phi = ptau["phi"][best]

    has_pair = tau_pt.counts > 0

    mu_pt = mu_pt[has_pair]
    mu_eta = mu_eta[has_pair]
    mu_phi = mu_phi[has_pair]

    tau_pt = tau_pt[has_pair]
    tau_eta = tau_eta[has_pair]
    tau_phi = tau_phi[has_pair]

    # ========================================================
    # VISIBLE MASS
    # ========================================================

    deta = mu_eta - tau_eta

    dphi = np.arctan2(
        np.sin(mu_phi - tau_phi),
        np.cos(mu_phi - tau_phi)
    )

    mvis = np.sqrt(
        2.0 *
        mu_pt *
        tau_pt *
        (
            np.cosh(deta) -
            np.cos(dphi)
        )
    )

    # ========================================================
    # MC WEIGHTS
    # ========================================================

    if is_mc:

        # Generator weight
        gen_weight = br["genWeight"][has_pair]

        # Number of true interactions
        pu_ntrue = br["Pileup_nTrueInt"][has_pair]

        # Cross-section normalization
        normalization = (
            xsec *
            luminosity /
            float(genEventSumw)
        )

        # ----------------------------------------------------
        # Final event weight:
        #
        # genWeight
        # × cross-section normalization
        # × pileup weight
        # ----------------------------------------------------

        weights = []

        for i in range(len(mvis)):

            pu_weight = get_pu_weight(
                pu_ntrue[i]
            )

            weight = (
                gen_weight[i] *
                normalization *
                pu_weight
            )

            weights.append(weight)

        weights = np.asarray(weights)

    else:

        weights = np.ones(len(mvis))

    # ========================================================
    # FLATTEN
    # ========================================================

    mvis = mvis.flatten()
    weights = weights.flatten()

    # ========================================================
    # HISTOGRAM
    # ========================================================

    ROOT.gDirectory.cd("PyROOT:/")

    hMvis = ROOT.TH1F(
        hist_name,
        "Visible Mass;M_{vis}(#mu,#tau_{h}) [GeV];Events",
        100,
        0.0,
        150.0
    )

    for i in range(len(mvis)):

        hMvis.Fill(
            float(mvis[i]),
            float(weights[i])
        )

    return hMvis
