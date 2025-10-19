import ROOT

ROOT.ROOT.EnableImplicitMT()

# load custom cpp tool box 
ROOT.gInterpreter.Declare('#include "HTT_analysis.h"')

dataPath = "/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root" # path to muon>
smcPath = "/data/smc/8544982A-CEC6-6B4D-B7F8-9B5AF213B725.root" # path to sign>
bmcPath = "/data/bmc/0718C107-8960-6B44-B96A-C60D53D52A95.root" # path to back>

#load data
dfData = ROOT.RDataFrame("Events", dataPath)
dfSignal = ROOT.RDataFrame("Events", smcPath)
dfBackground = ROOT.RDataFrame("Events", bmcPath)



def applyEventSelection(df):
    """
    This function takes a dataframe and applies all the cuts
    for the H->tau(mu)tau(had) analysis for 2016 data.
    """
    dfFiltered = df.Filter(
    "hasGoodMuon(Muon_pt) && hasGoodTau(boostedTau_pt, boostedTau_eta, boostedTau_idMVAnewDM2017v2, boostedTau_decayMode)",
    "Event has one good muon (pT > 20) and one good tau"
)
    dfOppositeCharge = dfFiltered.Filter("Muon_charge[0] * boostedTau_charge[0] < 0", "Opposite charge requirement")

    dfBjetVeto = dfOppositeCharge.Filter("Sum(Jet_btagCSVV2 > 0.8484) == 0", "B-jet veto")

    #calculate the final observable: Trnsverse Mass 
    # FINAL CORRECTION: Using the PuppiMET branches
    dfFinal = dfBjetVeto.Define(
        "transverse_mass",
        "computeTransverseMass(Muon_pt[0], Muon_phi[0], PuppiMET_pt, PuppiMET_phi)"
    )
    return dfFinal


df_data_final = applyEventSelection(dfData)
hist_data = df_data_final.Histo1D(("data_obs", "Transverse Mass;m_{T} [GeV];Events", 40, 0, 200), "transverse_mass")
#print(df_data_final)
df_signal_final = applyEventSelection(dfSignal)
weight_signal = 1
hist_signal = df_signal_final.Define("weight", str(weight_signal)).Histo1D(("signal", "Signal (H#rightarrow#tau#tau)", 40, 0, 200), "transverse_mass", "weight")

df_bkg_final = applyEventSelection(dfBackground)
weight_bkg = 1.23
hist_bkg = df_bkg_final.Define("weight", str(weight_bkg)).Histo1D(("background", "Background (Z#rightarrow#tau#tau)", 40, 0, 200), "transverse_mass", "weight")

print("Starting event loop for Data...")
h_data = hist_data.GetValue()

print("Starting event loop for Signal...")
h_signal = hist_signal.GetValue()

print("Starting event loop for Background...")
h_bkg = hist_bkg.GetValue()

print("All event loops finished!")

#plot(signal_mu, background_mu, data_mu, "m_{4#mu} (GeV)", "higgs_4mu.pdf")
print("plotting results")

#canvas 
c1 = ROOT.TCanvas("c", "Transverse Mass", 800,700)
# styling background 

h_bkg.SetMarkerStyle(20) # 20 is a circle
h_bkg.SetMarkerColor(ROOT.kBlue)
h_bkg.SetLineColor(ROOT.kBlue)

# styling signal 
h_signal.SetLineColor(ROOT.kRed)
h_signal.SetLineWidth(2)

#styling data 
h_data.SetMarkerStyle(29) # 29 is a filled star
h_data.SetMarkerSize(1.5)
h_data.SetMarkerColor(ROOT.kBlack)




h_bkg.Draw("P")
h_signal.Draw("HIST SAME")
h_data.Draw("P SAME")
c1.SaveAs("transverse_mass_minimal.pdf")
print("Plots saved as transverse_mass_minimal.pdf")
