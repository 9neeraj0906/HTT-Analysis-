import ROOT
import os

# Data, signal, and background were each produced from a single full run
# this time (not split into parts), so there's nothing to hadd anymore --
# just open each output file directly.
os.system("hadd hMtData_combined.root -f hMTData_combinedPart1.root hMTData_combinedPart2.root")

fD = ROOT.TFile("hMtData_combined.root")
fS = ROOT.TFile("hMtSignal_combined.root")
fB = ROOT.TFile("hMtBackground_combined21.root")

# Get histograms
hMtData = fD.Get("hMTDataTotal")
hMtSignal = fS.Get("hMtSignalTotalmc")
hMtBackground = fB.Get("hMtBackgroundTotalmc")

# Style settings
hMtData.SetLineColor(ROOT.kBlack)
hMtData.SetMarkerStyle(20)
hMtSignal.SetLineColor(ROOT.kRed)
hMtSignal.SetLineStyle(2)
hMtSignal.SetLineWidth(2)
hMtBackground.SetLineColor(ROOT.kBlue)
hMtBackground.SetLineStyle(2)
hMtBackground.SetLineWidth(2)

# Draw
c = ROOT.TCanvas("c", "MT Plots", 800, 700)
hMtData.Draw("E")  # Data points
hMtBackground.Draw("HIST SAME")
hMtSignal.Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
leg.AddEntry(hMtData, "Data", "lep")
leg.AddEntry(hMtSignal, "Signal MC (scaled)", "l")
leg.AddEntry(hMtBackground, "Background MC (scaled)", "l")
leg.Draw()

# Save plot
c.SaveAs("combined3.pdf")
print "Background raw entries:", hMtBackground.GetEntries()
print "Background weighted integral:", hMtBackground.Integral()
print "Data raw entries:", hMtData.GetEntries()
