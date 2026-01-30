import ROOT
import os


# Load files
#Combining parts of data files
os.system("hadd hMTData_combined -f hMTData_combinedPart2.root hMTData_combinedPart1.root")
fD = ROOT.TFile("hMtData_combined.root")

fS = ROOT.TFile("hMtSignal_combined1.root")

#Combining parts of background files
os.system("hadd hMtBackground_combined -f hMtBackground_combinedPart1.root hMtBackground_combinedPart2.root")
fB = ROOT.TFile("hMtBackground_combined.root")

# Get histograms
hMtData = fD.Get("hMtDataTotal")
hMtSignal = fS.Get("hMtSignalTotal")
hMtBackground = fB.Get("hMtBackgroundTotal")


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
hMtSignal.Draw("HIST SAME", "nostack")
hMtBackground.Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
leg.AddEntry(hMtData, "Data", "lep")
leg.AddEntry(hMtSignal, "Signal MC (scaled)", "l")
leg.AddEntry(hMtBackground, "Background MC (scaled)", "l")
leg.Draw()

# Save plots
c.SaveAs("combined.pdf")
