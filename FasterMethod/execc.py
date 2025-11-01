import ROOT

# Load files
fD = ROOT.TFile("hMtData.root")
fS = ROOT.TFile("hMtSignal.root")
fB = ROOT.TFile("hMtBackground.root")

# Get histograms
hMtData = fD.Get("hMtData")
hMtSignal = fS.Get("hMtSignal")
hMtBackground = fB.Get("hMtBackground")

# Scale MC to Data for visual comparison
data_int = hMtData.Integral()
signal_int = hMtSignal.Integral()
bkg_int = hMtBackground.Integral()

if signal_int > 0:
    hMtSignal.Scale(data_int / signal_int)

if bkg_int > 0:
    hMtBackground.Scale(data_int / bkg_int)

# Style settings
hMtData.SetLineColor(ROOT.kBlack)
hMtData.SetMarkerStyle(20)

hMtSignal.SetLineColor(ROOT.kRed)
hMtSignal.SetLineWidth(2)

hMtBackground.SetLineColor(ROOT.kBlue)
hMtBackground.SetLineStyle(2)
hMtBackground.SetLineWidth(2)

# Draw
c = ROOT.TCanvas("c", "MT Plots", 800, 700)
c.SetLogy()  # Important to view small background

hMtData.Draw("E")  # Data points
hMtSignal.Draw("HIST SAME")
hMtBackground.Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
leg.AddEntry(hMtData, "Data", "lep")
leg.AddEntry(hMtSignal, "Signal MC (scaled)", "l")
leg.AddEntry(hMtBackground, "Background MC (scaled)", "l")
leg.Draw()

# Save plots
c.SaveAs("combined.pdf")
