import ROOT
import os

# Combine the two data files
os.system(
    "hadd hMvisData_combined.root -f "
    "hMvisData_combinedPart1.root "
    "hMvisData_combinedPart2.root"
)

# Open files
fD = ROOT.TFile("hMvisData_combined.root")
fB = ROOT.TFile("hMvisBackground_combined.root")

# Get histograms
hMvisData = fD.Get("hMvisDataTotal")
hMvisBackground = fB.Get("hMvisBackgroundTotal")

# Style settings
hMvisData.SetLineColor(ROOT.kBlack)
hMvisData.SetMarkerStyle(20)

hMvisBackground.SetLineColor(ROOT.kBlue)
hMvisBackground.SetLineStyle(2)
hMvisBackground.SetLineWidth(2)

# Draw
c = ROOT.TCanvas("c", "Visible Mass Plot", 800, 700)

hMvisData.Draw("E")
hMvisBackground.Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.60, 0.72, 0.88, 0.85)

leg.AddEntry(
    hMvisData,
    "Data",
    "lep"
)

leg.AddEntry(
    hMvisBackground,
    "Background MC (scaled)",
    "l"
)

leg.Draw()

# Save plot
c.SaveAs("combined_visible_mass.pdf")

# Print statistics
print "Background raw entries:", hMvisBackground.GetEntries()
print "Background weighted integral:", hMvisBackground.Integral()

print "Data raw entries:", hMvisData.GetEntries()
print "Data integral:", hMvisData.Integral()
