# -*- coding: utf-8 -*-
import ROOT
from analysisfxn import Analysis, Analysis_up
import awkward as ak


TEST_FILE = "/data/bmc/18C499E2-9229-DC45-A68F-6A40600A5C8E.root"
IS_MC = True

# -----------------------------------------------
# Run both versions
# -----------------------------------------------
print "Running FOR-LOOP version..."
hLoop = Analysis(TEST_FILE, hist_name="hLoop", is_mc=IS_MC)

print "\nRunning VECTORIZED version..."
hVec = Analysis_up(TEST_FILE, hist_name="hVec", is_mc=IS_MC)

# -----------------------------------------------
# Save both to a ROOT file
# -----------------------------------------------
out = ROOT.TFile("comparison.root", "RECREATE")
hLoop.Write("hLoop")
hVec.Write("hVec")
out.Close()
print "\nSaved to comparison.root"

# -----------------------------------------------
# Numerical check: bin-by-bin differences
# -----------------------------------------------
print "\n--- Bin-by-bin check ---"
print "%-10s %-15s %-15s %-10s" % ("Bin", "Loop", "Vec", "Diff%")
max_diff = 0.0
for i in range(1, 101):
    loop_val = hLoop.GetBinContent(i)
    vec_val  = hVec.GetBinContent(i)
    if loop_val > 0:
        diff = abs(vec_val - loop_val) / loop_val * 100.0
        if diff > max_diff:
            max_diff = diff
        if diff > 1.0:
            print "%-10d %-15.4f %-15.4f %-10.2f%%" % (i, loop_val, vec_val, diff)

print "\nMax bin difference: %.4f%%" % max_diff
if max_diff < 0.01:
    print "PASS: histograms match within 0.01%%"
else:
    print "FAIL: discrepancy found - check logic above"

# -----------------------------------------------
# Plot: overlay + ratio
# -----------------------------------------------
hLoop.SetLineColor(ROOT.kBlue)
hLoop.SetLineWidth(2)
hLoop.SetLineStyle(1)

hVec.SetLineColor(ROOT.kRed)
hVec.SetLineWidth(2)
hVec.SetLineStyle(2)

c = ROOT.TCanvas("c", "Comparison", 800, 900)

top = ROOT.TPad("top", "top", 0, 0.35, 1, 1.0)
top.SetBottomMargin(0.02)
top.Draw()
top.cd()

hLoop.Draw("HIST")
hVec.Draw("HIST SAME")

leg = ROOT.TLegend(0.65, 0.70, 0.88, 0.88)
leg.AddEntry(hLoop, "For-loop", "l")
leg.AddEntry(hVec,  "Vectorized", "l")
leg.Draw()

c.cd()
bot = ROOT.TPad("bot", "bot", 0, 0.0, 1, 0.35)
bot.SetTopMargin(0.02)
bot.SetBottomMargin(0.3)
bot.Draw()
bot.cd()

hRatio = hVec.Clone("hRatio")
hRatio.Divide(hLoop)
hRatio.SetTitle("")
hRatio.GetYaxis().SetTitle("Vec / Loop")
hRatio.GetYaxis().SetRangeUser(0.8, 1.2)
hRatio.GetYaxis().SetNdivisions(5)
hRatio.GetXaxis().SetTitle("MT(Mu,MET) [GeV]")
hRatio.GetXaxis().SetTitleSize(0.1)
hRatio.GetXaxis().SetLabelSize(0.08)
hRatio.GetYaxis().SetTitleSize(0.08)
hRatio.GetYaxis().SetLabelSize(0.07)
hRatio.Draw("EP")

line = ROOT.TLine(0, 1.0, 150, 1.0)
line.SetLineColor(ROOT.kBlack)
line.SetLineStyle(2)
line.Draw("SAME")

c.SaveAs("comparison.pdf")
print "\nSaved comparison.pdf"
