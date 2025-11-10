import ROOT
import math
import time

from analFxn import Analysis

signalMCPath = [
    "/data/smc/03F42243-AFD4-F348-92E9-6E6AC3B36FD6.root",
    "/data/smc/43C9FA0D-05BF-314C-8128-9B77E519ABAB.root",
    "/data/smc/650B8FFC-7702-D042-BC13-6299339E0618.root",
    "/data/smc/8544982A-CEC6-6B4D-B7F8-9B5AF213B725.root",
    "/data/smc/85C5EC80-3FBE-BD42-9413-37985082F475.root",
    "/data/smc/89535AE7-0E8F-E442-867D-A67F300684A0.root",
    "/data/smc/A69C2B66-ADEF-4442-A710-780DA9FF78A3.root"
]


output_file = ROOT.TFile("hMtSignal_combined.root", "RECREATE")
hMtSignalTotal = None


for i in range(len(signalMCPath)):
    hMtSignal = Analysis(signalMCPath[i], hist_name="htemp%d" % i, is_mc=False)
    if hMtSignalTotal is None:
        hMtSignalTotal= hMtSignal.Clone("hMtSignalTotal")
    else:
        hMtSignalTotal.Add(hMtSignal)


output_file.cd()
hMtSignalTotal.Write()

output_file.Close()
