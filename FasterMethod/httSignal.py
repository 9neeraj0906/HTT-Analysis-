import ROOT
import math
import time

from analFxn import Analysis

signalMCPath = "/data/smc/8544982A-CEC6-6B4D-B7F8-9B5AF213B725.root"
hMtSignal = Analysis(signalMCPath, hist_name="hMtSignal", is_mc=True)
hMtSignal.SaveAs("hMtSignal.root")
