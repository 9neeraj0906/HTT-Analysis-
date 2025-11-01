import ROOT
import time
import math
import sys, os



from analFxn import Analysis

backgroundMCPath = "/data/bmc/0718C107-8960-6B44-B96A-C60D53D52A95.root"
hMtBackground = Analysis(backgroundMCPath, hist_name="hMtBackground", is_mc=True)
hMtBackground.SaveAs("hMtBackground.root")
