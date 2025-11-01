import ROOT
import math
import time



from analFxn import Analysis
dataPath = "/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root"
hMtData = Analysis(dataPath, hist_name="hMtData", is_mc = False)
hMtData.SaveAs("hMtData.root")
