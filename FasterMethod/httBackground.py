import ROOT
import time
import math
import sys, os



from analFxn import Analysis





backgroundMCPath = [


    "/data/bmc/0718C107-8960-6B44-B96A-C60D53D52A95.root",
    "/data/bmc/4578E947-084C-C946-9B8D-1B45A126DCED.root",
    "/data/bmc/4BD16DC1-5991-6549-9B08-17B67F92A986.root",
    "/data/bmc/0082C29D-E74C-024A-BE9B-97B29EE7A4A2.root",
    "/data/bmc/500A8824-5BD4-5741-AECD-FE1A7FD0098D.root",
    "/data/bmc/52DA0012-7ABA-F94A-A7D3-2FC06E9DF43B.root"
    

]
output_file = ROOT.TFile("hMtBackground_combined.root", "RECREATE")
hMtBackgroundTotal = None


for i in range(len(backgroundMCPath)):
    hMtBackground = Analysis(backgroundMCPath[i], hist_name="htemp%d" % i, is_mc=False)
    if hMtBackgroundTotal is None:
        hMtBackgroundTotal= hMtBackground.Clone("hMtBackgroundTotal")
    else:
        hMtBackgroundTotal.Add(hMtBackground)


output_file.cd()
hMtBackgroundTotal.Write()

output_file.Close()
