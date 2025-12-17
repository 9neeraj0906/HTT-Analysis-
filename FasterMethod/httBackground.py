import ROOT
import time
import math
import sys, os



from analFxn import Analysis





backgroundMCPath = [
    "/data/bmc/0082C29D-E74C-024A-BE9B-97B29EE7A4A2.root",
    "/data/bmc/0718C107-8960-6B44-B96A-C60D53D52A95.root",
    "/data/bmc/0A995C8D-0ACC-7E45-8C7A-16AC51398938.root",
    "/data/bmc/4578E947-084C-C946-9B8D-1B45A126DCED.root",
    "/data/bmc/4BD16DC1-5991-6549-9B08-17B67F92A986.root",
    "/data/bmc/500A8824-5BD4-5741-AECD-FE1A7FD0098D.root",
    "/data/bmc/52DA0012-7ABA-F94A-A7D3-2FC06E9DF43B.root",
    "/data/bmc/5544E66D-B271-1649-A85B-C08F01375C57.root",
    "/data/bmc/69C45E78-12E9-3D41-B5EF-DAF915460237.root",
    "/data/bmc/6C3CD8A5-A288-724E-BF9F-BAAC46A4C139.root",
    "/data/bmc/6D7ABB90-DD75-964F-A5CD-14D1B21FAB36.root",
    "/data/bmc/7169FB77-0571-7641-B0D4-B4270ABF7B39.root",
    "/data/bmc/74356482-D5C8-D546-85D3-416B9012D25D.root",
    "/data/bmc/7AAA7AAD-2378-E547-B014-759BBEAA4C42.root",
    "/data/bmc/9400079E-20CF-1449-BB3C-A7FCDE90F9BF.root",
    "/data/bmc/94C0AF7B-6727-F246-8E57-384384DADE87.root",
    "/data/bmc/969B7293-8E74-2A40-8D14-B33204BD6401.root",
    "/data/bmc/9D579D34-7A90-DF4E-BD6B-A0E9565CF455.root",
    "/data/bmc/9EDE17FC-CAB6-A146-9714-F1C496F1CBE2.root",
    "/data/bmc/A17A523B-378F-9A4F-89D6-71C80D612F4E.root"







]
output_file = ROOT.TFile("hMtBackground_combinedPart1.root", "RECREATE")
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
