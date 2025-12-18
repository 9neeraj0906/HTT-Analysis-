import ROOT
import time
import math
import sys, os



from analFxn import Analysis





backgroundMCPath = [
    "/data/bmc/18C499E2-9229-DC45-A68F-6A40600A5C8E.root",
    "/data/bmc/18DA8B20-842D-0741-BAFB-1A1C7C522921.root",
    "/data/bmc/196F07F8-8EC2-6A48-9036-8B7804B71A1F.root",
    "/data/bmc/1D85E077-414B-0349-A901-CD6B3B3C349D.root",
    "/data/bmc/226E10F4-29BC-CF45-B8BF-31B13B4EC1CE.root",
    "/data/bmc/4046329B-AF7C-6142-9B76-7EE37FB93DCD.root",
    "/data/bmc/440A3E55-6EE5-3642-B16C-F7FC5FAB5081.root",
    "/data/bmc/A2058F3B-C2A2-4E43-B963-6C5A4DC8084B.root",
    "/data/bmc/A20C5B1A-2CFF-7B42-82A4-06347D0856A7.root",
    "/data/bmc/AE5A627E-C8D7-974F-9698-C136D2BFFF1D.root",
    "/data/bmc/AF650968-C8A8-9043-9BBB-F83AE7726FEE.root",
    "/data/bmc/C3A045B1-F50D-C540-8AD5-2F32956D4AA1.root",
    "/data/bmc/D05EF786-4706-D045-A152-4A81F06EDB04.root",
    "/data/bmc/D20B2931-4C27-FF4E-B640-3404DFD65150.root",
    "/data/bmc/D6094FF2-1AAE-784B-9991-AA8441B9FB45.root",
    "/data/bmc/E0984EEB-ED4B-B84C-BDB7-845E8A2AB85C.root",
    "/data/bmc/E215B4B9-11EA-D84F-A6D7-9FE57515B2CE.root",
    "/data/bmc/E82ACE74-60DA-DC41-A2B0-EF750C4AB5E5.root",
    "/data/bmc/E94099CB-E3C8-EB43-B182-4EA89C7C5411.root",
    "/data/bmc/ED08B6D8-2823-D24D-B1BA-1141EA893E7B.root"

]
output_file = ROOT.TFile("hMtBackground_combinedPart2.root", "RECREATE")
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
