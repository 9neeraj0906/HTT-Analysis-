import ROOT
import math
import time

from analFxn import Analysis

dataPath = [

    "/data/mu/01479010-3B51-B04A-9C5C-7800085D11B8.root",
    "/data/mu/02C01469-552A-7E45-9894-6406F0030DB3.root",
    "/data/mu/120B99EA-7A22-094A-889D-E9049DB1A21C.root",
    "/data/mu/136621ED-C32B-B347-9FC1-1BE4AF4F8EF3.root",
    "/data/mu/13EF306B-E9EA-B24F-8094-DAC8A7BC4AE0.root",
    "/data/mu/148777C2-F638-714C-AF76-00A32CD6F32A.root",
    "/data/mu/14EA6AFB-593D-074F-AF31-BF588987BE99.root",
    "/data/mu/17F6E03B-6B69-4341-863B-538EED26B89B.root",
    "/data/mu/23953892-6295-B341-9D0B-6892EBDD67FC.root",
    "/data/mu/28415B3F-DFB1-5542-9217-7D134447AC48.root",
    "/data/mu/2AF26192-4B3A-324D-B14C-CFFFB1F69F4C.root",
    "/data/mu/31BF084F-D414-1542-B61D-0F58520A17E8.root",
    "/data/mu/388AB3E1-8708-7D42-91BA-83E52373E808.root",
    "/data/mu/3BCDC9CF-57F5-104F-8C20-2BA300439111.root",
    "/data/mu/3E76909B-34D1-DE46-8355-51C2ADD587D5.root",
    "/data/mu/40EA32B0-1284-7246-83A3-A703D908F9FC.root",
    "/data/mu/4AD09F91-D53E-E345-817A-BBC670875478.root",
    "/data/mu/4B985344-7CD7-6142-BCB6-B95E47279B9A.root",
    "/data/mu/579013D1-626D-D943-A8B0-A1A558C54F33.root",
    "/data/mu/594F3D70-ED17-3D42-A322-269D0A074F9A.root",
    "/data/mu/626C68F8-9DBD-194B-993E-8672B8E73BAE.root",
    "/data/mu/6322429F-B620-854D-B4AC-6892359CAB2C.root",
    "/data/mu/6351C999-8EF7-6345-822D-E2EF1C203603.root",
    "/data/mu/6C8BBA90-79DC-CA45-9E0F-4AB683E91D8F.root",
    "/data/mu/7BA44B30-9EBA-6E4A-8DD6-2A239BDDBA50.root",
    "/data/mu/84E83B69-73A2-BF4C-A59C-C3C10F71C76E.root",
    "/data/mu/874B6634-10CD-4541-97E6-7A97BA44F4C6.root",
    "/data/mu/8E33310F-1983-EE41-A97A-7D6853AF1BAD.root",
    "/data/mu/8F897B77-7EDA-B54B-812C-2185345CE97F.root",
    "/data/mu/977A32B2-1672-0646-8029-C12B37014B64.root",
    "/data/mu/9C99B9F7-7F24-FE4A-A28D-E31D70DBDF66.root",
    "/data/mu/A5B416ED-2D1C-4846-92FF-CA7E5C6DBA61.root",
    "/data/mu/B0CFD000-E84C-344D-BED7-8AB35F3DBA36.root",
    "/data/mu/B4E0C1C3-E86E-3F4B-9B1A-D17E3ADB0908.root",
    "/data/mu/C250DA77-36D8-2440-8F71-4A913F18EA6C.root",
    "/data/mu/C2A7D179-B263-704C-A8EB-B75F4D82F7B5.root",
    "/data/mu/C653EBCD-C75E-8441-BB1D-12A90F76D597.root",
    "/data/mu/C8653579-E6F0-5D44-ACEA-154A482DF4AD.root",
    "/data/mu/CD1172F5-A00E-0041-9A9C-5E1862026357.root",
    "/data/mu/D548310D-36A5-DE46-97C5-15933EAF7DC9.root",
    "/data/mu/D87A4CC6-6CF7-6943-838E-475A7533EAA9.root",
    "/data/mu/E037AB14-FEA9-4D40-AA30-9ADF56BD6270.root",
    "/data/mu/F250A551-A8DC-BC47-895C-3BE19459F5DC.root",
    "/data/mu/F302A865-17B0-064B-8154-41526BB38244.root",
    "/data/mu/F773638E-CE20-944A-A120-755C3E071B2F.root",
    "/data/mu/FBD685D1-F6C9-8946-BA61-14F8DB0B1B82.root",
    "/data/mu/FFED1DB7-3626-A04C-8438-BF055BD16672.root"
]


output_file = ROOT.TFile("hMTData_combinedPart1.root", "RECREATE")
hMTDataTotal = None
for i in range(len(dataPath)):
    hMTData = Analysis(dataPath[i], hist_name="htemp%d" % i, is_mc = False)
    if hMTDataTotal is None:
        hMTDataTotal = hMTData.Clone("hMTDataTotal")
    else:
        hMTDataTotal.Add(hMTData)


output_file.cd()
hMTDataTotal.Write()
output_file.Close()
