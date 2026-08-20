import ROOT
import time
import math
import sys, os
from analysisfxn import Analysis_up

dataPath = [
    "/data/mu/mu part1/0107961B-4308-F845-8F96-E14622BBA484.root",
    "/data/mu/mu part1/0DEE1709-0416-F24B-ACB2-C68997CB6465.root",
    "/data/mu/mu part1/0F03FF10-7F20-8A4E-A19A-418C090A3F97.root",
    "/data/mu/mu part1/0F6506E1-94CC-E14B-AE40-F72BD37BF8D4.root",
    "/data/mu/mu part1/1C08614E-0C0E-6044-966A-CAF630CAEF8F.root",
    "/data/mu/mu part1/1D87B4FB-E31C-9F43-AC21-C32469DE9FC6.root",
    "/data/mu/mu part1/1EB443F2-1230-8042-B8AE-FD50329CA59B.root",
    "/data/mu/mu part1/2045F967-9F0A-7C46-9946-787B27D56E88.root",
    "/data/mu/mu part1/236A04EE-C105-D947-8A2E-F8CC6731644F.root",
    "/data/mu/mu part1/370BE877-DA24-DB41-A875-07A86EAB6852.root",
    "/data/mu/mu part1/39B0B7DF-2C71-9E4E-952F-833475062880.root",
    "/data/mu/mu part1/3A8A8FD7-1FB0-FC49-93FD-20378419A49A.root",
    "/data/mu/mu part1/4C688018-90E1-5847-AEA3-D5F5DEB330EA.root",
    "/data/mu/mu part1/50ADF677-F507-D044-BCF2-40392F4A89DB.root",
    "/data/mu/mu part1/576759DA-4A35-534B-B926-2A9E4A5A7268.root",
    "/data/mu/mu part1/5CA4CE73-629C-2F48-939C-4274B369F112.root",
    "/data/mu/mu part1/60E957CB-2729-5D4C-B261-43406A639E8F.root",
    "/data/mu/mu part1/61FC1E38-F75C-6B44-AD19-A9894155874E.root",
    "/data/mu/mu part1/6CBB587F-1366-6F40-973C-7F98EA9C8335.root",
    "/data/mu/mu part1/715BF52A-7FE3-D848-A45F-99A0DDA80205.root",
    "/data/mu/mu part1/788CA5CB-13E0-544A-8427-33E4591CF20B.root",
    "/data/mu/mu part1/791500CE-B4AB-7C46-91D7-4E689B94F55A.root",
    "/data/mu/mu part1/8F3012A6-72CA-6C4F-81AA-64EEB640EA84.root",
    "/data/mu/mu part1/8FB46D51-A84B-8043-A86B-25743CE20170.root",
    "/data/mu/mu part1/BE433087-4E55-6F44-8EBF-EDD13E745BAE.root",
    "/data/mu/mu part1/CBFBEA23-7C46-624E-98B1-2CAEE35458C8.root",
    "/data/mu/mu part1/CEB3C22E-30A5-704F-BD07-BD6D0838DAFA.root",
    "/data/mu/mu part1/CEE8B116-DBE4-E54F-B038-A81090B2A147.root",
    "/data/mu/mu part1/D2F76766-1F43-264A-A39F-2DE68BF68752.root",
    "/data/mu/mu part1/DBA71BF9-1B93-EC42-A96C-BD2B35E3538E.root",
    "/data/mu/mu part1/E6C5C5AA-946B-A34C-ABB5-4B14AA0E9ED7.root",
    "/data/mu/mu part1/E8E0A6FB-E2D7-BE40-B098-7FEABFB340FF.root",
    "/data/mu/mu part1/EC4283CE-BB96-B845-8012-6ED0CC0FA7A2.root",
    "/data/mu/mu part1/F4101B00-7984-C54F-98E9-2EF72AD66F43.root",
    "/data/mu/mu part1/FE81EE6E-0E21-1847-9105-D700F2E762DB.root"
]


output_file = ROOT.TFile("hMTData_combinedPart1.root", "RECREATE")
hMTDataTotal = None
start = time.time()
for i in range(len(dataPath)):
    file_start = time.time()
    hMTData = Analysis_up(dataPath[i], hist_name="htemp%d" % i, is_mc = False)
    print "[%d/%d] %s done in %.1fs" % (
        i + 1, len(dataPath), dataPath[i], time.time() - file_start)

    if hMTDataTotal is None:
        hMTDataTotal = hMTData.Clone("hMTDataTotal")
    else:
        hMTDataTotal.Add(hMTData)
print "\nAll files processed in %.1fs" % (time.time() - start)


output_file.cd()
hMTDataTotal.Write()
output_file.Close()
