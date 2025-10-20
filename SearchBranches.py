import ROOT
# Open ROOT file
dataPath = "/data/mu/0DEE1709-0416-F24B-ACB2-C68997CB6465.root"

f = ROOT.TFile.Open(dataPath, "READ")

# Load the tree
tree = f.Get("Events")

# --- Verify file and tree ---
print("File opened:", f.GetName())   # File name
print("Tree loaded:", tree.GetName())  # Tree name

for b in tree.GetListOfBranches():
	if "HLT" in b.GetName():
		print(b.GetName())
