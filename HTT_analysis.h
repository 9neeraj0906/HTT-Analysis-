#ifndef HTT_ANALYSIS_H
#define HTT_ANALYSIS_H

#include <cmath>
#include "Math/VectorUtil.h"
#include "ROOT/RVec.hxx"

using namespace ROOT::VecOps;

// Transverse mass calculation
float computeTransverseMass(float muon_pt, float muon_phi, float met_pt, float met_phi) {
    // Manual dphi calculation, handling wrap-around at PI
    float dphi = muon_phi - met_phi;
    if (dphi > M_PI) dphi -= 2 * M_PI;
    if (dphi < -M_PI) dphi += 2 * M_PI;

    return std::sqrt(2 * muon_pt * met_pt * (1 - std::cos(dphi)));
}

// Muon selection function

bool hasGoodMuon(RVec<float> &pt){
    for (size_t i = 0; i < pt.size(); ++i) {
         
        if (pt[i] > 20) {
            return true;
        }
    }
    return false;
}
// Tau selection function for 2016 data using bitmask ID
bool hasGoodTau(
    RVec<float>& pt,
    RVec<float>& eta,
    RVec<unsigned char>& idMVA,
    RVec<int>& decayMode
) {
    for (size_t i = 0; i < pt.size(); ++i) {
        bool kinematics = pt[i] > 30 && std::abs(eta[i]) < 2.3;

        
        bool identification = (idMVA[i] & (1 << 4)) != 0;

        bool goodDecayMode = decayMode[i] != 5 && decayMode[i] != 6;
        
        if (kinematics && identification && goodDecayMode) {
            return true;
        }
    }
    return false;
}

#endif // HTT_ANALYSIS_H
