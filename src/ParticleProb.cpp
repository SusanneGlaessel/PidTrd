#include "ParticleProb.hpp"

#include <iostream>
using std::cout;
using std::endl;

ClassImp(PidTrd::ParticleProb)

namespace PidTrd {
    
  float ParticleProb::Eval(float mom, float dEdx)
  {
    Int_t binx = hprobabilities_->GetXaxis()->FindBin(mom);
    Int_t biny = hprobabilities_->GetYaxis()->FindBin(dEdx);
    return hprobabilities_->GetBinContent(binx,biny);
  }
}
    
