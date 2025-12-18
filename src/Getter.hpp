#ifndef PIDTRD_GETTER_HPP_
#define PIDTRD_GETTER_HPP_

#include "TObject.h"
#include "TString.h"
#include "ParticleProb.hpp"
#include "TrdContainer.hpp"

namespace PidTrd {
  class Getter : public TObject {
    
  public:
    Getter() = default;
    virtual ~ Getter() = default;

    void SetMinHits(int nhits_min) { nhits_min_ = nhits_min; } 
    void SetTruncationMode(int trunc_mode) { trunc_mode_ = trunc_mode; }  
    void SetProbabiltyMode(int prob_mode) { prob_mode_ = prob_mode; }
   
    int GetMinHits() { return nhits_min_; } 

    void AddParticlesProb(std::map<int, ParticleProb> particlesprob) {particles_prob_ = particlesprob;}
    void AddParticleProb(ParticleProb particleprob) {
      int type = particleprob.GetCharge() > 0 ? particleprob.GetType() : particleprob.GetType() + NumberOfPidsTrd-1;
      int nhits = particleprob.GetNhits();
      int truncmode = particleprob.GetTruncMode();
      int probmode = particleprob.GetProbMode();  
      particles_prob_[type + nhits * 100 + truncmode * 1000 + probmode * 10000] = particleprob;
    }
    
    std::map<int, ParticleProb> GetParticlesProb() { return particles_prob_; }
    ParticleProb GetParticleProb(int id) { return particles_prob_[id];}
    ParticleProb GetParticleProb(int type, int charge, int nhits, int truncmode, int probmode) {
      if (charge < 0) type += NumberOfPidsTrd-1;
      return particles_prob_[type + nhits * 100 + truncmode * 1000 + probmode * 10000];
    }
    
    std::array<float, NumberOfPidsTrd> GetTrdProbabilities(TrdContainer trdtrack, int ihit = -1);
    std::array<float, NumberOfPidsTrd> GetTrdProbabilitiesMulti(TrdContainer trdtrack);
    int GetTrdPid(std::array<float, NumberOfPidsTrd> prob, float purity, int charge);
    
  private:

    std::map<int, ParticleProb> particles_prob_{};
    
    int nhits_min_{1};
    int trunc_mode_{0};
    int prob_mode_{0};

    array<TString, NumberOfProbMode> probnames_ = {"probP", "probE"};
    
    ClassDef(Getter, 2);
  }; 
}// namespace PidTrd
#endif//PIDTRD_GETTER_HPP_
