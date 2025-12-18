#include "Getter.hpp"
#include <iostream>
using std::cout;
using std::endl;
using std::to_string;

ClassImp(PidTrd::Getter)

namespace PidTrd {
  std::array<float, NumberOfPidsTrd> Getter::GetTrdProbabilities(TrdContainer trdtrack, int ihit) {
    std::array<float, NumberOfPidsTrd> prob;
    prob.at(NumberOfPidsTrd-1) = 1.0;

    int nhits = trdtrack.GetNhitsTrd();
    float mom = trdtrack.GetP();
    float dEdx;
    if (prob_mode_ == 0) {
      trdtrack.CalculateEnergyLossTrack(trunc_mode_);
      dEdx = trdtrack.GetdEdxTrack(trunc_mode_);
    }
    if (prob_mode_ == 1) dEdx = trdtrack.GetdEdxHits().at(ihit);
    
    for(size_t ipid = 0; ipid < NumberOfPidsTrd - 1; ipid++) {
      Int_t ipid_pm = trdtrack.GetCharge() > 0 ? ipid : ipid + NumberOfPidsTrd-1;
      Int_t trunc_mode_getter;
      if (trunc_mode_ == 0) trunc_mode_getter = trdtrack.GetNhitsTrd() - 1;
      else if (trunc_mode_ > trdtrack.GetNhitsTrd()) trunc_mode_getter = trdtrack.GetNhitsTrd() - 1; // if (truncation > nhits) truncation = truncation of nhits
      else trunc_mode_getter = trunc_mode_ - 1;
      Int_t id = ipid_pm + (trdtrack.GetNhitsTrd() - 1) * 100 + trunc_mode_getter * 1000 + prob_mode_ * 10000;
      prob.at(ipid) = particles_prob_.find(id)->second.Eval(mom,dEdx);
      prob.at(NumberOfPidsTrd-1) -= prob.at(ipid);
    }
    return prob;
  }

  std::array<float, NumberOfPidsTrd> Getter::GetTrdProbabilitiesMulti(TrdContainer trdtrack) {
    std::array<float, NumberOfPidsTrd> pid_prob_trd;

    trdtrack.SelectHitIndices(trunc_mode_);
    
    for (int ipid = 0; ipid < NumberOfPidsTrd; ipid++) 
      pid_prob_trd.at(ipid) = 1.0;
    std::array<float, NumberOfPidsTrd> pid_prob_tmp;
    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++) {
      if (trdtrack.GetHitsSelIndex().at(ihit) == false) continue;   
      pid_prob_tmp = GetTrdProbabilities(trdtrack, ihit);
      for (int ipid = 0; ipid < NumberOfPidsTrd-1; ipid++) 
	pid_prob_trd.at(ipid) *= pid_prob_tmp.at(ipid);
    }

    float prob_tot_trd = 0.0;
    for (int ipid = 0; ipid < NumberOfPidsTrd-1; ipid++) 
      if (pid_prob_trd.at(ipid) >= 0. && pid_prob_trd.at(ipid) <= 1.) prob_tot_trd += pid_prob_trd.at(ipid);
    for (int ipid = 0; ipid < NumberOfPidsTrd-1; ipid++) {
      if (prob_tot_trd > 0) pid_prob_trd.at(ipid) /= prob_tot_trd;
      else pid_prob_trd.at(ipid) = 0.0;
    }
    pid_prob_trd.at(NumberOfPidsTrd-1) = 0.0;
    return pid_prob_trd;
  }

  int Getter::GetTrdPid(std::array<float, NumberOfPidsTrd> prob, float purity, int charge) {
    int pid;
    auto prob_max = std::max_element(std::begin(prob), std::end(prob));
    auto prob_max_index = std::distance(prob.begin(), prob_max);
    if (*prob_max >= purity)
      pid = pid_codes_trd_.at(prob_max_index).first*charge;
    else
      pid = 1*charge;
    return pid;
  }
}
