#ifndef TrdContainer_HPP
#define TrdContainer_HPP

#include "Constants.hpp"

#include "TMath.h"
using namespace std;
#include <stdexcept>

class TrdContainer {
 public:
  TrdContainer() = default;
  
  explicit TrdContainer(float mom, float pT, int charge, int nhits_trd, std::array<float,NumberOfTrdLayers> dEdx_hits, std::array<float, NumberOfTruncMode> dEdx_track, int nhits_sel, std::array<bool, NumberOfTrdLayers> hits_sel_index, int mc_pdg, bool dEdxIsScaled = false) : mom_(mom), pT_(pT), charge_(charge), nhits_trd_(nhits_trd), dEdx_hits_(dEdx_hits), dEdx_track_(dEdx_track), nhits_sel_(nhits_sel), hits_sel_index_ (hits_sel_index), mc_pdg_(mc_pdg), dEdx_is_scaled_(dEdxIsScaled) {}
  
  explicit TrdContainer(const float mom, float pT, int charge, int nhits_trd, std::array<float,NumberOfTrdLayers> dEdx_hits, int mc_pdg, bool dEdxIsScaled = false) : mom_(mom), pT_(pT), charge_(charge), nhits_trd_(nhits_trd), dEdx_hits_(dEdx_hits), mc_pdg_(mc_pdg), dEdx_is_scaled_(dEdxIsScaled) {
    nhits_sel_ = 0;
    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
      hits_sel_index_.at(ihit) = false;
  }
  
  explicit TrdContainer(float mom, float pT, int charge, int nhits_trd, std::array<float,NumberOfTrdLayers> dEdx_hits, bool dEdxIsScaled = false) : mom_(mom), pT_(pT), charge_(charge), nhits_trd_(nhits_trd), dEdx_hits_(dEdx_hits), dEdx_is_scaled_(dEdxIsScaled) {
    nhits_sel_ = 0;
    for (int ihit = 0; ihit < NumberOfTrdLayers; ihit++)
      hits_sel_index_.at(ihit) = false;
    mc_pdg_ = -2;
  }

  virtual ~TrdContainer() = default;

  void ScaleEnergyLossLength();
  void CalculateEnergyLossTrack(int trunc_mode);
  void CalculateEnergyLossTrackAllModes();
  void SelectHitIndices(int trunc_mode);    // Returns hit indices with dEdx > 0 selected in truncation mode
  int  GetNHitsSel(int trunc_mode);         // Returns number of hits with dEdx > 0 selected in truncation mode

  //void SetCharge(int charge) { charge_ = charge; }; //not used 
  //void SetNHitsSel ( int nhits_sel) { nhits_sel_ = nhits_sel; }; //not used
  //void SetHitsSelIndex ( std::array<bool, NumberOfTrdLayers> hits_sel_index ) { hits_sel_index_ = hits_sel_index; }; //not used
  
  float GetP()       const { return mom_;       }
  int GetCharge()    const { return charge_;    }
  int GetNhitsTrd()  const { return nhits_trd_; }
  int GetMcPdg()     const {
    if (mc_pdg_ == -2) throw std::runtime_error("MC PDG not set.");
    else return mc_pdg_;}
  
  std::array<float,NumberOfTrdLayers> GetdEdxHits() const { return dEdx_hits_; }
  float GetdEdxTrack(int trunc_mode) const { return dEdx_track_.at(trunc_mode); }

  std::array<float, NumberOfTrdLayers> GetdEdxHitsSorted (); // Orders hits from lowest to highest dEdx (hits with dEdx = 0 in last posisition)
  std::array<bool, NumberOfTrdLayers> GetHitsSelIndex ()   const {
    //if ( IndicesIsSet_  == false)
    // Warning("Exec", "Could not assign pz to the track, use pz=0. Energy loss will not be scaled.");
    //else
      return hits_sel_index_; }
  
protected:
  float mom_{0.0};
  float pT_{0.0};;
  int charge_{0};
  int nhits_trd_{0};
  int mc_pdg_{-1};

  std::array<float,NumberOfTrdLayers> dEdx_hits_ = {0.0, 0.0, 0.0, 0.0};
  std::array<float, NumberOfTruncMode> dEdx_track_ = {0.0, 0.0, 0.0, 0.0, 0.0};

  std::array<float, NumberOfTrdLayers> hits_sorted_ = {0.0, 0.0, 0.0, 0.0};
  int nhits_sel_{0}; 
  std::array<bool, NumberOfTrdLayers> hits_sel_index_= {false, false, false, false};

  bool dEdx_is_scaled_{false};
  bool indices_are_sel_{false};
  
};

#endif//TrdContainer_HPP
