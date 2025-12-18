#ifndef PIDTRD_INTERFACE_PIDFILLER_HPP_
#define PIDTRD_INTERFACE_PIDFILLER_HPP_

#include "Constants.hpp"
#include "AnalysisTree/Task.hpp"
#include "AnalysisTree/TaskManager.hpp"

#include "TH2F.h"

using std::make_pair;

#include "Getter.hpp"
#include "TrdContainer.hpp"

class PidTrdFiller : public AnalysisTree::Task {

 public:

  PidTrdFiller(const std::string& getter_file, TString getter_name);
    ~PidTrdFiller() override = default;

  void Init() override;
  void Exec() override;
  void Finish() override;
  /*void Finish() override {
    auto* man = AnalysisTree::TaskManager::GetInstance();
    }*/

  void SetRecTracksName(const std::string& name) { rec_tracks_name_ = name; }
  void SetTrdTracksName(const std::string& name) { trd_tracks_name_ = name; };
  void SetRichRingsName(const std::string& name) { rich_rings_name_ = name; };

  void SetMinHits(int nhits_min) { nhits_min_ = nhits_min; }            // Min. number of hits per track
  void SetTruncationMode(int trunc_mode) { trunc_mode_ = trunc_mode; }  // Calculation of energy loss for up to 4 layers:
                                                                        // =0: <dEdx> average over all hits
                                                                        // =1-4: Select hits with lowest dEdx:
                                                                        // =1: 1 hit, =2: 2 hits, =3: 3 hits, =4: 4 hits
  void SetProbabilityMode(int prob_mode) { prob_mode_ = prob_mode; }    // Probability for particle species i:
                                                                        // =0: total probability - probability based on particle multiplicites i
									// =1: likelihood - probability based on dEdx-distribution of particle
  
  void SetPurity(const float purity) { purity_ = purity; }

 protected:
  
  //void InitTrdMc();
  float GetMomentum(const AnalysisTree::BranchChannel& trd_particle);
  float GetPt(const AnalysisTree::BranchChannel& trd_particle);
  float GetPz(const AnalysisTree::BranchChannel& trd_particle, float mom);
  int GetCharge(const AnalysisTree::BranchChannel& trd_particle);
  void GetEnergyLossHits(const AnalysisTree::BranchChannel& trd_track, int &nhits_trd, std::array<float, NumberOfTrdLayers> &dEdx);
  bool IsRichElectron(const AnalysisTree::BranchChannel& rec_particle);
  
  AnalysisTree::Branch rec_tracks_;
  AnalysisTree::Branch trd_tracks_;
  AnalysisTree::Branch rich_rings_;
  AnalysisTree::Branch ana_tracks_;
  
  AnalysisTree::Matching* rec_to_trd_{nullptr};
  AnalysisTree::Matching* rec_to_rich_{nullptr};

  std::vector<AnalysisTree::Field> trd_dEdx_field_{};
  
  AnalysisTree::Field el_rich_field_;
  AnalysisTree::Field pid_trd_field_;
  AnalysisTree::Field nhits_trd_eloss_field_;
  std::vector<AnalysisTree::Field> pid_prob_trd_field_{};

  std::vector<AnalysisTree::Matching*> in_matches_{};
  std::vector<AnalysisTree::Matching*> out_matches_{};

  std::string rec_tracks_name_{"RecTracks"};    // Branch with input tracks
  std::string trd_tracks_name_{"TrdTracks"};
  std::string rich_rings_name_{"RichRings"}; 
  std::string out_branch_name_{"RecTracks"};
  
  TString pidtrd_path_;
  TFile *outFileTrd_;
  TFile *inFileTrd_;
  TH2F* h2dEdx_p_pos_[NumberOfTrdLayers];
  TH2F* h2dEdx_p_neg_[NumberOfTrdLayers];
  TH2F* h2dEdx_p_pdg_prob_[4][(NumberOfPidsTrd-1)*2];
  
  // field ids for input parameters

  int trunc_mode_{0};
  int prob_mode_{0};
  float purity_{0.5};
  int nhits_min_{1};

  std::vector<TrdContainer> trd_container_;
  
  std::unique_ptr<PidTrd::Getter> getter_{};

};

#endif//PIDTRD_AT_INTERFACE_PIDFILLER_HPP_
