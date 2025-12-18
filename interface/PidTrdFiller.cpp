#include "PidTrdFiller.hpp"

#include "AnalysisTree/Task.hpp"
#include "AnalysisTree/TaskManager.hpp"

using std::string;
using std::cout;
using std::endl;
using namespace AnalysisTree;

PidTrdFiller::PidTrdFiller(const std::string& getter_file, TString getter_name) {
  
  std::unique_ptr<TFile> pid_file( TFile::Open(getter_file.c_str()) );

  if ((!pid_file) || (pid_file->IsZombie())) {
    throw std::runtime_error("No file or file is zombie: " + getter_file);
  }

  getter_ = std::unique_ptr<PidTrd::Getter>(pid_file->Get<PidTrd::Getter>(getter_name));
    
}

float PidTrdFiller::GetMomentum(const AnalysisTree::BranchChannel& trd_track) {
  float momentum = 0.0; 
  if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("p")) > 0))
    momentum = trd_track.Value(trd_tracks_.GetField("p"));
  else if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("p_out")) > 0))
    momentum = trd_track.Value(trd_tracks_.GetField("p_out"));
  else
    Warning("Exec", "Could not assign any momentum to the track, use p=0.");
    
  return momentum;
}

float PidTrdFiller::GetPt(const AnalysisTree::BranchChannel& trd_track) { 
  float pT = 0; 
  if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("pT")) > 0))
    pT = trd_track.Value(trd_tracks_.GetField("pT"));
  else if (TMath::Abs(trd_track.Value(trd_tracks_.GetField("pT_out")) > 0))
    pT = trd_track.Value(trd_tracks_.GetField("pT_out"));
  else 
    Warning("Exec", "Could not assign any pT to the track, use pT=0.");
  return pT;
}

float PidTrdFiller::GetPz(const AnalysisTree::BranchChannel& trd_track, float mom = -1.) {

  float pT = GetPt(trd_track);
  if (mom == -1.){
    mom = GetMomentum(trd_track);

  }
  float pz = 0.0;
  if (TMath::Abs(mom) > TMath::Abs(pT))
    pz = TMath::Sqrt(mom*mom - pT*pT);
  else 
    Warning("Exec", "Could not assign pz to the track, use pz=0.");

  return pz;
}

int PidTrdFiller::GetCharge(const AnalysisTree::BranchChannel& rec_track) {
  return rec_track.Value(rec_tracks_.GetField("q"));
}

void PidTrdFiller::GetEnergyLossHits(const AnalysisTree::BranchChannel& trd_track, int &nhits_trd, std::array<float, NumberOfTrdLayers> &dEdx) {

  float mom = GetMomentum(trd_track);
  float pz = GetPz(trd_track, mom);
  
  nhits_trd = 0;
  
  for (int ihit = 0 ; ihit < NumberOfTrdLayers; ihit++) {
    if ( trd_track.Value(trd_dEdx_field_.at(ihit)) > 0) {
      dEdx.at(ihit) = trd_track.Value(trd_dEdx_field_.at(ihit));   
      if ( pz > 0 && mom > 0) dEdx.at(ihit) *= pz / mom;
      nhits_trd ++;
    }
    else
      dEdx.at(ihit) = 0;  
  }
}

bool PidTrdFiller::IsRichElectron(const AnalysisTree::BranchChannel& rec_track) {
  bool isElectron = false;
  int i_rich = rec_to_rich_->GetMatch(rec_track.GetId());
  if (i_rich > 0) {
    const auto& rich_ring = rich_rings_[i_rich];
    float Aaxis = rich_ring.Value(rich_rings_.GetField("axis_a"));
    float Baxis = rich_ring.Value(rich_rings_.GetField("axis_b"));
    Double_t dist  = 0;  // richRing->GetDistance();
	    
    Float_t mom = GetMomentum(rec_track);
    Double_t MeanA    = 4.95;
    Double_t MeanB    = 4.54;
    Double_t RmsA     = 0.30;
    Double_t RmsB     = 0.22;
    Double_t RmsCoeff = 3.5;
    Double_t DistCut  = 1.;
	    
    if (mom < 5.) {
      if (fabs(Aaxis - MeanA) < RmsCoeff * RmsA && fabs(Baxis - MeanB) < RmsCoeff * RmsB && dist < DistCut)
	isElectron = true;
    }
    else {                     
      ///2 sigma
      Double_t polAaxis = 5.64791 - 4.24077 / (mom - 3.65494);
      Double_t polBaxis = 5.41106 - 4.49902 / (mom - 3.52450);
      if (Aaxis < (MeanA + RmsCoeff * RmsA) && Aaxis > polAaxis && Baxis < (MeanB + RmsCoeff * RmsB) && Baxis > polBaxis && dist < DistCut)              
	isElectron = true;
    }      
  }
  return isElectron;
}

void PidTrdFiller::Init() {

  getter_->SetMinHits(nhits_min_);
  getter_->SetTruncationMode(trunc_mode_);
  getter_->SetProbabiltyMode(prob_mode_);

  auto man = TaskManager::GetInstance();
  auto chain = man->GetChain();

  chain->InitPointersToBranches({});

  rec_tracks_ = chain->GetBranchObject(rec_tracks_name_);
  trd_tracks_ = chain->GetBranchObject(trd_tracks_name_);
  rich_rings_ = chain->GetBranchObject(rich_rings_name_);
  rec_to_trd_ = chain->GetMatching(rec_tracks_name_, trd_tracks_name_);
  rec_to_rich_ = chain->GetMatching(rec_tracks_name_, rich_rings_name_);
  
  in_branches_.emplace(rec_tracks_name_);
  in_branches_.emplace(trd_tracks_name_);
  in_branches_.emplace(rich_rings_name_);

  auto conf = rec_tracks_.GetConfig().Clone(out_branch_name_, AnalysisTree::DetType::kParticle);
  
  conf.AddField<bool>("electron_rich", "rich electron hypothesis");
  conf.AddField<int>("pid_trd", "trd pid hypothesis");
  conf.AddField<int>("nhits_trd_eloss", "number of trd hits dEdx > 0");

  std::vector<std::string> names{};
  for (const auto& pid : pid_codes_trd_) {
    string name = Form("prob_trd_%s", pid.second.Data());
    names.push_back(Form("prob_trd_%s", pid.second.Data()));
  }
  conf.AddFields<float>(names, "trd probability to be proton, pion, kaon etc");

  ana_tracks_ = Branch(conf);
  ana_tracks_.SetMutable();
  ana_tracks_.Freeze();

  rec_tracks_.Freeze();

  man->AddBranch(&ana_tracks_);

  int imatch{0};
  auto match_br = {"SimParticles", "RichRings", "TofHits", "TrdTracks"};
  out_matches_.assign(match_br.size(), nullptr);

  for (const auto& br : match_br) {
    in_matches_.emplace_back(chain->GetMatchPointers().find({rec_tracks_name_ + "2" + br})->second);
    man->AddMatching(out_branch_name_, br, out_matches_.at(imatch));
    imatch++;
  }

  nhits_trd_eloss_field_ = ana_tracks_.GetField("nhits_trd_eloss");
  pid_trd_field_ = ana_tracks_.GetField("pid_trd");
  el_rich_field_ = ana_tracks_.GetField("electron_rich");
  
  for (const auto& pid : pid_codes_trd_) {
    pid_prob_trd_field_.push_back(ana_tracks_.GetField(Form("prob_trd_%s", pid.second.Data())));
  }
  
  for (int i = 0; i < NumberOfTrdLayers; i++)
    trd_dEdx_field_.push_back(trd_tracks_.GetField(("energy_loss_" + std::to_string(i)).c_str()));

}

void PidTrdFiller::Exec() {
  ana_tracks_.ClearChannels();
  for (int i_track = 0; i_track < rec_tracks_.size(); ++i_track) {
    const auto& rec_track = rec_tracks_[i_track];
    auto track_new = ana_tracks_.NewChannel();
    track_new.CopyContent(rec_track);
    int itrd = rec_to_trd_->GetMatch(rec_track.GetId());
    int nhits_trd = 0;
    
    if (itrd > -1) {
      const auto& trd_track = trd_tracks_[itrd];
      float mom = GetMomentum(trd_track);
      float pT = GetPt(trd_track);
      int charge = GetCharge(rec_track);
      
      std::array<float, NumberOfTrdLayers> dEdx_hits = {0.0, 0.0, 0.0, 0.0};
      GetEnergyLossHits(trd_track, nhits_trd, dEdx_hits);

      TrdContainer trdtrack(mom, pT, charge, nhits_trd, dEdx_hits);
      trdtrack.ScaleEnergyLossLength();     
      
      int pidtrd_hypo;
      std::array<float, NumberOfPidsTrd> pid_prob_trd;
      if (nhits_trd >= nhits_min_) {
	if (prob_mode_ == 0) 
	  pid_prob_trd = getter_->GetTrdProbabilities(trdtrack);
	if (prob_mode_ == 1) 
	  pid_prob_trd = getter_->GetTrdProbabilitiesMulti(trdtrack);

	for (int ipid = 0; ipid < NumberOfPidsTrd; ipid++) 
	  track_new.SetValue(pid_prob_trd_field_.at(ipid), pid_prob_trd.at(ipid));	
	
	pidtrd_hypo = getter_->GetTrdPid(pid_prob_trd, purity_, charge);
	track_new.SetValue(pid_trd_field_, pidtrd_hypo);
      }
      else { 
	track_new.SetValue(pid_trd_field_, -2);
	for (int ipid = 0; ipid < NumberOfPidsTrd; ipid++) 
	  track_new.SetValue(pid_prob_trd_field_.at(ipid), -1.f);	
      }
    }
    else { 
      track_new.SetValue(pid_trd_field_, -2);
      for (int ipid = 0; ipid < NumberOfPidsTrd; ipid++) 
	track_new.SetValue(pid_prob_trd_field_.at(ipid), -1.f);	
    }
    
    track_new.SetValue(nhits_trd_eloss_field_, nhits_trd);

    bool isElectron = IsRichElectron(rec_track);
    track_new.SetValue(el_rich_field_, isElectron);
  }
  int i{0};
  for (auto& match : out_matches_) {
    auto m1 = in_matches_[i]->GetMatches(false);
    auto m2 = in_matches_[i]->GetMatches(true);
    match->SetMatches(m1, m2);
    i++;
  }
}

void PidTrdFiller::Finish() {
  //getter_->Finish();
}

