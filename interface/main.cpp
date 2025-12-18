#include "PidTrdFiller.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

int main(int argc, char** argv) {

  if (argc < 6) {
    std::cout << "Wrong number of arguments! Please use: " << std::endl;
    std::cout << " ./main filelist.txt outputfile truncation_mode probability_mode min_hits\n";
    return EXIT_FAILURE;
  }

  const std::string& filelist = argv[1];
  const std::string& getter_file = "pidtrd_getter.root";
  const std::string& getter_name = "PidTrdGetter";

  int truncation_mode = atoi(argv[3]);
  int prob_mode = atoi(argv[4]);
  int trdhits_min = atoi(argv[5]);
  float purity = 0.1;

  std::cout << "truncation mode : " << truncation_mode << std::endl;
  std::cout << "probability mode: " << prob_mode       << std::endl;
  std::cout << "min hits trd    : " << trdhits_min     << std::endl;

  TString output_file = argv[2];
  auto* man = TaskManager::GetInstance();
  man->SetOutputName(output_file.Data(), "aTree");

  man->SetWriteMode(eBranchWriteMode::kCopyTree);;
  man->SetBranchesExclude({"RecTracks"});

  auto* pid_trd_filler = new PidTrdFiller(getter_file, getter_name);
  pid_trd_filler->SetRecTracksName("RecTracks");
  pid_trd_filler->SetTrdTracksName("TrdTracks");
  pid_trd_filler->SetRichRingsName("RichRings");

  pid_trd_filler->SetMinHits(trdhits_min);
  pid_trd_filler->SetTruncationMode(truncation_mode);
  pid_trd_filler->SetProbabilityMode(prob_mode);
  pid_trd_filler->SetPurity(purity);

  man->AddTask(pid_trd_filler);  
  man->Init({filelist}, {"pTree"});
  man->Run(-1);// -1 = all events
  man->Finish();
  man->ClearTasks();

  return EXIT_SUCCESS;
}
