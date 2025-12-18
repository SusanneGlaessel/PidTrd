#include <TROOT.h>

#include "PidTrdMcInput.hpp"
#include "AnalysisTree/TaskManager.hpp"

using namespace AnalysisTree;

int create_mcinput(const std::string& filelist, const std::string& outpath, const std::string& outfilename, const std::string& jobId) {

  auto start = std::chrono::system_clock::now();
  ROOT::EnableImplicitMT(4);

  gROOT -> SetBatch (true);
  gROOT -> SetEscape (true);

  const std::string& outfile = Form("%s/%s.%s", outpath.data(), jobId.data(), outfilename.data());
  
  Bool_t update_mchistos = kFALSE;
  
  auto* pid_trd_mcinput = new PidTrdMcInput(outfile);
  pid_trd_mcinput->SetRecTracksName("RecTracks"); //analysistree with tof-pid
  //pid_trd_mcinput->SetRecTracksName("VtxTracks"); // analysistree without tof-pid
  pid_trd_mcinput->SetSimTracksName("SimParticles");
  pid_trd_mcinput->SetTrdTracksName("TrdTracks");
  pid_trd_mcinput->SetUpdateMcHistos(update_mchistos);

  auto* man = TaskManager::GetInstance();
  man->AddTask(pid_trd_mcinput);  
  man->Init({filelist}, {"pTree"}); //analysistree with tof-pid
  //man->Init({filelist}, {"rTree"});  // analysistree without tof-pid
  man->Run(-1);// -1 = all events
  man->Finish();
  man->ClearTasks();
  
  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {
  
  if (argc < 5) {
    std::cout << "Wrong number of arguments! Please use: " << std::endl;
    std::cout << " ./create_mcinput filelist.txt pathname outfilename jobId\n";
    return EXIT_FAILURE;
  }

  const std::string& filelist    = argv[1];
  const std::string& outpath     = argv[2];
  const std::string& outfilename = argv[3];
  const std::string& jobId       = argv[4];
  create_mcinput(filelist, outpath, outfilename, jobId);
  
  return 0;
}
