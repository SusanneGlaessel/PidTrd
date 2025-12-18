#include <TROOT.h>
#include <iostream>
#include "PidTrdRunGetter.hpp"
#include "AnalysisTree/TaskManager.hpp"
#include "TStopwatch.h"

int create_getter(const std::string& path, const std::string& mcfile_name, const std::string& jobId) {

  Bool_t write_mchistos_out = kTRUE;

  TString macroname = "create_getter";
  TString inFile = Form("%s/%s.%s", path.data(), jobId.data(), mcfile_name.data());
  const std::string& getter_file = "pidtrd_getter.root";
  const std::string& getter_name = "PidTrdGetter";
  TString outFile;
  if (write_mchistos_out == kTRUE) outFile = Form("%s/%s.%s_probabilities.root", path.data(), jobId.data(), mcfile_name.data());
  
  std::cout << std::endl;
  std::cout << "-I- " << macroname << ": Using input file " << inFile << std::endl;

  TStopwatch timer;
  timer.Start();
  
  auto start = std::chrono::system_clock::now(); 
  ROOT::EnableImplicitMT(4);

  gROOT -> SetBatch (true);
  gROOT -> SetEscape (true);

  auto* pid_trd_rungetter = new PidTrdRunGetter(path, mcfile_name, jobId, getter_file, getter_name);
  pid_trd_rungetter->SetWriteMcHistogramsOut(write_mchistos_out);
  pid_trd_rungetter->Init();
  pid_trd_rungetter->Exec();
  pid_trd_rungetter->Finish();

  timer.Stop();
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished successfully." << std::endl;
  std::cout << "Getter is " << getter_file << std::endl;
  if (write_mchistos_out == kTRUE) std::cout << "Output file is " << outFile << std::endl;
  std::cout << "Real time " << rtime << " s, CPU time " << ctime << "s"
            << std::endl
            << std::endl;

  return EXIT_SUCCESS;
}

int main(int argc, char **argv) {

  
  if (argc < 4) {
    std::cout << "Wrong number of arguments! Please use: " << std::endl;
    std::cout << " ./create_getter path_mc filename_mc jobId\n";
    return EXIT_FAILURE;
  }

  const std::string& path = argv[1];
  const std::string& mcfile_name = argv[2];
  const std::string& jobId = argv[3];
  create_getter(path, mcfile_name, jobId);
  
  return 0;
}
