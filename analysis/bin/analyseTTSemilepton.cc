#include "delphes/analysis/interface/TTSemileptonAnalyser.h"

#include "TSystem.h"

int main(int argc, char* argv[]) {
  if (argc != 6) {
    // FIXME
    std::cerr << "Usage: analyseTTSemilepton isTT chargeAlgorithm exponent outpath inputFile" << std::endl;
    return 1;
  }

  bool is_tt = std::atoi(argv[1]);
  std::string chargeAlgorithm = argv[2];
  float exponent = std::stof(argv[3]);
  TString out_path(argv[4]);
  TString in_path(argv[5]);

  //Ssiz_t len=in_path.Length();
  //Ssiz_t start=in_path.Index("delphes342");

  //TSubString sub=in_path(start, len-start);
  //TString outfile(sub);

  //TString output = out_path+outfile;

  TTSemileptonAnalyser analyser(in_path, out_path, is_tt, chargeAlgorithm, exponent);
  analyser.loop();

  return 0;
}

