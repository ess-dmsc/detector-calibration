#include <TH1D.h>
#include <TPolyMarker.h>
#include <TList.h>
#include <TF1.h>
#include <H5Cpp.h>
#include <TCanvas.h>
#include <vector>
#include <nlohmann/json.hpp>



class CalibrationCalculator {

public:
    struct Straw {
      std::vector<double> measuredPeaks;
      std::vector<double> simulatedPeaks;
      bool goodStraw;
      std::vector<double> calibrationParameters;
      int measuredHitsCount;

      // Define a to_json method for Straw
    void to_json(nlohmann::json& j) const {
        j = {
            {"measuredPeaks", measuredPeaks},
            {"simulatedPeaks", simulatedPeaks},
            {"goodStraw", goodStraw},
            {"calibrationParameters", calibrationParameters},
            {"measuredHitsCount", measuredHitsCount}
        };
    }

    // Define a from_json method for Straw
    void from_json(const nlohmann::json& j) {
        j.at("measuredPeaks").get_to(measuredPeaks);
        j.at("simulatedPeaks").get_to(simulatedPeaks);
        j.at("goodStraw").get_to(goodStraw);
        j.at("calibrationParameters").get_to(calibrationParameters);
        j.at("measuredHitsCount").get_to(measuredHitsCount);
    }
    };

  std::map<int, Straw> strawInfo;

  // whether graphs get saved to file or not
  bool plottingGraphs = true;

  // number of pixels per straw
  int strawResolution = 512;
  int nStraws = 896;
  
  // once calibration parameters are calculated, go from straw n to (resolution - n) and check that calibration parameters don't push the values higher than resolution or lower than 0
  int rangeChecking = 20;
  int minimumMeasuredHitsCount = 100;

  /// \brief takes measured and simulated events as a vector of vectors of ints, where each vector of ints represents a straw, and each int represents a neutron event's position along the straw
  /// and calculates a 4th order polynomial relationship between the measured peaks of events and the simulated peaks of events, then saving it to a json file
  void calculateCalibration(std::vector<std::vector<int>> measuredEvents, std::vector<std::vector<int>> simulatedEvents);

  /// \brief takes the hits for a single straw and uses ROOTs ShowPeaks function, doesn't do any further refinement
  std::vector<double> getStrawPeaksSimple(std::vector<int> hits, int strawNum, std::string file_prefix);
  
  /// \brief takes the hits for a single straw and uses ROOTs ShowPeaks function for first pass, and then multi gaussian fitting to refine
  /// and returns the x values of the peaks on the straw
  std::vector<double> getStrawPeaksGaussian(std::vector<int> hits, int strawNum, std::string file_prefix);
  
  /// \brief inserts each hit into the histogram
  void fillHistogram1D(std::vector<int> hits, TH1D* histogram);
  
  /// \brief given the x and y values of pre-calculated peaks, and a histogram, refines those peaks using multi-gaussian fits
  /// x will be adjusted according to the newly refined peak locations, return is void
  void gaussianFit(std::vector<double> x, std::vector<double> y, TH1D* histogram);
  
  
  std::pair<std::vector<double>, std::vector<double>> findPeaks(TH1D* histogram);
    
  /// \brief single straw calibration parameter calculation
  std::vector<double> calculateStrawCalibrationParameters(std::vector<double> measuredPeaks, std::vector<double> simulatedPeaks, int strawId);


  void writeStrawInfoToFile(std::string filename);
  void loadStrawInfoFromFile(std::string filename);
  void saveCalibrationParametersToFile();

  // sorts the elements in array a from smallest to largest, and re-organises
  // array b the same way, regardless of the numerical order of elements in b
  // n is the length of both arrays a and b 
  void selectionSort(std::vector<double> a, std::vector<double> b, int n);

  /// \brief iterated over positions from rangeChecking to strawResolution - rangeChecking, and if any are thrown outside the bounds of 0 to strawResolution, return false
  bool checkRange(std::vector<double> strawCalibrationParams);

  std::vector<double> applyCalibrationParams(std::vector<double> measured, std::vector<double> params);
  
  void to_json(nlohmann::json& j, const std::map<int, Straw>& m) {
    for (auto& [key, value] : m) {
        j[std::to_string(key)] = value;
    }
  }

};