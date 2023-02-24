#include <loki/CalibrationCalculator.h>
#include <iostream>
#include <fstream> 
#include <TGraph.h>

void CalibrationCalculator::calculateCalibration(std::vector<std::vector<int>> measuredEvents, std::vector<std::vector<int>> simulatedEvents){
  std::map<int, std::vector<double>> calibrationPerStraw;
  for (int strawNumber = 0; strawNumber < measuredEvents.size(); strawNumber++){
  // for (int strawNumber = 0; strawNumber < 1; strawNumber++){
    strawInfo[strawNumber].measuredHitsCount = measuredEvents.size();
    if(strawInfo[strawNumber].measuredHitsCount >= minimumMeasuredHitsCount){
      strawInfo[strawNumber].measuredPeaks = getStrawPeaksGaussian(measuredEvents[strawNumber], strawNumber, "measured");
      strawInfo[strawNumber].simulatedPeaks = getStrawPeaksGaussian(measuredEvents[strawNumber], strawNumber, "measured");
      

      if (strawInfo[strawNumber].measuredPeaks.size() +1 != strawInfo[strawNumber].simulatedPeaks.size()){
        std::cout << "npeaks doesn't match, measured = " + std::to_string(strawInfo[strawNumber].measuredPeaks.size()) + " and simulated = " + std::to_string(strawInfo[strawNumber].simulatedPeaks.size()) << std::endl;
        //TODO, determine way to filter peaks so that they match
      }
      else{
        std::cout << "npeaks match!!! " + std::to_string(strawInfo[strawNumber].measuredPeaks.size())<< std::endl;
        strawInfo[strawNumber].calibrationParameters = calculateStrawCalibrationParameters(strawInfo[strawNumber].measuredPeaks, strawInfo[strawNumber].simulatedPeaks, strawNumber);
      }
    }
    else{

    }
    
  }
 saveCalibrationParametersToFile();

}

void CalibrationCalculator::saveCalibrationParametersToFile(){
  using json = nlohmann::json;

  // Create a JSON object
  json j;
  j["CaenCalibration"]["ntubes"] = 128;
  j["CaenCalibration"]["nstraws"] = 7;
  j["CaenCalibration"]["resolution"] = 512;
  
  std::vector<std::vector<double>> polynomials;
  // Create an array of polynomials
  for(double i = 0.0; i < 128*7; i++){
    if (strawInfo.find(i) == strawInfo.end()){
      std::vector<double> cal = {i, 0.0, 0.0, 0.0, 0.0};
      polynomials.push_back(cal);
    }
    else{
      std::vector<double> cal = {i};
      cal.insert(cal.end(), strawInfo[i].calibrationParameters.begin(), strawInfo[i].calibrationParameters.end());
      polynomials.push_back(cal);
    }
  }

  j["CaenCalibration"]["polynomials"] = polynomials;

  // Write the JSON data to a file
  std::ofstream file("calibration.json");
  file << j.dump(4) << std::endl;
  file.close();

}



std::vector<double> CalibrationCalculator::calculateStrawCalibrationParameters(std::vector<double> measuredPeaks, std::vector<double> simulatedPeaks, int strawId){
  std::vector<double> result;
  const int n = measuredPeaks.size();
  std::sort(measuredPeaks.begin(), measuredPeaks.end());
  std::sort(simulatedPeaks.begin(), simulatedPeaks.end());
  
  // specific to larmor 
  simulatedPeaks.erase(simulatedPeaks.begin()); 


  for (int i = 0; i < n; i++){
    std::cout << std::to_string(measuredPeaks[i]) << " " << std::to_string(simulatedPeaks[i]) << std::endl;
    simulatedPeaks[i] = measuredPeaks[i] - simulatedPeaks[i];
  }

  TGraph *graph = new TGraph(n, measuredPeaks.data(), simulatedPeaks.data());
  TF1 *function = new TF1("f", "[0] + [1]*x + [2]*x^2+ [3]*x^3", 0, strawResolution);
  graph->Fit(function, "Q");

  // getting each of the 4 polynomial values
  for (int i = 0; i < 4; ++i) {
    result.push_back(function->GetParameter(i));
    std::cout << std::to_string(function->GetParameter(i)) << std::endl;
  }

  if (plottingGraphs){
    TCanvas *canvas1 = new TCanvas("c1", "c1", 600, 400);
    graph->Draw("AP");
    std::string filename  = std::to_string(strawId) + "_calibrationfit.png";
    canvas1->SaveAs(filename.c_str());
  }

  return result;
}



std::pair<std::vector<double>, std::vector<double>> CalibrationCalculator::findPeaks(TH1D* histogram){
	TH1D* histogram_background = (TH1D*) histogram->ShowBackground(20, "same");   //Measure Background
	histogram->Add(histogram_background, -1.);                      //subtract Background
	int npeaks=histogram->ShowPeaks(2.0, "", 0.5);                           // find peaks
	std::vector<double> x;
  std::vector<double> y;

  //code to display peaks found by ShowPeaks function above
  TList *histogram_functions = histogram->GetListOfFunctions();
  TPolyMarker *amy_peaks = (TPolyMarker *)histogram_functions->FindObject("TPolyMarker");
  if (amy_peaks == nullptr) {
    std::cout << "Error getting TPolyMarker object!" << std::endl;
    return std::make_pair(x, y);
  }
  double *ax = amy_peaks->GetX();
  double *ay = amy_peaks->GetY();

   for (int a = 0; a < npeaks; a++) {
    x.push_back(ax[a]);
    y.push_back(ay[a]);
  }
	selectionSort(x, y, npeaks);
	return std::make_pair(x, y);
}

std::vector<double> CalibrationCalculator::getStrawPeaksSimple(std::vector<int> hits, int strawNum, std::string file_prefix){
  if(hits.size() == 0){
     std::vector<double> emptyVector;
    return emptyVector;
  }
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  std::pair<std::vector<double>, std::vector<double>> x_y_peaks;
  TH1D* histogram = new TH1D(("histogram_straw_" + std::to_string(strawNum)).c_str(), ("histogram_title" + std::to_string(strawNum)).c_str(), strawResolution, 0, strawResolution);
  fillHistogram1D(hits, histogram);
  x_y_peaks = findPeaks(histogram);
  histogram->Draw();
  if(plottingGraphs){
    canvas->Update();
    std::string filename = file_prefix + "_canvas_" + std::to_string(strawNum) + ".png";
    canvas->SaveAs(filename.c_str());
  }
  delete canvas;
  delete histogram;
  return x_y_peaks.first;
}

std::vector<double> CalibrationCalculator::getStrawPeaksGaussian(std::vector<int> hits, int strawNum, std::string file_prefix){
  if(hits.size() == 0){
     std::vector<double> emptyVector;
    return emptyVector;
  }
  TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 600);
  TH1D* histogram = new TH1D(("histogram_straw_" + std::to_string(strawNum)).c_str(), ("histogram_title" + std::to_string(strawNum)).c_str(), strawResolution, 0, strawResolution);
  fillHistogram1D(hits, histogram);
  std::pair<std::vector<double>, std::vector<double>> x_y_peaks = findPeaks(histogram);
  std::string filename;

  if (plottingGraphs){
    histogram->Draw();
    canvas->Update();
    filename = file_prefix + "_canvas_" + std::to_string(strawNum) + ".png";
    canvas->SaveAs(filename.c_str());
  }
	
  // refines the peaks found by findPeaks using multiple gaussians fitting
  gaussianFit(x_y_peaks.first, x_y_peaks.second, histogram);

  if (plottingGraphs){
    canvas->Update();
    filename = "gaussian_" + filename;
	  canvas->SaveAs(filename.c_str());
  }
  delete canvas;
  delete histogram;

  // returning just the x values of the peaks - y values no longer relevant
  return x_y_peaks.first;
}


void CalibrationCalculator::gaussianFit(std::vector<double> x, std::vector<double> y, TH1D* histogram){
	
  // setting up function string for length of multi gaussian needed, based on precalculated peaks
  char functionstring[2000];
	sprintf(functionstring, "gaus(0)");
  int npeaks = x.size();
	if (npeaks > 1) {
	  for (int k = 1; k < npeaks; k++) {
		  sprintf(functionstring, "%s + gaus(%d)", functionstring, k * 3);
	  }
	}
	sprintf(functionstring, "%s + pol2(%d)", functionstring, npeaks * 3);
  
  //TODO, set as parameters to the program, and store them in the calibration file
  TF1* function = new TF1("gaus", functionstring, 100, 450);

  // setting up function parameters, a y value, x value, and --- per straw, and then ...
	const int npar = npeaks * 3 + 3;
	double* par = new double[npar];
	for (int ii = 0; ii < npeaks; ii++) {
		par[ii * 3 + 0] = y[ii];
		par[ii * 3 + 1] = x[ii];
		par[ii * 3 + 2] = 1.1;
	}
	par[npeaks * 3 + 0] = 3.1;
	par[npeaks * 3 + 1] = 0.001;
	par[npeaks * 3 + 2] = 0.00001;

  // setting parameters
	function->SetParameters(par);
	function->SetLineColor(kOrange);
	function->SetNpx(100000);

  // fitting function and drawing it on the histogram canvas
	histogram->Fit(function, "E0", "", 100, 450);
  function->Draw("same");
  
  double chi2 = function->GetChisquare();
  double ndf = function->GetNDF();
  double reduced_chi2 = chi2 / ndf;
  std::cout << "reduced chi2 = " + std::to_string(reduced_chi2) << std::endl;
  // // getting the gaussian parameters - we don't use these for anything, could be removed
  // double* results = new double[npar];
  // function->GetParameters(results);


}


void CalibrationCalculator::writeStrawInfoToFile(std::string filename) {
 
  nlohmann::json j = strawInfo;

  // Save the JSON data to a file
  std::ofstream outfile(filename);
  outfile << j.dump(4);  // Use pretty-printing with 4 spaces for indentation
  outfile.close();
}

void CalibrationCalculator::loadStrawInfoFromFile(std::string filename) {
  // Read the JSON data from the file
  std::ifstream infile(filename);
  nlohmann::json j;
  infile >> j;
  infile.close();

  // Parse the JSON data into the strawInfo map
  j.get_to(strawInfo);
}

void CalibrationCalculator::selectionSort(std::vector<double> a, std::vector<double> b, int n) {
	int i, j, min;
	double  tempa, tempb;
	for (i = 0; i < n - 1; i++) {
		min = i;
		for (j = i + 1; j < n; j++)
			if (a[j] < a[min])
				min = j;
		tempa = a[i];
		tempb = b[i];
		a[i] = a[min];
		b[i] = b[min];
		a[min] = tempa;
		b[min] = tempb;
	}
}

void CalibrationCalculator::fillHistogram1D(std::vector<int> hits, TH1D* histogram) {
	std::cout << "Histogram contains " + std::to_string(hits.size()) + " hits \n";
  for (int hit : hits) {
		histogram->Fill(hit);
	}
}


bool CalibrationCalculator::checkRange(std::vector<double> strawCalibrationParams){  
  std::vector<double> range(strawResolution - 40);
  for(int i = 20; i < strawResolution-20; i++){
    range[i-20] = i; 
  }
  std::vector<double> result = applyCalibrationParams(range, strawCalibrationParams);
  for(int i = 0; i < range.size(); i++){
    if (not ((result[i] >= 0) && (result[i] <= strawResolution))){
      return false;
    }
  }
  return true;
}


// mirrors the way calibration parameters are applied in the EFU
std::vector<double> CalibrationCalculator::applyCalibrationParams(std::vector<double> measured, std::vector<double> params){
  std::vector<double> result;
  for(int i = 0; i < measured.size(); i++){
    double adjustedPosition = measured[i];
    for(int j = 0; j < params.size(); j++){
      adjustedPosition -= params[j] * pow(measured[i], j);
    }
    result.push_back(adjustedPosition);
  }
  return result;
}