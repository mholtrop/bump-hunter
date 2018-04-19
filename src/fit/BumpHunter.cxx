
/** 
 * @file BumpHunter.cxx
 * @brief
 * @author Omar Moreno <omoreno1@ucsc.edu>
 *         Santa Cruz Institute for Particle Physics
 *         University of California, Santa Cruz
 * @date January 14, 2015
 *
 */

#include <BumpHunter.h>

BumpHunter::BumpHunter(BkgModel model, int poly_order, int res_factor) 
    : _model(nullptr),
      signal(nullptr), 
      bkg(nullptr),
      ofs(nullptr),
      _res_factor(res_factor), 
      _poly_order(poly_order) {

    std::cout << "[ BumpHunter ]: Background polynomial: " << _poly_order << std::endl;
    std::cout << "[ BumpHunter ]: Resolution multiplicative factor: " << _res_factor << std::endl;

    // Turn off all messages except errors
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    // Independent variable
    mass_ = new RooRealVar("Invariant Mass", "Invariant Mass (GeV)", 0., 0.15);
    //mass_ = new RooRealVar("Invariant Mass", "Invariant Mass (GeV)", 0., 0.1);

    //   Signal PDF   //
    //----------------//   

    variable_map["A' mass"]  = new RooRealVar("A' Mass",  "A' Mass",  0.03);

    variable_map["A' mass resolution"] 
        = new RooRealVar("A' Mass Resolution", "A' Mass Resolution", this->getMassResolution(0.03));

    signal = new RooGaussian("signal", "signal", *mass_,
                             *variable_map["A' mass"], *variable_map["A' mass resolution"]);

    //   Bkg PDF   //
    //-------------//

    std::string name;
    for (int order = 1; order <= _poly_order; ++order) {
        name = "t" + std::to_string(order);
        if (order == 5) { 
            variable_map[name] = new RooRealVar(name.c_str(), name.c_str(), 0, -0.01, 0.01);
        } else { 
            variable_map[name] = new RooRealVar(name.c_str(), name.c_str(), 0, -1, 1);
        }
        arg_list.add(*variable_map[name]);
    } 
    
    switch(model) { 
        case BkgModel::POLY: {
            std::cout << "[ BumpHunter ]: Modeling the background using a polynomial of order " 
                      << poly_order << std::endl;
            bkg = new RooChebychev("bkg", "bkg", *mass_, arg_list);
        } break;
        case BkgModel::EXP_POLY: {
            std::cout << "[ BumpHunter ]: Modeling the background using an exp(poly of order "
                      << poly_order << ")" << std::endl;
            RooChebychev* exp_poly_bkg 
                = new RooChebychev("exp_poly_bkg", "exp_poly_bkg", *mass_, arg_list);
            bkg = new RooExponential("bkg", "bkg", *exp_poly_bkg, *(new RooRealVar("const", "const", 1))); 
        } break;
        case BkgModel::EXP_POLY_X_POLY: { 
            std::cout << "[ BumpHunter]: Modeling the background using an exp(-cx)*(poly of order "
                      << poly_order << ")" << std::endl;
            variable_map["c"] = new RooRealVar("const", "const", 0, -2, -0.000001);
            RooChebychev* poly_bkg 
                = new RooChebychev("poly_bkg", "poly_bkg", *mass_, arg_list);
            RooExponential* exponential 
                = new RooExponential("exp_bkg", "exp_bkg", *mass_, *variable_map["c"]);
            bkg = new RooProdPdf("bkg", "bkg", *poly_bkg, *exponential); 
        } break;
    }


    //   Composite Models   //
    //----------------------//
    std::cout << "[ BumpHunter ]: Creating composite model." << std::endl;

    //variable_map["signal yield"] = new RooRealVar("signal yield", "signal yield", 0, -1e13, 1e13);
    variable_map["signal yield"] = new RooRealVar("signal yield", "signal yield", 0, -1e13, 1e13);
    //variable_map["signal yield"] = new RooRealVar("signal yield", "signal yield", 0, -100000, 100000);
    //variable_map["bkg yield"] = new RooRealVar("bkg yield", "bkg yield", 30000000, -1e13, 1e13);
    variable_map["bkg yield"] = new RooRealVar("bkg yield", "bkg yield", 30000000, 0, 1e8);

    _model = new RooAddPdf("comp model", "comp model", RooArgList(*signal, *bkg), 
                               RooArgList(*variable_map["signal yield"], *variable_map["bkg yield"]));


    for (auto& element : variable_map) {
        default_values[element.first] = element.second->getVal(); 
        default_errors[element.first] = element.second->getError(); 
    }

    RooRandom::randomGenerator()->SetSeed(time(nullptr)); 
}

BumpHunter::~BumpHunter() {

    delete printer; 

    /*
    for (auto& element : variable_map) { 
       delete element.second; 
    }
    variable_map.clear();*/
    //delete signal;
    //delete bkg;
    //delete _model; 
}

void BumpHunter::initialize(TH1* histogram, double &mass_hypothesis) { 

    mass_hypothesis_ = mass_hypothesis; 

    // Shift the mass hypothesis so it sits in the middle of a bin.  
    std::cout << "[ BumpHunter ]: Shifting mass hypothesis to nearest bin center " 
              << mass_hypothesis << " ---> "; 
    int mbin = histogram->GetXaxis()->FindBin(mass_hypothesis); 
    mass_hypothesis = histogram->GetXaxis()->GetBinCenter(mbin);
    std::cout << mass_hypothesis << " GeV" << std::endl;

    // If the mass hypothesis is below the lower bound, throw an exception.  A 
    // search cannot be performed using an invalid value for the mass hypothesis.
    //if (mass_hypothesis < _lower_bound) {
    //    throw std::runtime_error("Mass hypothesis less than the lower bound!"); 
    //}

    // Set the mean of the Gaussian signal distribution
    //variable_map["A' mass"]->setVal(mass_hypothesis);

    // Correct the mass to take into account the mass scale systematic
    double corr_mass = this->correctMass(mass_hypothesis);

    // Get the mass resolution at the corrected mass 
    mass_resolution_ = this->getMassResolution(corr_mass);
    std::cout << "[ BumpHunter ]: Mass resolution: " << mass_resolution_ << " GeV" << std::endl;

    // Set the width of the Gaussian signal distribution
    //variable_map["A' mass resolution"]->setVal(mass_resolution); 

    // Calculate the fit window size
    window_size_ = mass_resolution_*_res_factor;
    this->printDebug("Window size: " + std::to_string(window_size_));

    // Find the starting position of the window. This is set to the low edge of 
    // the bin closest to the calculated value. If the start position falls 
    // below the lower bound of the histogram, set it to the lower bound.
    window_start_ = mass_hypothesis - window_size_/2;
    int window_start_bin = histogram->GetXaxis()->FindBin(window_start_);  
    window_start_ = histogram->GetXaxis()->GetBinLowEdge(window_start_bin);
    if (window_start_ < _lower_bound) { 
        std::cout << "[ BumpHunter ]: Starting edge of window (" << window_start_ 
                  << " MeV) is below lower bound." << std::endl;
        window_start_bin = histogram->GetXaxis()->FindBin(_lower_bound);
        window_start_ = histogram->GetXaxis()->GetBinLowEdge(window_start_bin);
    }
    std::cout << "[ BumpHunter ]: Setting starting edge of fit window to " 
              << window_start_ << " MeV." << std::endl;
    
    // Find the end position of the window.  This is set to the upper edge of 
    // the bin closest to the calculated value. If the window edge falls above
    // the upper bound of the histogram, set it to the upper bound.
    // Furthermore, check that the bin serving as the upper boundary contains
    // events. If the upper bound is shifted, reset the lower window bound.
    //window_end_ = mass_hypothesis + window_size_/2;
    window_end_ = window_start_ + window_size_;
    int window_end_bin = histogram->GetXaxis()->FindBin(window_end_);
    window_end_ = histogram->GetXaxis()->GetBinLowEdge(window_end_bin);
    /*if (window_end_ > _upper_bound) { 
        std::cout << "[ BumpHunter ]: Upper edge of window (" << window_end_ 
                  << " MeV) is above upper bound." << std::endl;
        window_end_bin = histogram->GetXaxis()->FindBin(_upper_bound);
        
        int last_bin_above = histogram->FindLastBinAbove(); 
        if (window_end_bin > last_bin_above) window_end_bin = last_bin_above; 
        
        window_end_ = histogram->GetXaxis()->GetBinLowEdge(window_end_bin);
        window_start_bin = histogram->GetXaxis()->FindBin(window_end_ - window_size_);  
        window_start_ = histogram->GetXaxis()->GetBinLowEdge(window_start_bin);
    }*/
    std::cout << "[ BumpHunter ]: Setting upper edge of fit window to " 
              << window_end_ << " MeV." << std::endl;

    // Set the range that will be used in the fit
    //range_name_ = "mass_" + std::to_string(mass_hypothesis) + "gev"; 
    //mass_->setRange(range_name_.c_str(), window_start_, window_end_); 
    //mass_->setRange(window_start_, window_end_); 

    // Estimate the background normalization within the window by integrating
    // the histogram within the window range.  This should be close to the 
    // background yield in the case where there is no signal present.
    integral_ = histogram->Integral(window_start_bin, window_end_bin);
    this->printDebug("Window integral: " + std::to_string(integral_)); 

    //variable_map["bkg yield"]->setVal(integral_);
    //variable_map["bkg yield"]->setError(round(sqrt(integral_))); 
    //default_values["bkg yield"] = integral_;
    //default_errors["bkg yield"] = sqrt(integral_);  

    //variable_map["signal yield"]->setError(sqrt(integral_));
    //default_errors["signal yield"] = sqrt(integral_); 

    // Calculate the size of the background window as 2.56*(mass_resolution)
    //_bkg_window_size = std::trunc(mass_resolution*2.56*10000)/10000 + 0.00005;

    // Find the starting position of the bkg window
    //double bkg_window_start = mass_hypothesis - _bkg_window_size/2;

    //int bkg_window_start_bin = histogram->GetXaxis()->FindBin(bkg_window_start);  
    //int bkg_window_end_bin = histogram->GetXaxis()->FindBin(bkg_window_start + _bkg_window_size);  

    //_bkg_window_integral = histogram->Integral(bkg_window_start_bin, bkg_window_end_bin); 
    
}

HpsFitResult* BumpHunter::performSearch(TH1* histogram, double mass_hypothesis, bool skip_bkg_fit = false) { 
   
    // Calculate all of the fit parameters e.g. window size, mass hypothesis
    this->initialize(histogram, mass_hypothesis);

    HpsFitResult* fit_result = new HpsFitResult(); 

    // Create a histogram object compatible with RooFit.
    //data_ = new RooDataHist("data", "data", RooArgList(*mass_), histogram);

    // Get the number of bins in the resulting histogram.  This will be used
    // when generating toys.
    //bins_ = data_->numEntries();
    //this->printDebug("Total number of bins: " + std::to_string(bins_)); 
   
    
    if (!skip_bkg_fit)  {
         
        // Start by performing a background only fit.  The results from this fit 
        // are used to get an intial estimate of the background parameters.  
        std::cout << "*************************************************" << std::endl;
        std::cout << "*************************************************" << std::endl;
        std::cout << "[ BumpHunter ]: Performing a background only fit." << std::endl;
        std::cout << "*************************************************" << std::endl;
        std::cout << "*************************************************" << std::endl;
   
        //
        BkgFunction bkg_func(mass_hypothesis, window_end_ - window_start_); 
        TF1* bkg = new TF1("bkg", bkg_func, -1, 1, 6);

        //
        bkg->SetParameters(4,0,0,0,0,0);
        bkg->SetParNames("pol0","pol1","pol2","pol3","pol4","pol5");

        TFitResultPtr result = histogram->Fit("bkg", "LES+", "", window_start_, window_end_); 
        fit_result->setBkgFitResult(result); 
        result->Print();

    }
         
    std::cout << "***************************************************" << std::endl;
    std::cout << "***************************************************" << std::endl;
    std::cout << "[ BumpHunter ]: Performing a signal+background fit." << std::endl;
    std::cout << "***************************************************" << std::endl;
    std::cout << "***************************************************" << std::endl;
        
    FullFunction full_func(mass_hypothesis, window_end_ - window_start_); 
    TF1* full = new TF1("full", full_func, -1, 1, 9);
    
    full->SetParameters(4,0,0,0,0,0,0,0,0);
    full->SetParNames("pol0","pol1","pol2","pol3","pol4","pol5","signal norm","mean","sigma");
    full->FixParameter(7,0.0); 
    full->FixParameter(8, mass_resolution_); 
       
    TFitResultPtr full_result = histogram->Fit("full", "LES+", "", window_start_, window_end_);
    fit_result->setCompFitResult(full_result);  
    full_result->Print();

    //if (!batch) { 
        printer->print(histogram, window_start_, window_end_, "test_print.pdf");
    //}

    this->calculatePValue(fit_result);
    //this->getUpperLimit(histogram, fit_result);
    
    //if (!_batch) { 
    //    std::string output_path = "fit_result_" + std::string(histogram->GetName()) 
      //                            + "_" + std::to_string(mass_hypothesis) + "gev_bkg_only.png";
     //   printer->print(mass_, data_, _model, range_name_, output_path); 
        /*if (_write_results) { 
     
            // Create the output file name string
            char buffer[100];
            std::string output_file = "fit_result_" + std::string(histogram->GetName()) 
                                      + "_" + std::to_string(mass_hypothesis) + "gev" + (bkg_only ? "_bkg" : "_full") + ".txt";
            std::cout << "[ BumpHunter ]: Writing results to " << output_file << std::endl;
            sprintf(buffer, output_file.c_str()); 

            // Create a file stream  
            ofs = new std::ofstream(buffer, std::ofstream::out); 
            result->getRooFitResult()->printMultiline(*ofs, 0, kTRUE, "");

            ofs->close();
        }*/
    //}
    
    //}

    // Persist the mass hypothesis used for this fit
    //result->setMass(mass_hypothesis); 

    // Set the window size 
    //result->setWindowSize(window_size_);

    // Set the total number of events within the window
    //result->setIntegral(default_values["bkg yield"]);

    // Calculate the size of the background window as 2.56*(mass_resolution)
    //result->setBkgWindowSize(_bkg_window_size); 

    //result->setBkgTotal(_bkg_window_integral); 

    //result->setNBins(mass_->getBinning().numBins()); 


    //result->setCorrectedMass(this->correctMass(mass_hypothesis)); 

    return fit_result; 
}

/*HpsFitResult* BumpHunter::fit(RooDataHist* data, std::string range_name = "") { 
   
    RooFitResult* result = _model->fitTo(*data, RooFit::Range(range_name.c_str()), RooFit::Extended(kTRUE), 
                RooFit::SumCoefRange(range_name.c_str()), RooFit::Save(kTRUE), RooFit::PrintLevel(-1000)); 
  
    // Return the saves result
    return new HpsFitResult(result); 
}*/

void BumpHunter::calculatePValue(HpsFitResult* result) {

    std::cout << "[ BumpHunter ]: Calculating p-value: " << std::endl;

    //
    double signal_yield = result->getCompFitResult()->Parameter(6); 
    this->printDebug("Signal yield: " + std::to_string(signal_yield)); 

    // In searching for a resonance, a signal is expected to lead to an 
    // excess of events.  In this case, a signal yield of < 0 is  
    // meaningless so we set the p-value = 1.  This follows the formulation
    // put forth by Cowen et al. in https://arxiv.org/pdf/1007.1727.pdf. 
    if (signal_yield <= 0) { 
        result->setPValue(1);
        this->printDebug("Signal yield is negative ... setting p-value = 1"); 
        return; 
    }

    // Get the NLL obtained by minimizing the composite model with the signal
    // yield floating.
    double mle_nll = result->getCompFitResult()->MinFcnValue(); 
    this->printDebug("NLL when mu = " + std::to_string(signal_yield) + ": " + std::to_string(mle_nll));

    // Get the NLL obtained from the Bkg only fit.
    double cond_nll = result->getBkgFitResult()->MinFcnValue(); 
    printDebug("NLL when mu = 0: " + std::to_string(cond_nll));

    // 1) Calculate the likelihood ratio which is chi2 distributed. 
    // 2) From the chi2, calculate the p-value.
    double q0 = 0; 
    double p_value = 0; 
    this->getChi2Prob(cond_nll, mle_nll, q0, p_value);  

    std::cout << "[ BumpHunter ]: p-value: " << p_value << std::endl;
    
    // Update the result
    result->setPValue(p_value);
    result->setQ0(q0);  
    
}

void BumpHunter::printDebug(std::string message) { 
    if (debug) std::cout << "[ BumpHunter ]: " << message << std::endl;
}

void BumpHunter::getUpperLimit(TH1* histogram, HpsFitResult* result) {

    FullFunction comp_func(mass_hypothesis_, window_end_ - window_start_); 
    TF1* comp = new TF1("comp_ul", comp_func, -1, 1, 9);   
    comp->Print(); 
    
    comp->SetParameters(4,0,0,0,0,0,0,0,0);
    comp->SetParNames("pol0","pol1","pol2","pol3","pol4","pol5","signal norm","mean","sigma");
    comp->FixParameter(7,0.0); 
    comp->FixParameter(8, mass_resolution_); 
    comp->Print(); 
    std::cout << "Mass resolution: " << mass_resolution_ << std::endl; 
    std::cout << "[ BumpHunter ]: Calculating upper limit." << std::endl;

    //  Get the signal yield obtained from the signal+bkg fit
    double signal_yield = result->getCompFitResult()->Parameter(6); 
    this->printDebug("Signal yield: " + std::to_string(signal_yield)); 

    // Get the minimum NLL value that will be used for testing upper limits.
    // If the signal yield (mu estimator) at the min NLL is < 0, use the NLL
    // obtained when mu = 0.
    double mle_nll = result->getCompFitResult()->MinFcnValue(); 

    if (signal_yield < 0) {
        this->printDebug("Signal yield @ min NLL is < 0. Using NLL when signal yield = 0");

        // Get the NLL obtained assuming the background only hypothesis
        mle_nll = result->getBkgFitResult()->MinFcnValue(); 
        
        signal_yield = 0;
    } 
    this->printDebug("MLE NLL: " + std::to_string(mle_nll));   

    signal_yield = floor(signal_yield) + 1; 

    double p_value = 1;
    double q0 = 0; 
    while(true) {
       
        this->printDebug("Setting signal yield to: " + std::to_string(signal_yield));
        std::cout << "[ BumpHunter ]: Current p-value: " << p_value << std::endl;
        comp->FixParameter(6, signal_yield); 
        
        TFitResultPtr full_result 
            = histogram->Fit("comp_ul", "LES+", "", window_start_, window_end_);
        double cond_nll = full_result->MinFcnValue(); 

        // 1) Calculate the likelihood ratio which is chi2 distributed. 
        // 2) From the chi2, calculate the p-value.
        this->getChi2Prob(cond_nll, mle_nll, q0, p_value);  

        this->printDebug("p-value after fit : " + std::to_string(p_value)); 
        
        if ((p_value <= 0.0455 && p_value > 0.044)) { 
            
            std::cout << "[ BumpHunter ]: Upper limit: " << signal_yield << std::endl;
            std::cout << "[ BumpHunter ]: p-value: " << p_value << std::endl;

            result->setUpperLimit(signal_yield);
            result->setUpperLimitPValue(p_value); 
     
            break; 
        } else if (p_value <= 0.044) {
            signal_yield -= 1;
        } else if (p_value <= 0.059) signal_yield += 1;
        else if (p_value <= 0.10) signal_yield += 20;
        else if (p_value <= 0.2) signal_yield += 40; 
        else signal_yield += 100;  
    }
    /*
    this->printDebug("p-value from result: " + std::to_string(p_value));
    double q0 = 0;
    int fit_counter = 1;
    bool fell_below_threshold = false; 


        this->resetParameters();
        variable_map["signal yield"]->setConstant(kFALSE); 
        variable_map["signal yield"]->setVal(signal_yield);
        variable_map["signal yield"]->setConstant(kTRUE);
        std::cout << "[ BumpHunter ]: Setting signal to " << variable_map["signal yield"]->getValV() << std::endl; 

        HpsFitResult* current_result = this->fit(data, range_name); 
        
        double cond_nll = current_result->getRooFitResult()->minNll(); 

        result->addSignalYield(signal_yield); 
        result->addLikelihood(cond_nll); 

        this->getChi2Prob(cond_nll, mle_nll, q0, p_value);  

    

        ++fit_counter; 

        
        delete current_result; 
    }*/
}

std::vector<TH1*> BumpHunter::generateToys(TH1* histogram, double n_toys) { 

    /*
    variable_map["signal yield"]->setConstant(kFALSE);
    this->resetParameters();

    RooRealVar* toy_mass = new RooRealVar("Invariant Mass", "Invariant Mass (GeV)", 0., 0.15);
    RooDataHist* data = new RooDataHist("toy_data", "toy_data", RooArgList(*toy_mass), histogram);
    
    // Set the range that will be used in the fit
    std::string range_name = "toy_mass_" + std::to_string(mass_hypothesis_) + "gev"; 
    toy_mass->setRange(range_name.c_str(), window_start_, window_end_); 
    toy_mass->setRange(window_start_, window_end_);
    toy_mass->setBins(bins_);  

    // Fix the signal yield at 0.
    variable_map["signal yield"]->setConstant(kTRUE);
    bkg_only_result_ = this->fit(data, range_name);

    std::vector<HpsFitResult*> toy_results;
    for (int toy_n = 0; toy_n < n_toys; ++toy_n) { 
            //std::cout << "Toy: " << toy_n << std::endl;
            RooDataHist* hist = _model->generateBinned(RooArgSet(*toy_mass), 
                                 integral_, RooFit::Range(range_name.c_str()));  
            hists.push_back(hist->createHistogram(("toy_" + std::to_string(toy_n)).c_str(), *toy_mass)); 
            //this->resetParameters(); 
            //TH1* rhist = hist->createHistogram(("toy_" + std::to_string(toy_n)).c_str(), *mass_); 
            //delete hist;
            //RooDataHist* data = new RooDataHist("toy_data", "toy_data", RooArgList(*mass_), rhist);
            //toy_results.push_back(this->fit(data, range_name_));
            //delete data; 

    }
    variable_map["signal yield"]->setConstant(kFALSE);
   
    delete toy_mass;
    delete data; 

    //return toy_results;
    */
    std::vector<TH1*> hists;
    return hists;
}


void BumpHunter::getChi2Prob(double cond_nll, double mle_nll, double &q0, double &p_value) {
   

    this->printDebug("Cond NLL: " + std::to_string(cond_nll)); 
    this->printDebug("Uncod NLL: " + std::to_string(mle_nll));  
    double diff = cond_nll - mle_nll;
    this->printDebug("Delta NLL: " + std::to_string(diff));
    
    q0 = 2*diff;
    this->printDebug("q0: " + std::to_string(q0));
    
    p_value = ROOT::Math::chisquared_cdf_c(q0, 1);
    this->printDebug("p-value before dividing: " + std::to_string(p_value));  
    p_value *= 0.5;
    this->printDebug("p-value: " + std::to_string(p_value)); 
}

void BumpHunter::setBounds(double lower_bound, double upper_bound) {
    _lower_bound = lower_bound; 
    _upper_bound = upper_bound;
    printf("Fit bounds set to [ %f , %f ]\n", _lower_bound, _upper_bound);   
}

double BumpHunter::correctMass(double mass) { 
    double offset = -1.19892320e4*pow(mass, 3) + 1.50196798e3*pow(mass,2) 
                    - 8.38873712e1*mass + 6.23215746; 
    offset /= 100; 
    this->printDebug("Offset: " + std::to_string(offset)); 
    double cmass = mass - mass*offset; 
    this->printDebug("Corrected Mass: " + std::to_string(cmass)); 
    return cmass;
}

/**
 * BkgFunction
 */

BkgFunction::BkgFunction(double mass_hypothesis, double window_size)
    : mass_hypothesis_(mass_hypothesis), 
      window_size_(window_size) { 
}

double BkgFunction::operator() (double* x, double* par) { 
    
    double xp = (x[0] - mass_hypothesis_)/(window_size_*2.0); 
  
    // Chebyshevs between given limits
    double t0 = par[0];
    double t1 = par[1]*xp;
    double t2 = par[2]*(2*xp*xp - 1);
    double t3 = par[3]*(4*xp*xp*xp - 3*xp);
    double t4 = par[4]*(8*xp*xp*xp*xp - 8*xp*xp + 1);
    double t5 = par[5]*(16*xp*xp*xp*xp*xp - 20*xp*xp*xp + 5*xp);
  
    double pol = t0+t1+t2+t3+t4+t5;

    return TMath::Power(10,pol);
}

FullFunction::FullFunction(double mass_hypothesis, double window_size)
    : mass_hypothesis_(mass_hypothesis), 
      window_size_(window_size) { 
}


double FullFunction::operator() (double* x, double* par) { 
    
    double xp = (x[0] - mass_hypothesis_)/(window_size_*2.0); 
  
    // Chebyshevs between given limits
    double t0 = par[0];
    double t1 = par[1]*xp;
    double t2 = par[2]*(2*xp*xp - 1);
    double t3 = par[3]*(4*xp*xp*xp - 3*xp);
    double t4 = par[4]*(8*xp*xp*xp*xp - 8*xp*xp + 1);
    double t5 = par[5]*(16*xp*xp*xp*xp*xp - 20*xp*xp*xp + 5*xp);
  
    double pol = t0+t1+t2+t3+t4+t5;

    double gauss = (1.0)/(sqrt(2.0*TMath::Pi()*pow(par[8],2))) *
        TMath::Exp( - pow((xp-par[7]),2)/(2.0*pow(par[8],2)) );
  
    return TMath::Power(10,pol)+0.0001*par[6]*gauss;
}
