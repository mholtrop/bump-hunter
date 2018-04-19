/**
 * @file FitPrinter.cxx
 * @brief Class used to print a RooFit result to a file.
 * @author Omar Moreno
 *         SLAC National Accelerator Laboratory 
 */

#include "FitPrinter.h"

FitPrinter::FitPrinter() {
}

FitPrinter::~FitPrinter() {
    delete _canvas;
}

void FitPrinter::print(TH1* histogram, double window_start, double window_end, std::string output_path) {

    TColor* fte_blue   = new TColor(5001, 0, 143, 213);
    TColor* fte_orange = new TColor(5002, 252, 79, 48);
    TColor* fte_green  = new TColor(5003, 109, 144, 79); 

    histogram->SetLineColor(kAzure+2); 
    histogram->SetMarkerColor(kAzure+2); 
    histogram->SetTitle(""); 
    histogram->GetXaxis()->CenterTitle(); 
    histogram->Draw("hist e0");

    TF1* bkg = histogram->GetFunction("bkg"); 
    bkg->SetLineColor(kGreen-1);
    bkg->Draw("same");

    TF1* comp = histogram->GetFunction("full"); 
    comp->SetLineColor(kRed+3);  
    comp->Draw("same"); 

    histogram->GetXaxis()->SetRangeUser(window_start, window_end);  

    std::cout << "[ FitPrinter ]: Saving file to: " << output_path << std::endl;
    _canvas->SetLogy(); 
    _canvas->SaveAs(output_path.c_str());

    /*

    RooPlot* plot = var->frame(); 
    data->plotOn(plot, 
                 RooFit::MarkerSize(.4),
                 RooFit::MarkerColor(kAzure+2),
                 RooFit::LineColor(kAzure+2), 
                 RooFit::Name("data")
                 );
    model->plotOn(plot,
                  RooFit::Range(range.c_str()), 
                  RooFit::NormRange(range.c_str()), 
                  RooFit::LineColor(kRed + 2), 
                  RooFit::Name("model")
                  ); 
    model->plotOn(plot, 
                  RooFit::Range(range.c_str()), 
                  RooFit::NormRange(range.c_str()),  
                  RooFit::Components("bkg"), 
                  RooFit::LineStyle(kDashed),
                  RooFit::LineColor(kGreen + 2));
    model->plotOn(plot, 
                  RooFit::Range(range.c_str()), 
                  RooFit::NormRange(range.c_str()),  
                  RooFit::Components("signal"), 
                  RooFit::LineColor(kOrange+10));
    plot->GetYaxis()->SetTitleOffset(1.6);
    plot->SetTitle(range.c_str());

    RooHist* hist = (RooHist*) plot->findObject("data"); 
    RooCurve* curve = (RooCurve*) plot->findObject("model"); 
    double sum = 0; 
    for (int i = 0; i < hist->GetN(); ++i) { 
        double x, y; 
        hist->GetPoint(i, x, y); 
        double y_func = curve->interpolate(x); 
        sum += (y - y_func)*(y-y_func)/y_func;
    }
    std::cout << "Chi2: " << sum << std::endl;
    std::cout << "Chi2/dof: " << sum/(hist->GetN()-6) << std::endl;

    RooPlot* tmp_plot = var->frame(); 
    data->plotOn(tmp_plot);//, RooFit::Binning(int(n_bins/4)));
    model->plotOn(tmp_plot, 
                  RooFit::Range(range.c_str()), 
                  RooFit::NormRange(range.c_str()),  
                  RooFit::Components("bkg"), 
                  RooFit::LineStyle(kDashed),
                  RooFit::LineColor(kGreen));

    RooHist* residuals = tmp_plot->residHist();
    residuals->SetMarkerSize(.4);

    RooPlot* res_plot = var->frame();
    res_plot->addPlotable(residuals, "P");
    */
    /*model->plotOn(res_plot, 
                  RooFit::Range(range.c_str()), 
                  RooFit::NormRange(range.c_str()),  
                  RooFit::Components("signal"), 
                  RooFit::LineStyle(kDashed),
                  RooFit::LineColor(kRed));*/
    /*
    res_plot->SetMinimum(-1000);
    res_plot->SetMaximum(1000);
    res_plot->SetTitle("");
    res_plot->GetYaxis()->SetTitle("Residual");
    res_plot->GetYaxis()->SetTitleSize(0.1);
    res_plot->GetYaxis()->SetLabelSize(0.1);
    res_plot->GetYaxis()->SetTitleOffset(.5);
    
    _main_pad->cd();
    plot->Draw();

    _res_pad->cd(); 
    res_plot->Draw(); 

    */
}
