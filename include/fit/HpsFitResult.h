/**
 *
 *
 *
 */

#ifndef __HPS_FIT_RESULT_H__
#define __HPS_FIT_RESULT_H__

#include <TFitResult.h>
#include <TFitResultPtr.h>

//------------//
//   RooFit   //
//------------//
#include <RooFitResult.h>
#include <RooRealVar.h>

class HpsFitResult { 

    public: 

        /** Default constructor */
        HpsFitResult(); 

        ~HpsFitResult(); 
      
        TFitResultPtr getBkgFitResult() { return bkg_result_; }; 
        
        TFitResultPtr getCompFitResult() { return comp_result_; }; 

        void addLikelihood(double likelihood) { _likelihoods.push_back(likelihood); }

        void addSignalYield(double signal_yield) { _signal_yields.push_back(signal_yield); }

        double getCorrectedMass() const { return _cmass; }; 

        /** @return _integral The integral within the fit window. */
        double getIntegral() { return _integral; };

        /** @return The mass hypothesis used for this fit. */
        double getMass() const { return _mass; };

        /** */
        double getQ0() { return q0; };

        /** */
        double getPValue() { return p_value; };

        /** @return The signal yield obstained from the sig+bkg fit. */
        float getSignalYield() { return comp_result_->Parameter(poly_order_ + 1);  };
         
        /** @return The error on the signal yield. */
        float getSignalYieldErr() { return comp_result_->ParError(poly_order_ + 1); }; 

        /** */
        double getUpperLimit() { return upper_limit; };

        double getUpperLimitPValue() { return _upper_limit_p_value; }; 

        /** @return The size of the fit window used. */
        double getWindowSize() { return this->window_size; }; 

        std::vector<double> getLikelihoods() { return _likelihoods; }
        
        std::vector<double> getSignalYields() { return _signal_yields; }

        //-------------//
        //   Setters   //
        //-------------//

        void setIntegral(double integral) { _integral = integral; };

        /** */
        double setQ0(double q0) { this->q0 = q0; };

        /** */
        void setPValue(double p_value) { this->p_value = p_value; };  
        
        /** 
         * Set the result from the background only fit.
         *
         * @param result Result from the background only fit.
         */
        void setBkgFitResult(TFitResultPtr bkg_result) { bkg_result_ = bkg_result; };


        /** 
         * Set the result from the signal+background only fit.
         *
         * @param result Result from the signal+background only fit.
         */
        void setCompFitResult(TFitResultPtr comp_result) { comp_result_ = comp_result; };

        /** 
         * Set mass hypothesis used for this fit. 
         * 
         * @param mass The mass hypothesis. 
         */
        void setMass(double mass) { _mass = mass; };

        /**
         *
         *
         */
        void setCorrectedMass(double cmass) { _cmass = cmass; }; 

        /** 
         * Set the order polynomial used by the fitter.
         *
         * @param poly_order Polynomial order used by the fitter.
         */
        void setPolyOrder(int poly_order) { poly_order_ = poly_order; };

        /**
         * Set the 2 sigma upper limit.
         *
         * @param upper_limit The 2 sigma upper limit.
         */
        void setUpperLimit(double upper_limit) { this->upper_limit = upper_limit; };

        /** Set the p-value at the end of the upper limit calculation. */
        void setUpperLimitPValue(double upper_limit_p_value) { _upper_limit_p_value = upper_limit_p_value; }

        /**
         * Set the size of the fit window used to get this results.
         *
         * @oaram window_size The size of the fit window.
         */
        void setWindowSize(double window_size) { this->window_size = window_size; }; 

    private: 

        /** Result from fit using a background model only. */
        TFitResultPtr bkg_result_{nullptr}; 

        /** Result from fit using a signal+background model. */
        TFitResultPtr comp_result_{nullptr}; 

        std::vector<double> _likelihoods; 

        std::vector<double> _signal_yields; 

        /** Total number of events within the fit window. */
        double _integral{0};

        /** Mass hypothesis. */
        double _mass{0}; 

        double _cmass{0}; 

        /** q0 value */
        double q0;
        
        /** p-value. */
        double p_value;

        /** Order polynomial used by the fitter. */
        double poly_order_{0}; 

        /** 2 sigma upper limit on the signal. */
        double upper_limit; 

        /** p-value at the upper limit. */
        double _upper_limit_p_value{-9999};

        /*** Size of the fit window. */
        double window_size;

}; // HpsFitResult

#endif // __HPS_FIT_RESULT_H__
