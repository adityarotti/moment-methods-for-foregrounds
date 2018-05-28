import sys
import numpy as np
from scipy.optimize import curve_fit
import constants as cnst
import analytic_sed as analytic_sed


class moment_fit(object):
    
    def __init__(self,n):
        self.n=n
        self.ana_sed=analytic_sed.analytic_sed()
        self.ana_sed.create_fn_dir(self.n)
        self.ana_sed.fn_dir



    def moment_expansion_function(self,x,*params):
        '''Pass a minimum to 3 parameters: (T, slope,A) to *params else this function will fail.\n Pass 3 parameters to fit a MBB sed \n Pass 6 parameters to fit a MBB with second derivative functions. \n Pass 10 parameters to fit a MBB with thridd derivative functions.'''
 
        A=params[0] ; T=params[1] ; slope=params[2]
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz  # Pivot frequency in GHz
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const # This defines h*nu/k in the exponential
        
        y=A*self.ana_sed.fn_dir[0](x,1./T,slope,nu0,c0)
        for i, arg in enumerate(params[3:]):
            # 3+i is to exclude 1st derivative terms. Remember i starts from 0.
            y=y+A*arg*self.ana_sed.fn_dir[3+i](x,1./T,slope,nu0,c0)
        return y



    def perturbative_moment_expansion_function(self,x,A,T,slope,*params):
        '''This function is used to fit for moments perturbatively around a MBB sed defined by parameters T,Slope,A \n The zeroth parameter fits for the perturbation in the Amplitude about the mean. \n In the perturbative fits we also fit for the first order derivatives, unlike in the case of exact fits.'''
        
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz # Pivot frequency in GHz
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const # This defines h*nu/k in the exponential
        y=self.ana_sed.fn_dir[0](x,1./T,slope,nu0,c0)*A
        for i, par in enumerate(params):
            y=y+A*par*self.ana_sed.fn_dir[i](x,1./T,slope,nu0,c0)
        return y
    
# Fitting carried out in the function defined below. Uses the scipy curve_fit routine.
    
    def fit_monopole_sed(self,nu,Inu,der_order=0,guess=[],lb=[],ub=[],prange=1000.,maxfev=20000,bounds_true=False,flat_sensitivity=True):

        if der_order==1:
            der_order=0
            print "Warning: The coefficients of the first derivative terms are expeced to be zero. Resetting to default der_order=0"
            print "Effectively only fitting for a single modified black body"
            print "To perform moment fits please provide a der_order > 1"
        

        # Parameter guesses and bounds.
        if np.size(guess)==0:
            par_size=max(self.ana_sed.calc_num_vec(der_order),3)
            guess=np.append((1.,10.,0.),np.zeros(par_size-3)) # A,T,slope + moments
            print "You provided no guesses"
            print "Fitting for", par_size, " parameters"
        elif np.size(guess)!=self.ana_sed.calc_num_vec(der_order):
            par_size=max(self.ana_sed.calc_num_vec(der_order),3)
            guess=np.append((1.,10.,0.),np.zeros(par_size-3))
            print "For derivative order =",der_order,"the code requires that you pass a guess array of size ", self.ana_sed.calc_num_vec(der_order), "but you only provided a guess array of size ",np.size(guess)
            print "Resetting the guesses to the default guess values ", guess
        else:
            par_size=np.size(guess)

        # Need to incorporate user provided bounds.
        if bounds_true:
            if np.size(lb)==0 or np.size(ub)==0 or np.size(lb)!=np.size(ub):
                print "Bounds not provided or they lower and upper bound arrays are of unequal size"
                print "Setting the default bounds"
                lb=np.append([1e-6,0.1,-10],np.ones(par_size-3)*(-prange)) # T & A cannot be negative.
                ub=np.append([10,50.,10],np.ones(par_size-3)*(prange))

        # Assumption on sensitivity.
        if flat_sensitivity:
            # This sets a flat error across different frequencies
            # This results in relatively worse fits in regions where the SED is low as compared to where the SED is high.
            err=np.ones(np.size(nu))
        else:
            # This sets the error to be proportional to the SED. This is ideal for working with simulations where one desires the fits to have constant relative error on the SED at all frequencies.
            # The quality of the fit is not affected by the amplitude of the SED.
            err=Inu

        if bounds_true:
            fit_par,fit_par_cov=curve_fit(self.moment_expansion_function ,nu,Inu,p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
        else:
            fit_par,fit_par_cov=curve_fit(self.moment_expansion_function,nu,Inu,p0=guess,sigma=err,maxfev=maxfev)

        return fit_par


    def fit_sed_perturbatively(self,nu,Inu,der_order,A,T,slope,prange=1000.,maxfev=20000,bounds_true=True,flat_sensitivity=True):

        par_size=self.ana_sed.calc_num_vec(der_order)
        print "Fitting for", par_size, " parameters"

        guess=np.zeros(par_size,np.float64) # The length of "guess" is what informs curve fit how many parameters to fit.
        # Need to incorporate user provided bounds.
        if bounds_true:
            lb=np.ones(par_size)*(-prange) # T & A cannot be negative.
            ub=np.ones(par_size)*(prange)

        # Assumption on sensitivity.
        if flat_sensitivity:
            # This sets a flat error across different frequencies
            # This results in relatively worse fits in regions where the SED is low as compared to where the SED is high.
            err=np.ones(np.size(nu))
        else:
            # This sets the error to be proportional to the SED. This is ideal for working with simulations where one desires the fits to have constant relative error on the SED at all frequencies.
            # The quality of the fit is not affected by the amplitude of the SED.
            err=Inu

        def fitting_fn(x,*params):
            return self.perturbative_moment_expansion_function(x,A,T,slope,*params)

        if bounds_true:
            fit_par,fit_par_cov=curve_fit(fitting_fn,nu,Inu,p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
        else:
            fit_par,fit_par_cov=curve_fit(fitting_fn,nu,Inu,p0=guess,sigma=err,maxfev=maxfev)

        return fit_par
