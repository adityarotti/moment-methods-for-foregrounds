import numpy as np
import constants as cnst
import analytic_sed as analytic_sed
from scipy.optimize import curve_fit
import sys

class moment_fit(object):
    
    def __init__(self,n):
        self.n=n
        self.ana_sed=analytic_sed.analytic_sed()
        self.ana_sed.create_fn_dir(self.n)
        self.ana_sed.fn_dir

    def moment_exact(self,x,*params):
        '''Pass a minimum to 3 parameters: (T, slope,A) to *argv else this function will fail.\n Pass 3 parameters to fit a MBB sed \n Pass 6 parameters to fit a MBB with second derivative functions. \n Pass 10 parameters to fit a MBB with thridd derivative functions.'''
 
        A=params[0] ; T=params[1] ; slope=params[2]
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz 
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const
        
        y=A*self.ana_sed.fn_dir[0](x,1./T,slope,nu0,c0)
        for i, arg in enumerate(params[3:]):
            # 3+i is to exclude 1st derivative terms.
            y=y+A*arg*self.ana_sed.fn_dir[3+i](x,1./T,slope,nu0,c0)
        return y

    def moment_perturbative(self,x,A,T,slope,*params):
        '''This function is used to fit for moments perturbatively around a MBB sed defined by parameters T,Slope,A'''
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz 
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const
        y=self.ana_sed.fn_dir[0](x,1./T,slope,nu0,c0)*A
        for i, par in enumerate(params):
            y=y+A*par*self.ana_sed.fn_dir[i](x,1./T,slope,nu0,c0)
        return y
    
    def fit_monopole(self,nu,Inu,der_order=3,guess=[],lb=[],ub=[],prange=1000.,maxfev=20000,bounds_true=False,flat_sensitivity=True):

        if der_order==1:
            print "Warning: The derivative order has to less than 1 or greater than one"

        # Parameter guesses and bounds.
        if np.size(guess)==0:
            par_size=max(self.ana_sed.calc_num_vec(der_order),3)
            guess=np.append((1.,10.,0.),np.zeros(par_size-3)) # A,T,slope
            print "You provided no guesses"
            print "Fitting for", par_size, " parameters"
        else:
            par_size=np.size(guess)

        # Need to incorporate user provided bounds.
        if bounds_true:
            lb=np.append([1e-6,0.1,-10],np.ones(par_size-3)*(-prange)) # T & A cannot be negative.
            ub=np.append([10,50.,10],np.ones(par_size-3)*(prange))

        # Assumption on sensitivity.
        if flat_sensitivity:
            err=np.ones(np.size(nu))
        else:
            err=Inu

        if bounds_true:
            fit_par,fit_par_cov=curve_fit(self.moment_exact,nu,Inu,p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
        else:
            fit_par,fit_par_cov=curve_fit(self.moment_exact,nu,Inu,p0=guess,sigma=err,maxfev=maxfev)

        return fit_par


    def fit_around_monopole(self,nu,Inu,der_order,A,T,slope,prange=1000.,maxfev=20000,bounds_true=True,flat_sensitivity=True):

        par_size=self.ana_sed.calc_num_vec(der_order)
        print "Fitting for", par_size, " parameters"

        guess=np.zeros(par_size,np.float64) # The length of "guess" is what informs curve fit how many parameters to fit.
        # Need to incorporate user provided bounds.
        if bounds_true:
            lb=np.ones(par_size)*(-prange) # T & A cannot be negative.
            ub=np.ones(par_size)*(prange)

        # Assumption on sensitivity.
        if flat_sensitivity:
            err=np.ones(np.size(nu))
        else:
            err=Inu

        def fitting_fn(x,*params):
            return self.moment_perturbative(x,A,T,slope,*params)

        if bounds_true:
            fit_par,fit_par_cov=curve_fit(fitting_fn,nu,Inu,p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
        else:
            fit_par,fit_par_cov=curve_fit(fitting_fn,nu,Inu,p0=guess,sigma=err,maxfev=maxfev)

        return fit_par
