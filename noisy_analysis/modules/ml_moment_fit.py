import sys
import numpy as np
from scipy.optimize import curve_fit
import constants as cnst
import analytic_sed as analytic_sed


class ml_moment_fit(object):
    
    def __init__(self,n):
        self.n=n
        self.ana_sed=analytic_sed.analytic_sed()
        self.ana_sed.create_fn_dir(self.n)

    def get_ml_moment_soln(self,nu,T,slope,n,Inu,err_Inu=[],n_is_der_order=True,inc_1der=True,return_only_rec_inu=False):
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz # Pivot frequency in GHz
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const

        if inc_1der:
            vec_idx_shift=0
        else:
            vec_idx_shift=2

        if n_is_der_order:
            num_vec=self.ana_sed.calc_num_vec(n)
        else:
            num_vec=n
        num_vec=num_vec-vec_idx_shift

        if np.size(err_Inu)==0:
            err_Inu=np.ones(np.size(nu),float)

        vecs=np.zeros((num_vec,np.size(nu)),float)
        cov_vecs=np.zeros((num_vec,num_vec),float)

        vdoti=[]
        for i in range(num_vec):
            k=i
            if i>0:
                k=i+vec_idx_shift
            vecs[i,:]=self.ana_sed.fn_dir[k](nu,1./T,slope,nu0,c0)
            vdoti=np.append(vdoti,np.dot(Inu/err_Inu**2.,vecs[i,:]))
            for j in range(i+1):
                cov_vecs[i,j]=np.dot(vecs[i,:],vecs[j,:]/(err_Inu**2.))
                cov_vecs[j,i]=cov_vecs[i,j]
        
        moments=np.asarray(np.matmul(np.linalg.inv(np.matrix(cov_vecs)),np.transpose(np.matrix(vdoti)))).flatten()
        moments_cov=np.linalg.inv(np.matrix(cov_vecs))
        rec_inu=0
        for i in range(num_vec):
            rec_inu=rec_inu+moments[i]*vecs[i,:]

        if return_only_rec_inu:
            return rec_inu
        else:
            return moments,moments_cov,rec_inu

    
    def get_best_fit_param(self,nu,Inu,err_Inu,n,guess=[10.,0.],inc_1der=False,maxfev=200000):
        fn = lambda nu,T,slope: self.get_ml_moment_soln(nu,T,slope,n=n,Inu=Inu,err_Inu=err_Inu,inc_1der=inc_1der,return_only_rec_inu=True)
        par,par_cov=curve_fit(fn,nu,Inu,sigma=err_Inu,p0=guess,bounds=([0.1,-1.],[20.,3.]),absolute_sigma=True,maxfev=maxfev)
        moments,rec_Inu=self.get_ml_moment_soln(nu,T=par[0],slope=par[1],n=n,Inu=Inu,err_Inu=err_Inu,inc_1der=inc_1der,return_only_rec_inu=False)
        return par,moments,rec_Inu
