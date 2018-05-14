import numpy as np
import constants as cnst
import analytic_sed as analytic_sed
from scipy.optimize import curve_fit
import sys

class gram_schmidt_fitting(object):
    
    def __init__(self,n):
        self.n=n
        self.ana_sed=analytic_sed.analytic_sed()
        self.ana_sed.create_fn_dir(self.n)
        self.ana_sed.fn_dir

    def gen_vectors(self,nu,T,slope):
        self.vectors=list()
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz 
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const
        for i in range(len(self.ana_sed.fn_dir)):
            v=self.ana_sed.fn_dir[i](nu,1./T,slope,nu0,c0)
            self.vectors.append(v/max(abs(v)))
        
    def gram_schmidt(self,vectors):
        basis = []
        un_basis = []
        for v in vectors:
            w = v - np.sum( np.dot(v,b)*b  for b in basis )
            un_basis.append(w)
        
            #Calculate the norm for w-vector.
            norm=np.sqrt(np.dot(w,w)) #; print norm
        
            #Normalize vector.
            if (norm > 1e-12):  
                basis.append(w/norm)         

        return np.array(un_basis),np.array(basis)

    def gram_schmidt_iterative(self,tol=1e-8,iter=10):
        vectors=self.vectors

        normalize=True ; counter=0
        while normalize:
            print counter
            un_basis,self.basis=self.gram_schmidt(vectors) ; normalize=False
        
            if iter>counter:
                un_basis,basis_new=self.gram_schmidt(self.basis)
            
                norm=[]
                for i in range(np.shape(un_basis)[0]):
                    norm=np.append(norm,np.dot(un_basis[i,:],un_basis[i,:]))    
                    #print abs(norm-1.)/tol
            
                if (abs(norm-1.)>tol).any() or counter==iter:
                    vectors=self.basis
                    normalize=True
                    counter=counter+1   

    def get_gram_schmidt_param(self,Inu):
        for i in range()
