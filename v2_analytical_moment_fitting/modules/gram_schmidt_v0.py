import sys
import numpy as np
import constants as cnst
import analytic_sed as analytic_sed

class gram_schmidt_fitting(object):
    
    def __init__(self,n):
        self.n=n
        self.ana_sed=analytic_sed.analytic_sed()
        self.ana_sed.create_fn_dir(self.n)
        self.ana_sed.fn_dir
        self.num_basis=0

    def gen_vectors(self,nu,T,slope):
        nu0=cnst.boltzman_const*T/cnst.planck_const/cnst.ghz2hz # Pivot frequency in GHz
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const # This defines h*nu/k in the exponential
        vectors=list()
        for i in range(len(self.ana_sed.fn_dir)):
            v=self.ana_sed.fn_dir[i](nu,1./T,slope,nu0,c0)
            vectors.append(v/max(abs(v)))
        return vectors

    def gram_schmidt(self,vectors,min_vec_norm=1e-16):
        basis = [] # Directory of normalized basis functions
        un_basis = []  # Directory of un-normalized basis functions
        for v in vectors:
            w = v - np.sum( np.dot(v,b)*b  for b in basis )
            un_basis.append(w)
            
            #Calculate the norm for the un-normalized basis vector w.
            norm=np.dot(w,w)
                
            #Normalize the basis vector.
            if (np.sqrt(norm) > min_vec_norm):
                basis.append(w/np.sqrt(norm))
            # Note: If one has sampled at only N frequencies, one cannot have more than N basis.
        return np.array(un_basis),np.array(basis)

    def gram_schmidt_iterative(self,nu,T,slope,tol=1e-8,max_iter=10):
        vectors=self.gen_vectors(nu,T,slope)

        normalize=True ; counter=0
        while normalize:
            #print counter
            un_basis,self.basis=self.gram_schmidt(vectors) ; normalize=False
        
            if max_iter>counter:
                un_basis,basis_new=self.gram_schmidt(self.basis)
            
                norm=[]
                for i in range(np.shape(un_basis)[0]):
                    norm=np.append(norm,np.dot(un_basis[i,:],un_basis[i,:]))    
                    #print abs(norm-1.)/tol
            
                if (abs(norm-1.)>tol).any() or counter==max_iter:
                    vectors=self.basis
                    normalize=True
                    counter=counter+1   
        self.num_basis=np.shape(self.basis)[0]


    def get_basis_coeffs(self,Inu,n,n_is_der_order=True):

        if n_is_der_order and  n>self.n:
            n=self.n
            print "The derivatives are only set up to order ", self.n
            print "Resetting n =", self.n
        
        par=[]
        if n_is_der_order:
            num_par=self.ana_sed.calc_num_vec(n)
        else:
            num_par=n

        if self.num_basis<num_par:
            num_par=self.num_basis
            print "Have only", num_par, "basis function, so will fit only as many parameters."

        for i in range(num_par):
        #for i in range(der_order): # These are not derivative order, but number of basis
            par.append(np.dot(self.basis[i],Inu))
        return par

    def reconstruct_sed(self,*coeffs):
        sed=0.
        for i,par in enumerate(coeffs):
            sed=sed+par*self.basis[i]
        return sed
        # It has to be the gram-schmidt basis function here and not the raw vectors.
