import numpy as np
import constants as cnst
import analytic_sed as analytic_sed
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import sys


#Notes: 
# In this module, we construct Gram-Schmidt functions at some very high sampling between nu_min and nu_max.
# We interpolate these function frequencies at which SED measurements are available.
# We find the coeffiecients of expansion in the Gram-Schmidt basis by dotting these sparse basis vectors to the SED.
# Note that these sparse function are no longer orthonormal and hence we need to invert the covariance of this sparse basis
# arrive at the effective coefficient of expansion in the fine basis.

class gram_schmidt_fitting(object):
    
    def __init__(self,n):
        self.n=n
        self.ana_sed=analytic_sed.analytic_sed()
        self.ana_sed.create_fn_dir(self.n)
        self.ana_sed.fn_dir
        self.num_basis=0

    def gen_vectors(self,T,slope,nu_min,nu_max,sampling=1e5,nu=[],logspace=True):
        self.T=T
        self.slope=slope

        nu0=cnst.boltzman_const*self.T/cnst.planck_const/cnst.ghz2hz 
        c0=cnst.planck_const*cnst.ghz2hz/cnst.boltzman_const

        if np.size(nu)==0:
            if logspace:
                self.nu=np.logspace(np.log10(nu_min),np.log10(nu_max),sampling)
                self.nu[0]=nu_min ; self.nu[-1]=nu_max
            else:
                self.nu=np.linspace(nu_min,nu_max,sampling)
                self.nu[0]=nu_min ; self.nu[-1]=nu_max
                
        else:
            self.nu=nu

        self.vectors=list()
        for i in range(len(self.ana_sed.fn_dir)):
            v=self.ana_sed.fn_dir[i](self.nu,1./self.T,self.slope,nu0,c0)
            self.vectors.append(v)
        
    def gram_schmidt(self,vectors,min_vec_norm=1e-8):
        basis = []
        un_basis = []
        fewer_basis_than_vectors=False

        for ivec,v in enumerate(vectors):
            w = v - np.sum( np.dot(v,b)*b  for b in basis )
        
            #Calculate the norm for w-vector.
            norm=np.sqrt(np.dot(w,w))
            if ivec==0:
                norm0=norm
            #print norm/norm0

            #Normalize vector.
            # If one has sampled at only N frequencies, one cannot have more than N basis
            if (norm/norm0 > min_vec_norm):  
                un_basis.append(w)
                basis.append(w/norm)         
            else:
                fewer_basis_than_vectors=True

        if fewer_basis_than_vectors:
            print "Number of vectors:",len(vectors)
            print "Number of basis vectors:",len(basis)
            print "There are fewer basis functions than the number of vectors passed"

        return np.array(un_basis),np.array(basis)

    def gram_schmidt_iterative(self,tol=1e-8,iter=10,min_vec_norm=1e-6):
        vectors=self.vectors

        normalize=True ; counter=0
        while normalize:
            print counter
            un_basis,self.basis=self.gram_schmidt(vectors,min_vec_norm=min_vec_norm) ; normalize=False
        
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
        self.num_basis=np.shape(self.basis)[0]

    def get_gram_schmidt_param(self,nu,Inu,n,n_is_der_order=True):
    #def get_gram_schmidt_param(self,nu,Inu,n,n_is_der_order=True,rcond=1e-6):
        par=[]
        if n_is_der_order:
            num_par=self.ana_sed.calc_num_vec(n)
        else:
            num_par=n

        if self.num_basis<num_par:
            num_par=self.num_basis
            print "Have only", num_par, "basis function, so will fit only as many parameters."

        self.cov=np.zeros((num_par,num_par),np.float64)
        for i in range(num_par):
            calc_sparse_basis_i=interp1d(self.nu,self.basis[i],kind="cubic")
            sparse_basis_i=calc_sparse_basis_i(nu)
            for j in range(num_par):
                calc_sparse_basis_j=interp1d(self.nu,self.basis[j],kind="cubic")
                sparse_basis_j=calc_sparse_basis_j(nu)
                self.cov[i,j]=np.dot(sparse_basis_i,sparse_basis_j)
            par.append(np.dot(sparse_basis_i,Inu))

        # Computing the inverse of the covariance to return the coefficients of expansion.
        gs_par=np.asarray(np.matmul(np.linalg.inv(np.matrix(self.cov)),np.transpose(np.matrix(par)))).flatten()
        #gs_par=np.asarray(np.matmul(np.linalg.pinv(np.matrix(self.cov),rcond=rcond),np.transpose(np.matrix(par)))).flatten()
        return gs_par

    def reconstruct_sed(self,nu,*coeffs):
        sed=0.
        for i,par in enumerate(coeffs):
            calc_sparse_basis=interp1d(self.nu,self.basis[i],kind="cubic")
            sparse_basis=calc_sparse_basis(nu)
            sed=sed+par*sparse_basis
        return sed
        # It has to be the gram-schmidt basis function here and not the raw vectors.
