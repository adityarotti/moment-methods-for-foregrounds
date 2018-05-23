import sympy as sp
import numpy as np

class analytic_sed(object):

    def __init__(self):
        self.nu,self.beta,self.s,self.nu0,self.c0=sp.symbols("nu beta s nu0 c0")
        self.mbb=((self.nu/self.nu0)**self.s)*(self.nu**3.)/(sp.exp(self.c0*self.nu*self.beta)-1.)

    def mbb_der(self,der_beta=0,der_s=0):
        expr=self.mbb.diff(self.beta,der_beta)*self.mbb.diff(self.s,der_s)*(self.beta**der_beta)/self.mbb/(sp.factorial(der_beta)*sp.factorial(der_s))
        #expr=self.mbb.diff(self.beta,der_beta).diff*(self.s,der_s)*(self.beta**der_beta)/(sp.factorial(der_beta)*sp.factorial(der_s))
        return expr

    
    def create_fn_dir(self,n):
        '''This function creates a directory of functions\n The function in the directories are the various derivatives of the base SED \n 0 : MBB \n 1 : d(MBB)/ds \n 2 : beta*d(MBB)/dbeta \n 3 : d^2(MBB)/ds^2 \n 4 : beta*d^2(MBB)/ds dbeta \n 5 : (beta^2) d^2(MBB)/dbeta^2 \n ... \n ... '''

        self.fn_dir={}

        # Here we increase the derivate order in slope and Temperature together.
        self.num_of_vectors=0
        for der in range(n+1):
            for der_T in range(der+1):
                der_s=abs(der-der_T)
                #print counter,der_T,der_s,self.mbb_der(der_T,der_s)/self.mbb_der()
                fn=sp.lambdify((self.nu,self.beta,self.s,self.nu0,self.c0),self.mbb_der(der_T,der_s),modules="numpy")
                self.fn_dir[self.num_of_vectors]=fn
                self.num_of_vectors=self.num_of_vectors+1

    def create_sym_fn_dir(self,n):
        '''This function creates a directory of functions\n The function in the directories are the various derivatives of the base SED \n 0 : MBB \n 1 : d(MBB)/ds \n 2 : beta*d(MBB)/dbeta \n 3 : d^2(MBB)/ds^2 \n 4 : beta*d^2(MBB)/ds dbeta \n 5 : (beta^2) d^2(MBB)/dbeta^2 \n ... \n ... '''

        self.sym_fn_dir={}
        self.num_of_sym_vectors=0
        # Here we increase the derivate order in slope and Temperature together.
        for der in range(n+1):
            for der_T in range(der+1):
                der_s=abs(der-der_T)
                self.sym_fn_dir[self.num_of_sym_vectors]=self.mbb_der(der_T,der_s)
                self.num_of_sym_vectors=self.num_of_sym_vectors+1

    def calc_num_vec(self,n):
        return int(((n+1.)*(n+2.)/2.))
