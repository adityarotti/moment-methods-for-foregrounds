import numpy as np
from scipy.interpolate import interp1d
import healpy as h
from matplotlib import pyplot as plt
import os

plt.ioff()

from matplotlib import rcParams,rc
params = {'backend': 'pdf',
          'savefig.dpi': 200,
          'axes.labelsize': 15,
          'axes.linewidth' : 2,
          'lines.linewidth' : 2,
          'font.size': 15,
          'xtick.labelsize': 15,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 15,
          'text.usetex': True,
          'font.family':'sans-serif',
          'font.sans-serif':'FreeSans'}
rc('text.latex', preamble='\usepackage{sfmath}')
rcParams.update(params)


class sim_methods(object):

    def __init__(self,npix,fmin,fmax,num_channel):

        self.h=6.62607004e-34 #Plancks constant #m2 kg / s
        self.k=1.38064852e-23 #Boltzmann constant #m2 kg s-2 K-1

        self.npix=npix
        self.fmin=fmin ; self.fmax=fmax ; self.num_channel=num_channel
        
        #self.nu=np.linspace(self.fmin,self.fmax,self.num_channel)        
        self.nu=np.logspace(np.log10(self.fmin),np.log10(self.fmax),self.num_channel)
        self.data=np.zeros((self.num_channel,self.npix),np.float64)
        self.nemit=np.zeros(self.npix,np.float64)
        self.mean_alpha=np.zeros(self.npix,np.float64)
        self.mean_T=np.zeros(self.npix,np.float64)

#   These function generate the frequency dependent spectra at each map pixel
    def mbb(self,nu,T,alpha,A,nu0=1.):
        x=self.h*nu*1e9/(self.k*T) 
        Inu=A*((nu/nu0)**alpha)*(nu**3.)/(np.exp(x)-1.)
        return Inu

    def integrate_mbb_spectra(self,nu,T_dist,alpha_dist,nu0=1.):
        Inu=np.zeros(np.size(nu),np.float64)
        for i in range(np.size(T_dist)):
            Inu=Inu + self.mbb(nu,T_dist[i],alpha_dist[i],A=1.,nu0=nu0)
        return Inu #/np.size(T_dist)

    def gen_mf_data(self,mu_T1,sig_T1,mu_T2,sig_T2,mu_alpha1,sig_alpha1,mu_alpha2,sig_alpha2,mu_N,sig_N,seed_init=0):
        for i in range(self.npix):


            np.random.seed(i+seed_init)
            # Sampling from a Chi^2 distribution
            if sig_N==0:
                num_emitters=mu_N
            else:
                num_emitters=int(max(1.,np.random.chisquare(mu_N,1)))
            self.nemit[i]=num_emitters
            num_emitters_1=int(np.floor(num_emitters/2.))
            num_emitters_2=int(np.ceil(num_emitters/2.))

            # Sampling the T for sources from a Gaussian distribution
            np.random.seed(i+seed_init+175464)

            if sig_T1==0.:
                T_dist1=np.ones(num_emitters_1,np.float64)*mu_T1
            else:
                T_dist1=np.random.normal(mu_T1,sig_T1,num_emitters_1) ; T_dist1[T_dist1<0]=mu_T1

            if sig_T2==0.:
                T_dist2=np.ones(num_emitters_2,np.float64)*mu_T2
            else:
                T_dist2=np.random.normal(mu_T2,sig_T2,num_emitters_2) ; T_dist2[T_dist2<0]=mu_T2

            T_dist=np.append(T_dist1,T_dist2) #; print T_dist
            self.mean_T[i]=np.sum(T_dist)/float(num_emitters)
                
            # Sampling the spectra slope alpha for sources from a Gaussian distribution
            np.random.seed(i+seed_init+188754)
            if sig_alpha1==0.:
                alpha_dist1=np.ones(num_emitters_1,np.float64)*mu_alpha1
            else:
                alpha_dist1=np.random.normal(mu_alpha1,sig_alpha1,num_emitters_1)

            if sig_alpha2==0.:
                alpha_dist2=np.ones(num_emitters_2,np.float64)*mu_alpha2
            else:
                alpha_dist2=np.random.normal(mu_alpha2,sig_alpha2,num_emitters_2)

            alpha_dist=np.append(alpha_dist1,alpha_dist2) #; print alpha_dist
            self.mean_alpha[i]=np.sum(alpha_dist)/float(num_emitters)

            self.data[:,i]=self.integrate_mbb_spectra(self.nu,T_dist,alpha_dist,nu0=1.)

    def normalize_data(self,A0,freq):
        monopole=np.zeros(self.num_channel,np.float64)
        for i in range(self.num_channel):
            monopole[i]=np.mean(self.data[i,:])

        fn=interp1d(self.nu,monopole,kind="quadratic")
        norm=A0/fn(freq) #; print norm

        for i in range(self.num_channel):
            self.data[i,:]=self.data[i,:]*norm


    def gen_data_plots(self,figpath,doall=False):
        datapath=figpath+"datafig/"
        if not os.path.exists(datapath):
            os.makedirs(datapath)

        h.mollview(self.nemit,min=np.min(self.nemit),max=np.max(self.nemit),title="Number of emitters")
        plt.savefig(datapath + "num_emitters.pdf",bbox_inches="tight",dpi=100)
        plt.close()

        h.mollview(self.mean_T,min=np.min(self.mean_T),max=np.max(self.mean_T),title="Mean T")
        plt.savefig(datapath + "mean_temperature.pdf",bbox_inches="tight",dpi=100)
        plt.close()
        
        h.mollview(self.mean_alpha,min=np.min(self.mean_alpha)-1e-5,max=np.max(self.mean_alpha)+1e-5,title=r"Mean slope $\alpha$")
        plt.savefig(datapath + "mean_slope.pdf",bbox_inches="tight",dpi=100)
        plt.close()

        if doall:
            for idx in range(self.num_channel):
                h.mollview(self.data[idx,:],title=str(round(self.nu[idx],0)) + " GHz")
                plt.savefig(datapath + "data" + str(i).zfill(3) + ".jpeg",bbox_inches="tight",dpi=100)
                plt.close()
        else: 
            x=[10,30,100,217,353,500,800,1200,1500,1800,2500,3000]
            for i,nu in enumerate(x):
                idx = (np.abs(self.nu - nu)).argmin()
                h.mollview(self.data[idx,:],title=str(round(self.nu[idx],0)) + " GHz")
                plt.savefig(datapath + "data" + str(i).zfill(3) + ".jpeg",bbox_inches="tight",dpi=100)
                plt.close()
    
        # This generates an animated gif image of the data at multi-frequencies.
        workdir=os.getcwd()
        os.chdir(datapath)
        os.system("convert -quality 99 -density 150 -delay 50 -loop 0 ./data*.jpeg data_animate.gif")
        os.system("rm *.jpeg")
        os.chdir(workdir)
