import healpy as h
import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
#from matplotlib.ticker import FormatStrFormatter
import time

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

class mbb_moment_fits(object):
        
    def __init__(self,npix,channels,data,figpath):
        self.nu=channels # In GHz
        self.npix=npix
        self.data=data

        self.h=6.62607004e-34 #Plancks constant #m2 kg / s
        self.k=1.38064852e-23 #Boltzmann constant #m2 kg s-2 K-1s
        self.T=[] ; self.alpha=[] ; self.A=[]
        
        self.figpath=figpath


    def get_data_monopole(self):
        self.monopole=np.zeros(np.size(self.nu),float)
        for i in range(np.size(self.nu)):
            self.monopole[i]=np.mean(self.data[i,:])

#   Fitting the Global monopole.
    def fit_mbb_monopole(self,nu,Inu,prange=2.,guess=[30.,0.,1e-6],lb=[0.,-3,0.],ub=[100.,3.,10.],maxfev=20000,makeplot=True):
        guess=np.append(guess,np.zeros(7,float))
        lb=np.append(lb,np.ones(7,float)*(-prange)) ; ub=np.append(ub,np.ones(7,float)*(prange))

        self.par_monopole=np.zeros(np.size(guess),float)

        err=Inu

        fit_param,fit_param_cov=curve_fit(self.mbb_moments,nu,Inu,p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
        self.T=fit_param[0] ; self.alpha=fit_param[1] ; self.A=fit_param[2]
        #print fit_param
        print "T =",self.T, r" alpha=",self.alpha, " A=",self.A

        self.par_monopole[:]=fit_param[:]        

        self.cumulative_err_order0()

        if makeplot:
            self.make_monopole_fit_plot()

#   Fitting the spatial moments, around the Global monopole fits.
    def fit_mbb_order1(self,guess=[0.,0.,0.],prange=2.,maxfev=20000,makeplot=True):
        self.par_o1=np.zeros((np.size(guess),self.npix),float)
        lb=np.ones(np.size(guess),float)*(-prange) ; ub=np.ones(np.size(guess),float)*(prange)
        #lb[0]=lb[0]*50 ; ub[0]=ub[0]*50

        start=time.clock()
        for i in range(self.npix):
            err=self.data[:,i]
            fit_param,fit_param_cov=curve_fit(self.mbb_order1,self.nu,self.data[:,i],p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
            self.par_o1[:,i]=fit_param[:]
        print "This took:", (time.clock()-start)/self.npix, "seconds per pixel."

        self.cumulative_err_order1()

        if makeplot:
            self.make_plot_order1()

    def fit_mbb_order2(self,guess=[0.,0.,0.,0.,0.,0.],prange=2.,maxfev=20000,makeplot=True):
        self.par_o2=np.zeros((np.size(guess),self.npix),float)
        lb=np.ones(np.size(guess),float)*(-prange) ; ub=np.ones(np.size(guess),float)*(prange)
        #lb[0]=lb[0]*50 ; ub[0]=ub[0]*50

        start=time.clock()
        for i in range(self.npix):
            err=self.data[:,i]
            fit_param,fit_param_cov=curve_fit(self.mbb_order2,self.nu,self.data[:,i],p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
            self.par_o2[:,i]=fit_param[:]
        print "This took:", (time.clock()-start)/self.npix, "seconds per pixel."

        self.cumulative_err_order2()
        
        if makeplot:
            self.make_plot_order2()

    def fit_mbb_order3(self,guess=[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.],prange=2.,maxfev=20000,makeplot=True):
        self.par_o3=np.zeros((np.size(guess),self.npix),float)
        lb=np.ones(np.size(guess),float)*(-prange) ; ub=np.ones(np.size(guess),float)*(prange)
        #lb[0]=lb[0]*50 ; ub[0]=ub[0]*50

        start=time.clock()
        for i in range(self.npix):
            err=self.data[:,i]
            fit_param,fit_param_cov=curve_fit(self.mbb_order3,self.nu,self.data[:,i],p0=guess,sigma=err,bounds=(lb,ub),maxfev=maxfev)
            self.par_o3[:,i]=fit_param[:]
        print "This took:", (time.clock()-start)/self.npix, "seconds per pixel."

        self.cumulative_err_order3()
 
        if makeplot:
            self.make_plot_order3()

#   Fitting functions
    def I0(self,nu,T,alpha,A):
        x=self.h*nu*1e9/self.k/T
        nu0=self.k*T/self.h/1e9
        return A*(((nu/nu0)**alpha)*(nu**3.))/(np.exp(x)-1.)

    def mbb_moments(self,nu,T,alpha,A,p22,p23,p33,p222,p223,p233,p333):
        x=self.h*nu*1e9/self.k/T
        nu0=self.k*T/self.h/1e9
        I0=A*(((nu/nu0)**alpha)*(nu**3.))/(np.exp(x)-1.)
        temp=I0*(1. + 0.5*p22*(np.log(nu/nu0))**2. + p23*np.log(nu/nu0)*self.y1(x)+ 0.5*p33*self.y2(x))
        temp = temp +  I0*(p222*(np.log(nu/nu0)**3.)/6. + 0.5*p223*(np.log(nu/nu0)**2.)*self.y1(x) + 0.5*p233*np.log(nu/nu0)*self.y2(x) + p333*self.y3(x)/6.)
        return temp

    # Functions for fitting around the Global fit.
    def mbb_order1(self,nu,p0,p2,p3):
        x=self.h*nu*1e9/self.k/self.T
        nu0=self.k*self.T/self.h/1e9
        I0=self.A*(((nu/nu0)**self.alpha)*(nu**3.))/(np.exp(x)-1.)
        temp=I0*(1. + p0 + p2*(np.log(nu/nu0)) + p3*self.y1(x))
        return temp

    def mbb_order2(self,nu,p0,p2,p3,p22,p23,p33):
        x=self.h*nu*1e9/self.k/self.T
        nu0=self.k*self.T/self.h/1e9
        I0=self.A*(((nu/nu0)**self.alpha)*(nu**3.))/(np.exp(x)-1.)
        temp=I0*(1. + p0 + p2*(np.log(nu/nu0)) + p3*self.y1(x))
        temp=temp + I0*(0.5*p22*(np.log(nu/nu0))**2. + p23*np.log(nu/nu0)*self.y1(x)+ 0.5*p33*self.y2(x))
        return temp

    def mbb_order3(self,nu,p0,p2,p3,p22,p23,p33,p222,p223,p233,p333):
        x=self.h*nu*1e9/self.k/self.T
        nu0=self.k*self.T/self.h/1e9
        I0=self.A*(((nu/nu0)**self.alpha)*(nu**3.))/(np.exp(x)-1.)
        temp=I0*(1. + p0 + p2*(np.log(nu/nu0)) + p3*self.y1(x))
        temp=temp + I0*(0.5*p22*(np.log(nu/nu0))**2. + p23*np.log(nu/nu0)*self.y1(x)+ 0.5*p33*self.y2(x))
        temp = temp +  I0*(p222*(np.log(nu/nu0)**3.)/6. + 0.5*p223*(np.log(nu/nu0)**2.)*self.y1(x) + 0.5*p233*np.log(nu/nu0)*self.y2(x) + p333*self.y3(x)/6.)
        return temp
       
# Definition of the derivates of the Planckian function
    def y1(self,x):
        return x*np.exp(x)/(np.exp(x)-1.)

    def y2(self,x):
        return self.y1(x)*x*np.cosh(0.5*x)/np.sinh(0.5*x)

    def y3(self,x):
        return self.y1(x)*x*x*(np.cosh(x)+2)/(np.cosh(x)-1)

    def y4(self,x):
        return self.y2(x)*0.5*x*x*(np.cosh(x)+5)/(np.sinh(0.5*x)**2.)

    def y5(self,x):
        return self.y1(x)*(x**4.)*(33.+26*np.cosh(x)+np.cosh(2.*x))/(8.*np.sinh(0.5*x)**4.)


#   Evaluating the cumulative error maps.
    def cumulative_err_order0(self):
        par=self.par_monopole
        self.cerr_o0=np.zeros(self.npix,float)
        for i in range(self.npix):
            self.cerr_o0[i]=np.sum(abs(self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9])-self.data[:,i])/self.data[:,i])/np.size(self.nu)

    def cumulative_err_order1(self):
        self.cerr_o1=np.zeros(self.npix,float)
        for i in range(self.npix):
            par=self.par_o1[:,i]
            y=self.mbb_order1(self.nu,p0=par[0],p2=par[1],p3=par[2])
            self.cerr_o1[i]=np.sum(abs(y[:]-self.data[:,i])/self.data[:,i])/np.size(self.nu)

    def cumulative_err_order2(self):
        self.cerr_o2=np.zeros(self.npix,float)
        for i in range(self.npix):
            par=self.par_o2[:,i]
            y=self.mbb_order2(self.nu,p0=par[0],p2=par[1],p3=par[2],p22=par[3],p23=par[4],p33=par[5])
            self.cerr_o2[i]=np.sum(abs(y[:]-self.data[:,i])/self.data[:,i])/np.size(self.nu)

    def cumulative_err_order3(self):
        self.cerr_o3=np.zeros(self.npix,float)
        for i in range(self.npix):
            par=self.par_o3[:,i]
            y=self.mbb_order3(self.nu,p0=par[0],p2=par[1],p3=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9])
            self.cerr_o3[i]=np.sum(abs(y[:]-self.data[:,i])/self.data[:,i])/np.size(self.nu)


    def make_plot_order1(self):
        #plt.figure()
        f, (ax1, ax2) = plt.subplots(2, sharex=True) ; plt.subplots_adjust(hspace=0.1)
        
        for i in range(self.npix):
            par=self.par_o1[:,i]
            ax1.plot(self.nu,self.data[:,i],"k-",lw=2,alpha=0.2)
            ax1.plot(self.nu,abs(self.data[:,i]-self.mbb_order1(self.nu,p0=par[0],p2=par[1],p3=par[2])),"c-",lw=2,alpha=0.2)
            ax2.plot(self.nu,abs(self.data[:,i]-self.mbb_order1(self.nu,p0=par[0],p2=par[1],p3=par[2]))/self.data[:,i],"c-",lw=2,alpha=0.2)
        
        par=self.par_monopole
        ax2.plot(self.nu,abs(self.monopole-self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]))/self.monopole,"r--",lw=2)
        ax1.plot(self.nu,self.monopole,"y-",lw=2,label="Data monopole")
        ax1.plot(self.nu,self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]),"m--",lw=2,label="Fit to monopole")
        
        ax1.loglog() ; ax1.set_ylim(1e-6,5e9) ; ax1.grid(which="both",alpha=0.4)
        ax2.loglog() ; ax2.set_ylim(1e-10,1e1) ; ax2.grid(which="both",alpha=0.4)
        ax2.set_yticks([1e0,1e-2,1e-4,1e-6,1e-8])
        #; ax2.set_yticklabels([1e0,1e-2,1e-4,1e-6,1e-8])
        #ax2.yaxis.set_major_formatter(FormatStrFormatter('%0.1e'))
        
        
        ax1.set_ylabel(r"$\Delta I_{\nu}$")
        ax2.set_ylabel(r"$\Delta I_{\nu}/I_{\nu}$")
        ax2.set_xlabel("Frequency GHz")
        ax1.set_title("1 Order fits")
        ax1.legend(loc=0,fontsize=10)
        plt.savefig(self.figpath + "assess_1order_fits.pdf",dpi=150,bbox_inches="tight")
        plt.close()


    def make_plot_order2(self):
        #plt.figure()
        f, (ax1, ax2) = plt.subplots(2, sharex=True) ; plt.subplots_adjust(hspace=0.1)
        
        for i in range(self.npix):
            par=self.par_o2[:,i]
            ax1.plot(self.nu,self.data[:,i],"k-",lw=2,alpha=0.2)
            ax1.plot(self.nu,abs(self.data[:,i]-self.mbb_order2(self.nu,p0=par[0],p2=par[1],p3=par[2],p22=par[3],p23=par[4],p33=par[5])),"b-",lw=2,alpha=0.2)
            ax2.plot(self.nu,abs(self.data[:,i]-self.mbb_order2(self.nu,p0=par[0],p2=par[1],p3=par[2],p22=par[3],p23=par[4],p33=par[5]))/self.data[:,i],"b-",lw=2,alpha=0.2)
            
        par=self.par_monopole
        ax2.plot(self.nu,abs(self.monopole-self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]))/self.monopole,"r--",lw=2)
        ax1.plot(self.nu,self.monopole,"y-",lw=2,label="Data monopole")
        ax1.plot(self.nu,self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]),"m--",lw=2,label="Fit to monopole")
        
        ax1.loglog() ; ax1.set_ylim(1e-6,5e9) ; ax1.grid(which="both",alpha=0.4)
        ax2.loglog() ; ax2.set_ylim(1e-10,1e1) ; ax2.grid(which="both",alpha=0.4)
        ax2.set_yticks([1e0,1e-2,1e-4,1e-6,1e-8])

        ax1.set_ylabel(r"$\Delta I_{\nu}$")
        ax2.set_ylabel(r"$\Delta I_{\nu}/I_{\nu}$")
        ax1.set_title("2 Order fits")
        ax2.set_xlabel("Frequency GHz")
        ax1.legend(loc=0,fontsize=10)
        plt.savefig(self.figpath + "assess_2order_fits.pdf",dpi=150,bbox_inches="tight")
        plt.close()


    def make_plot_order3(self):
        #plt.figure()
        f, (ax1, ax2) = plt.subplots(2, sharex=True) ; plt.subplots_adjust(hspace=0.1)
        
        for i in range(self.npix):
            par=self.par_o3[:,i]
            ax1.plot(self.nu,self.data[:,i],"k-",lw=2,alpha=0.2)
            ax1.plot(self.nu,abs(self.data[:,i]-self.mbb_order3(self.nu,p0=par[0],p2=par[1],p3=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9])),"g-",lw=2,alpha=0.2)
            ax2.plot(self.nu,abs(self.data[:,i]-self.mbb_order3(self.nu,p0=par[0],p2=par[1],p3=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]))/self.data[:,i],"g-",lw=2,alpha=0.2)

        par=self.par_monopole
        ax2.plot(self.nu,abs(self.monopole-self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]))/self.monopole,"r--",lw=2)
        ax1.plot(self.nu,self.monopole,"y-",lw=2,label="Data monopole")
        ax1.plot(self.nu,self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]),"m--",lw=2,label="Fit to monopole")

        ax1.loglog() ; ax1.set_ylim(1e-6,5e9) ; ax1.grid(which="both",alpha=0.4)
        ax2.loglog() ; ax2.set_ylim(1e-10,1e1) ; ax2.grid(which="both",alpha=0.4)
        ax2.set_yticks([1e0,1e-2,1e-4,1e-6,1e-8])

        ax1.set_ylabel(r"$\Delta I_{\nu}$")
        ax2.set_ylabel(r"$\Delta I_{\nu}/I_{\nu}$")
        ax1.set_title("3 Order fits")
        ax2.set_xlabel("Frequency GHz")
        ax1.legend(loc=0,fontsize=10)
        plt.savefig(self.figpath + "assess_3order_fits.pdf",dpi=150,bbox_inches="tight")
        plt.close()

        h.mollview(self.par_o3[0,:],title=r"$P_{A}$") ; plt.savefig(self.figpath + "p0_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[1,:],title=r"$P_{\alpha}$") ; plt.savefig(self.figpath + "p2_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[2,:],title=r"$P_{T}$") ; plt.savefig(self.figpath + "p3_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[3,:],title=r"$P_{\alpha \alpha}$") ; plt.savefig(self.figpath + "p22_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[4,:],title=r"$P_{\alpha T}$") ; plt.savefig(self.figpath + "p23_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[5,:],title=r"$P_{T T}$") ; plt.savefig(self.figpath + "p33_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[6,:],title=r"$P_{\alpha \alpha \alpha}$") ; plt.savefig(self.figpath + "p222_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[7,:],title=r"$P_{\alpha \alpha T}$") ; plt.savefig(self.figpath + "p223_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[8,:],title=r"$P_{\alpha T T}$") ; plt.savefig(self.figpath + "p233_moment.pdf",dpi=150,bbox_inches="tight")
        h.mollview(self.par_o3[9,:],title=r"$P_{T T T}$") ; plt.savefig(self.figpath + "p333_moment.pdf",dpi=150,bbox_inches="tight")


    def make_monopole_fit_plot(self):
        f, (ax1, ax2) = plt.subplots(2, sharex=True) ; plt.subplots_adjust(hspace=0.1)
        
        for i in range(self.npix):
            ax1.plot(self.nu,self.data[:,i],"k-",lw=2,alpha=0.2)

        ax1.plot(self.nu,self.monopole,"y-",lw=2,label="Data monopole")
        par=self.par_monopole
        ax1.plot(self.nu,self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]),"m--",lw=2,label="Fit to monopole")
        ax1.plot(self.nu,abs(self.monopole-self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9])),"r--",lw=2,label="Difference")
        ax2.plot(self.nu,abs(self.monopole-self.mbb_moments(self.nu,T=par[0],alpha=par[1],A=par[2],p22=par[3],p23=par[4],p33=par[5],p222=par[6],p223=par[7],p233=par[8],p333=par[9]))/self.monopole,"r--",lw=2)

        ax1.loglog() ; ax1.set_ylim(1e-6,5e9) ; ax1.grid(which="both",alpha=0.4)
        ax2.loglog() ; ax2.set_ylim(1e-10,1e1) ; ax2.grid(which="both",alpha=0.4)
        ax2.set_yticks([1e0,1e-2,1e-4,1e-6,1e-8])

        ax1.set_ylabel(r"$\Delta I_{\nu}$")
        ax2.set_ylabel(r"$\Delta I_{\nu}/I_{\nu}$")
        ax1.set_title("3 Order fit to Monopole")
        ax2.set_xlabel("Frequency GHz")
        ax1.legend(loc=0,fontsize=10)
        plt.savefig(self.figpath + "monopole_fits.pdf",dpi=150,bbox_inches="tight")
        plt.close()
