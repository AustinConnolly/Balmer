# Balmer Prac Calibrations
# Austin Connolly
# 18/8/2021

## Importing libraries -----------------------------------------------------------------------------------------------------------------------------------------

from matplotlib import rc
import numpy as np
import scipy as sp
import scipy.stats as stats
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.special import gamma
import pandas as pd
## -------------------------------------------------------------------------------------------------------------------------------------------------------------

## Reading in the files ----------------------------------------------------------------------------------------------------------------------------------------

path1 ='C:\\Users\Austin Connolly\\OneDrive - University of Cape Town\\University Work\\2021\\PHY3004W\\Labs\\Balmer Prac\\' # Path of files
filenm = ['015um-100ms-6334-6375-0.2A.csv','020um-100ms-6332-6375-0.5A.csv','025um-100ms-6271-6425-1A.csv'] # Names of files

# Stuff in for loop --------------------------------------------------------------------------------------------------------------------------------------------
sigmas=[]
mus=[]
#h=0
c=0
# --------------------------------------------------------------------------------------------------------------------------------------------------------------

## Looping through all files -----------------------------------------------------------------------------------------------------------------------------------
for h in range(len(filenm)):
    
    nm = path1 + filenm[h] # path and name of the file
        
    data = pd.read_csv(nm,header = None,skiprows=13).values # Collects the values as a 2-D array
    
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------
    endps = []
    waveln = data[:,0]
    counts = data[:,1]
    N = len(counts)
    rr=[]
    for yy in range(len(waveln)): # starts the counting after roughly 200 hits
        if counts[yy] >=0:
            endps.append(waveln[yy])
            rr.append(yy)
    ends=[endps[0],endps[-1]] # ends of the Gaussian for each data set
    
    print(ends)
    
    udata=[]
    ydata=[]
    xdata=[]
    tdata=[]
    
    c=0
    g=0
    while g<N:
        if waveln[g]>=ends[c] and waveln[g]<=ends[c+1]:
            udata.append(np.sqrt(counts[g]))
            ydata.append(counts[g])
            tdata.append(waveln[g])
        if waveln[g]==ends[c+1]:
            break
        g+=1
    print(udata)
    print(len(ydata))
    print(len(tdata))
    
    for uu in range(len(ydata)):
        if ydata[uu]==np.max(counts):
            muval = tdata[uu]
    
    def f(x, A, mu, sigma, D):
        return A*np.exp(-((x-mu)**2)/(2*sigma**2)) + D
    
    A0 = np.max(counts)
    mu0 = muval  #(ends[c+1]+ends[c])/2
    sigma0 = (ends[c+1]-ends[c])/(0.40*len(tdata))    #(2*np.sqrt(2))
    D0 = ydata[0]
    p0 = [A0, mu0, sigma0, D0]
    name = ['A', 'mu', 'sigma', 'Vertical Shift']
    
    print(sigma0)
    print(mu0)
    tmodel = np.linspace(ends[c], ends[c+1], 1000)
    
    ystart=f(tmodel,*p0)
    
    
    ## Curve Fit Function --------------------------------------------------------------------------------------------------------------------------------------
    
    popt,pcov=curve_fit(f,tdata,ydata,p0,sigma=udata,absolute_sigma=True)
    dymin = (ydata-f(tdata,*popt))/udata
    min_chisq = sum(dymin*dymin)
    dof=len(tdata) - len(popt)
    
    print('Chi square: ',min_chisq)
    print('Number of degrees of freedom: ',dof)
    print('Chi square per degree of freedom: ',min_chisq/dof)
    print()
    
    mn=min_chisq / dof
    
    print('Fitted parameters with 68% C.I.: ')
    
    for i, pmin in enumerate(popt):
        print('%2i %-10s %12f +/- %10f'%(i, name[i], pmin, np.sqrt(pcov[i,i])*np.sqrt(mn)))
    
    print()
    
    perr=np.sqrt(np.diag(pcov))
    print(perr)
    
    print("Correlation matrix")
    print("               ")
    
    for i in range(len(popt)): print('%-10s'%(name[i],)),
    print()
    
    for i in range(len(popt)):
        print('%10s'%(name[i])),
        for j in range(i+1):
            print('%10f'%(pcov[i,j]/np.sqrt(pcov[i,i]*pcov[j,j]),)),
        print()
    
    mus.append(popt[1])
    sigmas.append(popt[2])
    yfit=f(tmodel,*popt)
    
    
    ## Plotting ------------------------------------------------------------------------------------------------------------------------------------------------
    plt.scatter(tdata,ydata)
    plt.errorbar(tdata, ydata, xerr = None, yerr = udata, fmt = '', marker='.', ls = 'None',capsize=2.3, ecolor = 'b',label='Data points') # Plots errorbar
    plt.xlabel('Wavelength (A)') # Labels x-axis
    #plt.plot(tmodel,ystart,'-g') # Plots estimate line
    plt.ylabel('Counts') # Labels y-axis    
    xmin,xmax,ymin,ymax = plt.axis([ends[0]-1,ends[-1]+1,0,np.max(yfit)+100])
    plt.plot(tmodel,yfit,'-r', label='Line of best fit') # Plots 
    plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside
    plt.grid()
    plt.legend()
    plt.show()    
    ## ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Getting the calibration value -------------------------------------------------------------------------------------------------------------------------------

waveval = np.sum(mus)/3
u_wav = np.sqrt((sigmas[0])**2+(sigmas[1])**2+(sigmas[2])**2)/3
print(waveval,'+/-',u_wav,'A')
Difval = waveval - 6328.1646
u_wav = np.sqrt(u_wav**2+(0.0004)**2)
print('The Calibration factor is: ',round(Difval,2),'+/-',round(u_wav,2),'A')
## -------------------------------------------------------------------------------------------------------------------------------------------------------------