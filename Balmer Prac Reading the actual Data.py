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
from os import listdir
from os.path import isfile, join
## -------------------------------------------------------------------------------------------------------------------------------------------------------------

## Reading in the files ----------------------------------------------------------------------------------------------------------------------------------------
mypath ='C:\\Users\\Austin Connolly\\OneDrive - University of Cape Town\\University Work\\2021\\PHY3004W\\Labs\\Balmer Prac\\Hydrogen Data\\' # Path of files

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))] # Gets the names of the files in the path
o=0
c=3
filenm = []
for i in range(len(onlyfiles)):
    num='H'+str(c)
    filenam = onlyfiles[i]
    if filenam.find(num)==0:
        filenm.append(onlyfiles[i])
        o=o+1
        if o==2:
            o=0
            c=c+1

path1 ='C:\\Users\\Austin Connolly\\OneDrive - University of Cape Town\\University Work\\2021\\PHY3004W\\Labs\\Balmer Prac\\Hydrogen Data\\' # Path of files

# Names of files

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
    #print(udata)
    #print(len(ydata))
    #print(len(tdata))
    
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
    tmodel = np.linspace(ends[c], ends[c+1], 10000)
    
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
    #plt.title(filenm[h])
    plt.plot(tmodel,yfit,'-r', label='Line of best fit') # Plots 
    xmin,xmax,ymin,ymax = plt.axis([ends[0]-1,ends[-1]+1,0,np.max(yfit)+100])
    #plt.plot(f(waveln[rr[0]:rr[-1]],*p0),counts[rr[0]:rr[-1]])
    plt.grid()
    plt.legend()
    plt.tick_params(direction='in',top=True,right=True) # Assigns lines to be on the inside
    
    plt.show()    
## -------------------------------------------------------------------------------------------------------------------------------------------------------------

## Changing to meters ------------------------------------------------------------------------------------------------------------------------------------------
for i in range(len(mus)):
    mus[i]=(mus[i])*(10**(-10))
    sigmas[i] = (sigmas[i])*(10**(-10))

## -------------------------------------------------------------------------------------------------------------------------------------------------------------

## Getting the calibration value -------------------------------------------------------------------------------------------------------------------------------

print(mus)
print(sigmas)

newmus=[]
newsigmas=[]
c=0
for i in range(int(len(mus)/2)):
    newmus.append((mus[c]+mus[c+1])/2-(22.36)*(10**(-10)))
    newsigmas.append((np.sqrt((np.sqrt(sigmas[c]**2+sigmas[c+1]**2)/2)**2+((0.43)*(10**(-10)))**2)))
    c+=2
print(newmus)
print(newsigmas)

actval = [6562.85175*10**(-10),4861.283363*10**(-10),4340.472*10**(-10),4101.707462*10**(-10)]
u_actval = [0.00007*10**(-10),0.000024*10**(-10),0.006*10**(-10),0.000021*10**(-10)]

df1 = pd.DataFrame({'Wavelength (m)':newmus,'Uncertainty of Wavelength (m)':newsigmas, 'Actual Wavelength values (m)':actval,'Uncertainty of Wavelength for Actual (m)':u_actval})

df1.to_excel(path1+'Table with H data wavelengths.xlsx',index=False)

df1.to_csv(path1+'Table with H data wavelengths.csv',index=False)

## -------------------------------------------------------------------------------------------------------------------------------------------------------------

## Weighted Linear Fit -----------------------------------------------------------------------------------------------------------------------------------------

ms = [3,4,5,6]

xs = []
ys = []
us = []

for i in range(len(ms)):
    xs.append(((1/(2)**2)-(1/(ms[i])**2)))
    ys.append(1/newmus[i])
    us.append((newsigmas[i]/newmus[i])*(1/newmus[i]))


print(xs,'\n',ys,'\n',us)


i=0

xaxis=xs
yaxis=ys
uncertainty=us


N=len(yaxis)

udata=np.zeros(N)
ydata=np.zeros(N)
xdata=np.zeros(N)
tdata=np.zeros(N)

while i<N:
    udata[i]=uncertainty[i]
    ydata[i]=yaxis[i]
    tdata[i]=xaxis[i]
    i=i+1    


delta=np.sum(1/(udata**2))*np.sum((tdata**2)/(udata**2))-(np.sum((tdata)/(udata**2)))**2

grad=(np.sum(1/(udata**2))*np.sum((tdata*ydata)/(udata**2))-np.sum(tdata/(udata**2))*np.sum(ydata/(udata**2)))/delta

u_grad = np.sqrt(((np.sum(1/(udata**2)))/delta))

cval = ((np.sum((tdata**2)/(udata**2)))*(np.sum((ydata)/(udata**2)))-(np.sum((tdata)/(udata**2)))*(np.sum((tdata*ydata)/(udata**2))))/delta

u_cval = np.sqrt((np.sum((tdata**2)/(udata**2)))/delta)

print(delta)
print(grad,"+/-",u_grad)
print('10973731.6')

ts=np.linspace(xaxis[0],xaxis[-1],1000)
plt.scatter(xaxis,yaxis,label='Data Points')
#plt.errorbar(tdata, ydata, xerr = None, yerr = udata, fmt = '', marker='.', ls = 'None',capsize=2.3, ecolor = 'b',label='Data points')
plt.plot(ts,np.asarray(grad)*ts+cval,color='r',label='Weighted line of best fit')
plt.tick_params(direction='in',top=True,right=True)
plt.xlabel('$1/4 - 1/m^2$')
plt.ylabel('$1/\\lambda (m^{-1})$')
plt.grid()
plt.legend()
plt.show()