# 4/6/12 B. D. Malone
# The purpose of this script is to test the convergence of various parameters
# to the QP corrections. 


import sys,os,glob
import matplotlib 
matplotlib.use("agg")
from matplotlib import rc
from matplotlib import rcParams
from pylab import *
from matplotlib.font_manager import FontProperties

rcParams.update({'font.size':14})

topval=4  # top of the valence band (hack, was 20, but need this to work)

thefiles=glob.glob('sigma.*')
# we want to order the files such that the files are looked at in 
# decreasing order of accuracy
orderfiles=[]
xarray=[]
for thefile in thefiles:
    number=int(thefile.split('.')[1])
    xarray.append(number)
    index=(number,thefile)
    orderfiles.append(index)

orderfiles.sort(reverse=True)
xarray.sort()  # this needs to be sorted smallest to lowest, because it gets
               # passed to matplotlib

globaldata=[]
# the global data array is made up of things that look like 
# [label,[eqp1,eqp2,eqp3.....],[eqp1,eqp2,eqp3....], for all iknum]
for item in orderfiles:
    kcounter=0
    localdata=[]
    thefilename=item[1]
    label=item[0]
    localdata.append(label)
    infile=open(thefilename,'r')
    while True:
        line=infile.readline()
        if line=='': # EOF
            break
        if line.find('band_index')!=-1:
            minband=int(line.split()[1])
            maxband=int(line.split()[2])
            numbands=maxband-minband+1
        if line.find('k =')!=-1 and line.find('spin =')!=-1:
            kcounter=kcounter+1
            vallist=[]
            # found a kpoint
            kindex=int(line.split()[7])   
            kx=float(line.split()[2])
            ky=float(line.split()[3])
            kz=float(line.split()[4])
            if abs(kx)<1E-10 and abs(ky)<1E-10 and abs(kz)<1E-10:
                GAMindex=int(kindex)
            if abs(kx)<1E-10 and abs(ky-0.5)<1E-10 and abs(kz-0.5)<1E-10:
                Lindex=int(kindex)
            if abs(kx)<1E-10 and abs(ky)<1E-10 and abs(kz-0.5)<1E-10:
                Xindex=int(kindex)
            infile.readline()  # just blank line
            infile.readline()  # just header stuff we don't care about
            for i in range(0,numbands):
                sigmaval=infile.readline().split()[9]
                vallist.append(float(sigmaval))
            localdata.append(vallist)
    infile.close()
    globaldata.append(localdata)
    totalk=kcounter
    
print('There are {0} points in the files'.format(totalk))

# okay we have all the data now. Now what we want to do is plot it, which 
# requires us to set up some value arrays to plot with. Basically for every
# band and every kpoint we have line to plot, so it'll need it's own
# subset of an array. 

fullyarray=[]   # all y values for plotting will be subsets of this array
                # we will simply go through each kpoint and do all bands
                # and then move onto the next

maxerror=-1000.

for ik in range(0,totalk):
    for ib in range(0,numbands):
        values=[]
        # we need to get this number for this band and kpoint from each
        # file and put it into this array, but we need to subtract off
        # the converged value which is in the first set of data (since 
        # the list is reverse ordered)
        # we take all of these values relative to the mode at $\Gamma$, since
        # we really probably just want to care about the topology of the 
        # bandstructure changing, not the absolute value
        for number in range(len(xarray)-1,-1,-1):
            convvalue=globaldata[0][ik+1][ib]-globaldata[0][GAMindex][topval-1]
            curvalue=globaldata[number][ik+1][ib]-globaldata[number][GAMindex][topval-1]
            values.append(abs(curvalue-convvalue)*1000)
            if 1000*abs(curvalue-convvalue)>maxerror:
                maxerror=abs(curvalue-convvalue)*1000
#            if number==(len(xarray)-1):
#                if abs(curvalue-convvalue)>0.03:
#                    print ik,ib
        fullyarray.append(values)

GAMerror=[]
Lerror=[]
Xerror=[]
dirLerror=[]
dirXerror=[]
for number in range(len(xarray)-1,-1,-1):
    convGAM=globaldata[0][GAMindex][topval]-globaldata[0][GAMindex][topval-1]
    curGAM=globaldata[number][GAMindex][topval]-globaldata[number][GAMindex][topval-1]
    GAMerror.append(abs(curGAM-convGAM)*1000)
    convL=globaldata[0][Lindex][topval]-globaldata[0][GAMindex][topval-1]
    curL=globaldata[number][Lindex][topval]-globaldata[number][GAMindex][topval-1]
    Lerror.append(abs(curL-convL)*1000)
    convX=globaldata[0][Xindex][topval]-globaldata[0][GAMindex][topval-1]
    curX=globaldata[number][Xindex][topval]-globaldata[number][GAMindex][topval-1]
    Xerror.append(abs(curX-convX)*1000)
    dirLconv=globaldata[0][Lindex][topval]-globaldata[0][Lindex][topval-1]
    dirLcur=globaldata[number][Lindex][topval]-globaldata[number][Lindex][topval-1]
    dirLerror.append(abs(dirLcur-dirLconv)*1000)
    dirXconv=globaldata[0][Xindex][topval]-globaldata[0][Xindex][topval-1]
    dirXcur=globaldata[number][Xindex][topval]-globaldata[number][Xindex][topval-1]
    dirXerror.append(abs(dirXconv-dirXcur)*1000)


# plot
print len(fullyarray)
for i in range(0,len(fullyarray)):
    plot(xarray,fullyarray[i],'y-s')

plot(xarray,GAMerror,'k-o',label=r'$\Gamma-\Gamma$')
plot(xarray,Lerror,'b-v',label=r'$\Gamma-L$')
plot(xarray,Xerror,'g-s',label=r'$\Gamma-X$')
plot(xarray,dirLerror,'c-p',label=r'$L-L$')
plot(xarray,dirXerror,'r-D',label=r'$X-X$')

ylabel('Energy error (meV)')
#xlabel(r'# states in $\Sigma_{CH}$')
#xlabel(r'# states in $\epsilon$')
xlabel(r'$\Sigma$ cutoff')
plt.axhline(y=20.0,color='k',linewidth=2.0,alpha=0.3)
axis([xarray[0],xarray[-1],0,100])
legend(prop=FontProperties(size='small'))
plt.savefig("Sig-cutoff_conv.png")
show()

