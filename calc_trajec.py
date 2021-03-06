#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from math import factorial
import numpy as np
import sys
import matplotlib.pyplot as plt
import glob



def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def process_single_file(fname,head_name,oldav,step,name):
    
    fin = open(fname)
    count = 0
    natoms = 0
    col    = -1
    av = 0.0
    i_boxlow = 6
    xc = []
    
    if head_name == "y": i_boxlow += 1
    if head_name == "z": i_boxlow += 2
    for line in fin:
        count += 1
        if count == 4: 
            natoms = int(line)
            continue
        if count == i_boxlow:
            line = line.strip().split()
            lower = float(line[0])
            upper = float(line[1])
            boxlength = upper-lower
            boxcenter = (upper+lower)/2.0
        if count == 9: 
            header = line.strip().split()
            col = header.index(head_name)-2
            continue
        if count > 9:
            tmp = float(line.strip().split()[col])
            xc.append(tmp)
                
    av = np.average(xc)
        
    fin.close()
    return av

pers = ArgumentParser()

# Add options
pers.add_argument("fin",
                  help="input file", metavar="FILE")
pers.add_argument("prop",
                  help="data header name",  default="x", metavar="PROPERTY")
pers.add_argument("units",
                  help="data header units", default="A", metavar="UNITS")

pers.add_argument("-r", "--range",
                  help="multiple files", nargs=2, type=int, metavar=("START", "END"))
pers.add_argument("-o","--output",
                  help="output data", type=str, default="traj.dat")
pers.add_argument("--dt",
                  help="timestep", nargs=2, metavar=("VALUE","UNITS"))
pers.add_argument("--cont","-c",help="Make trajectory continuos",default=False,
                  action="store_true")

args = pers.parse_args()
print args
fout = open(args.output, "w")
if args.dt is None:
    tunits = "step"
else:
    tunits = args.dt[1]

out_header = "#time/step"+"("+tunits+") "     #col1
out_header += args.prop+"("+args.units+") "   #col2
out_header += "ysmooth" +"("+args.units+") "  #col3
out_header += " v("+args.units+"/"+tunits+")" #col4
fout.write(out_header)

print(args)
iwild = args.fin.find("*")
if iwild == -1:
    print("Warning: ONLY ONE FILE SPECIFIED")
    args.range = [0,1,1]
    keys = [args.fin]
else:    
    flist = glob.glob(args.fin) #list of filenames
    print("... %1.0d files found"%(len(flist)))
    print(iwild+1)
    if iwild+1==len(args.fin): iend = None
    else: iend = s.find(args.fin[iwild+1:])
    keys = [int( s[iwild:iend] ) for s in flist]
    keys.sort()
    flist = sorted(flist, key=lambda s: int( s[iwild:iend] )) 

step = -1
load = ["\\","|","/","-"]

if args.range is not None:
    N = args.range[1]-args.range[0]+1
    flist = flist[args.range[0]:args.range[1]+1]
else:
    N = len(flist)
print("... using %i files from original"%(N))
y = []

t = 0
tmp = None

for i in range(N):
    fname = flist[i]
    tmp = process_single_file(fname,args.prop,tmp,i,fname)
    y.append(tmp)
    t += 1
    prcnt = float(t)/N*100
    #print(load[i%4]+ "     %d %%"%(prcnt))
    #sys.stdout.write("\033[F")

y = np.array(y)

if args.cont:
    y -= y[0]
    for iy in range(1,len(y)):
        if y[iy]>y[iy-1]+20:
            y[iy:] -= y[iy]-y[iy-1]
        elif y[iy]<y[iy-1]-20:
            y[iy:] += y[iy-1]-y[iy]
        

#print y
#print(np.amin(y))
#print(np.argmin(y))
#y[np.argmin(y)-1:] -= np.amin(y)
#print y
ysmooth = savitzky_golay(y,31,3,deriv=0)
#dy = savitzky_golay(y,31,3,deriv=1)/float(args.dt[0])
dy = np.diff(ysmooth)/float(args.dt[0])
    
for i in range(N):
    if args.dt is None:
        fout.write("\n%1.0d %1.6f %1.6f"%(i,y[i],dy[i]))
    else:
        fout.write("\n%1.6f %1.6f %1.6f %1.6f"%(i*float(args.dt[0]), y[i], ysmooth[i], dy[1]))
print("... done writing to "+args.output)
        
        
        
    

