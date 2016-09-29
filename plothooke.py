#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from math import factorial
import numpy as np
import sys
import matplotlib.pyplot as plt

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


persy = ArgumentParser()

# Add more options
persy.add_argument("fin",help="input file", metavar="FILE")
persy.add_argument("x_col",help="column 1-# to plot on x-axis", type=int, metavar="X")
persy.add_argument("y_col",help="column 1-# to plot on y-axis", type=int, metavar="Y") 

persy.add_argument("-l", "--label", help="label for data (no spaces)", metavar="LABEL", default=" ")

persy.add_argument("-p","--png", dest="fout", 
                    help="plot results to PNG file",metavar="FILE.png")

persy.add_argument("-x", dest="x2", help="Transform bottom x-axis and set a label", 
                   nargs=2, metavar=("FAC","LABEL"))

persy.add_argument("-y", dest="y2", help="Right y-axis column 1-#", type=int, metavar="COL")

persy.add_argument("-m", "--movie", help="Create multiple pngs for movie creation",
                    nargs="?", const="frame")

persy.add_argument("-c","--constant",  dest="factor", help="multiply y-values by constant factor",
                     default=" 1.0", metavar="FACTOR")

persy.add_argument("-0", "--zero", help="substract y[0] from y-values",
                    action="store_true")

persy.add_argument("-s", "--smooth", 
                   help="smooth the data using Savistzky-Golay algorithm", action="store_true")

persy.add_argument("-z", help="Savistzky-Golay smoothing factors: WINDOW-size polynomial-ORDER", nargs=2, 
                   metavar=("WINDOW", "ORDER"), type=int, default=[51,3])

persy.add_argument("-d", "--header", type=int, nargs=1,
                   help="line number to get axes titles from", default=2)

persy.add_argument("-t","--title", nargs="*", help="title to print on top of the plot")
persy.add-argument("-v","--average",nargs=1,help="average data using a window", type=int)

args = persy.parse_args()
print args

x = []; y=[]
x2 = []; y2= []
fin = open(args.fin)
print("... opening "+args.fin)
nline = 0

for line in fin:
    nline +=1
    if line[0]=="#":
        if nline == args.header:
            xlbl = line.split()[args.x_col-1]
            ylbl = line.split()[args.y_col-1]
            if args.y2: y2lbl = line.split()[args.y2-1]
            print("... headers: "+xlbl+" "+ylbl)
        continue

    line = line.strip().split()
    x.append(float(line[args.x_col-1]))
    y.append(float(line[args.y_col-1]))
    if args.y2:
        y2.append(float(line[args.y2-1]))

x = np.array(x); y = np.array(y)   #change to np arrays
if args.average is not None:
    s = len(y)/args.average #number of windows
    r = len(y)%args.average #remainder
    y = y[:len(y)-r].reshape(-1,s)
    y = np.mean(y,1)
    x = x[:len(x)-r]
    x = [::s]
if args.zero:
    print("... zeroing y-data")
    y = y-y[0]           #zero results
if args.factor != " 1": 
    print("... scaling y-data by"+args.factor)
    y *= float(args.factor)                   #multiply by constant factor
yhat = savitzky_golay(np.array(y), args.z[0], args.z[1]) # window size 51, polynomial order 3
print("... %d datapoints"%(len(x)))


fig = plt.figure()
if args.title is not None:
    plt.title(args.title[0])
ax1 = fig.add_subplot(111)
ax1.plot(x,y,label=args.label, marker='.', linestyle='None')
if args.smooth:
    print("... smoothing")
    ax1.plot(x,yhat, color="red", label="smoothed Savitzky-Golay "+str(args.z))

#ax1.legend(loc="upper right", frameon=True)
ax1.set_xlabel(xlbl)
ax1.set_ylabel(ylbl)

marker_style = dict(color='cornflowerblue', linestyle=':', marker='o',
                    markersize=7, markerfacecoloralt='gray')

if args.x2 is not None:
    ax2 = ax1.twiny()
    ax2.set_xticks(ax1.get_xticks())
    ax2.set_xticklabels(np.array(ax1.get_xticks())*float(args.x2[0]))
    ax2.set_xlabel(args.x2[1])

if args.y2 is not None:
    ax3 = ax1.twinx()
    ax3.plot(x,y2,color="black",marker="s", linestyle="None")
    ax3.set_ylabel(y2lbl)

if args.movie is not None:
    load=["\\","|","/","-"]
    print("... making snapshots")
    a1 = ax1.plot(x[0],y[0], fillstyle='full', **marker_style)
    a2 = ax1.axhline(y=y[0],linestyle="--",color="grey")
    a3 = ax1.axvline(x=x[0], linestyle="--",color="grey")
    for i in range(len(x)):
        a1[0].set_data(x[i],y[i])
        a2.set_ydata(y[i])
        a3.set_xdata(x[i])
        plt.savefig(args.movie+"."+str(i)+".png")
        prcnt = float(i)/len(x)*100
        print(load[i%4]+ "     %d %%"%(prcnt))
        sys.stdout.write("\033[F")

elif args.fout is not None:
    print("... saving to file "+args.fout)
    plt.savefig(args.fout)
else:
    plt.show()
                    
