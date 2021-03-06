import numpy as np
import matplotlib.pyplot as plt
import sys
from math import factorial


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


fin = open(sys.argv[1])
print("... opening "+sys.argv[1])
x = []
y = []
h1 = int(sys.argv[2])-1
h2 = int(sys.argv[3])-1
fac = 1
zero_results = False

if len(sys.argv)>4: 
	for i in range(4,len(sys.argv)):
		if sys.argv[i]=="-f":
			fac = float(sys.argv[i+1])
		if sys.argv[i]=="-0":
			zero_results = True

n = 0
for line in fin:
    try:
        if line[0] == "#" : continue
        line = line.strip().split()
        if len(line) == 0: continue
        
        if len(line)<h2: break
        n +=1
        x.append(float(line[h1]))
        y.append(fac*float(line[h2]))
    except ValueError:
        print line
        break

if zero_results:
	for i in range(n-1,-1,-1):
    		y[i] =(y[i] - y[0]) #temporary remove it next time
print("... %d datapoints"%(n))

yhat = savitzky_golay(np.array(y), 51, 3) # window size 51, polynomial order 3

plt.plot(x,y,label="Zr_3.eam",marker='.', linestyle='None')
plt.plot(x,yhat, color='red',label="smoothed savitzky-golay")
plt.legend(loc="upper right", frameon=False)
plt.title(r'Edge $b=[11\bar 20], (1\bar 100)$ Serra Bacon')
plt.xlabel("strain %")
plt.ylabel("Fx/lx/ly (MPa)")

marker_style = dict(color='cornflowerblue', linestyle=':', marker='o',
                    markersize=7, markerfacecoloralt='gray')
plt.plot(0.05,17, fillstyle='full', **marker_style)
plt.text(0.05,17,"(0.05,17)")
plt.text(0.15,-10,
"system relaxed using NVT only\n different dimensions as in article", fontsize=10, horizontalalignment='center',
         verticalalignment='top',
         multialignment='center')
plt.show()



