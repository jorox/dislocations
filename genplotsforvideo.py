import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

fin = open(fname)
#print("... opening "+sys.argv[1])
x = []
y = []
h1 = 1
h2 = 2
fac = 1
zero_results = True
fac = -1

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
        print(line)
        break

if zero_results:
    for i in range(n-1,-1,-1):
        y[i] =(y[i] - y[0]) #temporary remove it next time
    #print("... %d datapoints"%(n))
for i in range(start-1,end+1):

    plt.plot(x,y,label="data")
    plt.legend(loc="upper left", frameon=False)
    plt.title(r'Edge $b=[11\bar 20], (\bar 1100)$ Khater Bacon, Acta Materialia, 2010')
    plt.xlabel("strain %")
    plt.ylabel("Fx/lx/ly (MPa)")

    marker_style = dict(color='cornflowerblue', linestyle=':', marker='o',
                    markersize=15, markerfacecoloralt='gray')
    plt.plot(x[i],y[i], fillstyle='full', **marker_style)
    plt.savefig("frame."+str(i)+".png")
    plt.clf()
    print(i)




