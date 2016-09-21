#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from math import factorial
import numpy as np
import sys
import matplotlib.pyplot as plt


def process_single_file(fname,dsp_1,dsp_3,n1,n2,nd):

# calculate the displacement along dir_1 between atoms across dsp_3
#    dsp_1: slip direction/burgers vector
#    dsp_3: normal to slip direction
#    n1   : number of atoms along dsp_1
#    n2   : number of atoms along dsp_2
#    nd   : number of atoms deleted from lower plane
    fin = open(fname)
    count = 0
    natoms = 0
    midplane = 0
    col_1    = -1
    col_2    = -1
    atoms = []
# read data from file and place in array    
    for line in fin:
        count += 1
        if count == 4: 
            natoms = int(line)
            continue
        if count == 9: 
            header = line.strip().split()
            col_1 = header.index(dsp_1)-2
            s = ["x","y","z"]; s.pop(s.index(dsp_1)); s.pop(s.index(dsp_3))
            col_2 = header.index(s[0])-2
            col_3 = header.index(dsp_3)-2
            continue
        if count > 9:
            line = line.strip().split()
            atoms.append( [int(line[0]), 
                           float(line[col_1]), float(line[col_2]), float(line[col_3])] )
            midplane += float(line[col_3])

# sort atoms according to z
    if n1*n2*2-nd*n2 == natoms: print("   ... perfect count") # ensure count is correct
    midplane /= natoms # z-midplane coordinate
    N = n1*n2          # number of atoms in upper
    nr = n2*nd         # number of atoms missing from lower because of the dislocation
    
    atoms = sorted(atoms, key=lambda x: x[3]) # sort atoms in each block according to z
    mdi_3 = 0                                 # find index of element that marks the midplane

    for i in range(len(atoms)):
        if atoms[i][3]>=midplane:
            mdi_3 = i
            print("#######DEBUG#####\n first atom upper plane = "+ str(i) + " " +
                  str(atoms[i])+"\n######DEBUG####\n");break
    if mdi_3<N-nr: 
        print("warning in file "+fname+" "+": midplane ill-defined")
        mdi_3 += N-nr-mdi_3
    # take 2nxny-ndny elements - this removes any outliers --> 0:N-nd*n2=lower; N-nd*n2:2N-nd=upper
    atoms = atoms[nr-N+mdi_3:mdi_3+N] #total number = 2N-nr

# sort lower/upper by x
    atoms[:N-nr] = sorted(atoms[:N-nr], key=lambda x:x[1])
    atoms[N-nr:] = sorted(atoms[N-nr:], key=lambda x:x[1])

# sort each x-plane upper and lower by y
    for i in range(2*n1-nd):
        atoms[i*n2:(i+1)*n2] = sorted(atoms[i*n2:(i+1)*n2],key=lambda x:x[2])
    atoms = np.array(atoms)
# atoms now sorted by y, then x, then z
    
    # DEBUG
    #print("\n====DEBUG-sorted atoms====")
    #print(atoms)
    tmpout = open("dirty.dat",'w')
    for jj in range(len(atoms)): tmpout.write("%i %i %1.8f %1.8f %1.8f\n"%(
            jj,atoms[jj,0],atoms[jj,1], atoms[jj,2], atoms[jj,3]))
    #print("====DEBUG====")

# calculate displacements
    delta_u = []                           # displacements
    x = []
    nplower = n1-nd                        # number of planes in lower 
    L = atoms[(nplower-1)*n2,1]-atoms[0,1] # length of lower plane
    L = L + L/nplower                      # the periodic boundary contains a ghost atom
    print("lower length = "+str(L))
    for i in range(n1):
        s = i*n2            # index first atom along D1 - upper
        t = s+n2            # number of atoms along D2 - upper
        sp = s%(nplower*n2) # index atoms in lower - wrapped 
        tp = sp+n2          # number of atoms in lower along D2
        
        #print("###DEBUG -sp,tp####\n"+str(sp)+", "+str(tp)+"\n#####DEBUG###\n")
        if i==0:
            print("#### DEBUG - CASE 0 ####\n"+str(atoms[sp:tp,1])+"\n"+str(atoms[N-nr+s:N-nr+t,1])+
                  "\n###### DEBUG #####\n")
        tmp = atoms[sp:tp,1]-atoms[N-nr+s:N-nr+t,1] #lower-upper
        delta_u.append(np.abs(tmp)%L*np.sign(tmp))
        x.append(atoms[N-nr+s:N-nr+t,1])

    delta_u=np.array(delta_u)
    x = np.array(x)
    #print("\n====DEBUG====")
    #print(delta_u)
    #print("====DEBUG====")

    delta_u = np.average(delta_u,1)
    x = np.average(x,1)
    #print("\n====DEBUG====")
    #print(delta_u)
    #print("====DEBUG====")


    fin.close()
    return zip(x,delta_u)

def derv(x,a):
    x = np.array(x)
    a = np.array(a)
    return np.diff(a)/np.diff(x)
####################################################################################

pers = ArgumentParser()

# Add options
pers.add_argument("fin", help="input file", metavar="FILE")
pers.add_argument("D1", help="lateral direction", default="x", metavar="D1")
pers.add_argument("np", help="number of atomic planes along and perp to D1", 
                  metavar=("Nprll","Nperp"), nargs=2, type=int)
pers.add_argument("D3", help="vertical direction", default="z", metavar="D3")
pers.add_argument("nd", help="difference in number of atoms", metavar="N_D",
                  type=int)
pers.add_argument("-u","--units", help="data header units", default="A", metavar="UNITS")

pers.add_argument("-m", "--multi", help="multiple files", nargs=3, type=int,
                  metavar=("START", "END", "INCR"))
pers.add_argument("-o","--output", help="output data", type=str, 
                  default="delta_u.dat")

args = pers.parse_args()
print args

if args.multi is not None:
    i_file = range(args.multi[0],args.multi[1]+args.multi[2],args.multi[2])
else:
    i_file = [1]

fout = open(args.output, "w")
out_header="#"+args.D1+"(atomic plane)" 
for i in i_file:
    out_header+= " Du_"+str(i)+"("+args.units+")"
fout.write(out_header)

load = ["\\","|","/","-"]
Du = []
DDu = []
N = len(i_file)
cnt=0
if args.multi is not None:
    print("... processing multi-files")
    wild = args.fin.index("*")
    for i in i_file:
        fname = args.fin[:wild]+str(i)+args.fin[wild+1:]
        Du.append(process_single_file(fname,args.D1,args.D3,args.np[0],args.np[1],args.nd))
        DDu.append(derv(zip(*Du[-1])[0],zip(*Du[-1])[1]))
        cnt +=1
        prcnt = float(cnt)/N*100
        print(load[cnt%4]+ "     %d %%"%(prcnt))
        sys.stdout.write("\033[F")
else:
    print("... processing single file")
    Du.append(process_single_file(args.fin,args.D1,args.D3,args.np[0],args.np[1],args.nd))
    DDu.append(derv(zip(*Du[-1])[0],zip(*Du[-1])[1]))

for j in range(args.np[0]):
    fout.write("\n")
    for i in range(len(i_file)):
        fout.write("%1.6f %1.6f %1.6f "%(Du[i][j][0],Du[i][j][1],DDu[i][j%(args.np[0]-1)]))

fout.close()
print("... Done writing to data to "+args.output)
