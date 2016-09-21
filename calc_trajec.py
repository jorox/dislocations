#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from math import factorial
import numpy as np
import sys
import matplotlib.pyplot as plt


def process_single_file(fname,head_name):
    fin = open(fname)
    count = 0
    natoms = 0
    col    = -1
    av = 0.0

    for line in fin:
        count += 1
        if count == 4: 
            natoms = int(line)
            continue
        if count == 9: 
            header = line.strip().split()
            col = header.index(head_name)-2
            continue
        if count > 9:
            av += float(line.strip().split()[col])
    
    av /= natoms
    fin.close()
    return av

pers = ArgumentParser()

# Add options
pers.add_argument("fin", help="input file", metavar="FILE")
pers.add_argument("prop", help="data header name", default="x", metavar="PROPERTY")
pers.add_argument("units", help="data header units", default="A", metavar="UNITS")

pers.add_argument("-m", "--multi", help="multiple files", nargs=2, type=int,
                  metavar=("START", "END"))
pers.add_argument("-o","--output", help="output data", type=str, 
                  default="traj.dat")
pers.add_argument("--dt", help="timestep", nargs=2, metavar=("VALUE","UNITS"))

args = pers.parse_args()
print args
fout = open(args.output, "w")
if args.dt is None:
    out_header = "#step "
else:
    out_header = "#timestep"+"("+args.dt[1]+")"
out_header += args.prop+"("+args.units+")"
fout.write(out_header)

load = ["\\","|","/","-"]
print("... reading data")
if args.multi is not None:
    N = args.multi[1]-args.multi[0]+1
    y = []
    wild = args.fin.index("*")

    for i in range(args.multi[0],args.multi[1]+1):
        fname = args.fin[:wild]+str(i)+args.fin[wild+1:]
        y.append(process_single_file(fname,args.prop))
        if args.dt is None: fout.write("\n%1.0d %1.6f"%(i,y[-1]))
        else: fout.write("\n%1.4f %1.4f"%(i*float(args.dt[0]), y[-1]))
        
        prcnt = float(i)/N*100
        print(load[i%4]+ "     %d %%"%(prcnt))
        sys.stdout.write("\033[F")
    print("... done writing to"+args.output)

else:
    print("Not implemented yet...")
        
        
        
    

