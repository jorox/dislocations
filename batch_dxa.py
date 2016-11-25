#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
import numpy as np
import sys
import math
import glob

# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.vis import *
import os

# This python script is meant to be run with Ovitos to create 4 movies 
# representing the 4 main windows of the ovito screen. The input data
# is the snapshots deform.*.dat files and perhaps the png snaphots of
# a stress-strain curve

pers = ArgumentParser()

pers.add_argument("fin", help="Input file", metavar="File-IN motif")
pers.add_argument("-o", "--output", metavar="File-OUT motif", default="batchdxa")
pers.add_argument("-r","--range", help="Range", nargs=3, type=int, metavar=("START","END","INCR"))
pers.add_argument("-p","--particles", help="output particle position instead of dislocations",
                  action="store_true")
args = pers.parse_args()

#Load file
#-------------------------------------------------------------------
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

node = import_file(args.fin)
print("... "+str(node.source.source_path)+ " loaded")
print("... "+str(node.source.num_frames)+" snapshots found")


#Modifiers
#-------------------------------------------------------------------
if args.particles:
    modf = CommonNeighborAnalysisModifier(mode   = CommonNeighborAnalysisModifier.Mode.FixedCutoff,
                                         cutoff = 3.86)
    node.modifiers.append(modf)
    select_normal = SelectExpressionModifier(
    expression = 'StructureType>0 || abs(Position.Z)>CellSize.Z/4')
    node.modifiers.append(select_normal)
    node.modifiers.append(DeleteSelectedParticlesModifier())
    print("... Common-Neighbor Analysis cutoff = 3.86")
else:
    modf = DislocationAnalysisModifier()
    modf.input_crystal_structure = DislocationAnalysisModifier.Lattice.HCP
    node.modifiers.append(modf)
    print("... Dislocation analysis for %s crystal"%(modf.input_crystal_structure))
    
#print(node.modifiers[0])

print("... modifiers added")


# COMPUTE
#----------------------------------------------------------

print("... writing to %s"%(os.getcwd()+"/"+args.output+".*"))
load = ["\\","|","/","-"]
N_snaps = (args.range[1]-args.range[0])/args.range[2]  #number of selected snapshots
print("... %1.0d video frames\n    processing modifiers:"%(N_snaps))
tmp = -1
for i in range(args.range[0],args.range[1],args.range[2]):
    tmp += 1
    ovito.dataset.anim.current_frame = i
    try:
        node.compute()
    except RuntimeError as err:
        print("Runtime Error: Cannot compute file %s"%(flist[i]))
        continue
    if not args.particles:
        total_line_length = node.output.attributes['DislocationAnalysis.total_line_length']
    else:
        total_line_length = 0
        
    prcnt = float(tmp)/N_snaps*100
    fname = os.getcwd()+"/"+args.output+"."+str(keys[i])
    print(load[tmp%4]+ "     %d %% Dislocation density: %g m^-2 --> %s"%(prcnt,
                                                                         total_line_length/node.output.cell.volume*1e20,
                                                                         fname))
    #sys.stdout.write("\033[F")

    if not args.particles:
        export_file(node, fname, "ca")
    else:
        export_file(node, fname, "lammps_dump",
                    columns = ["Particle Identifier", "Position.X", "Position.Y", "Position.Z"])
        
print("done")



