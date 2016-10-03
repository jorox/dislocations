#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
import numpy as np
import sys
import math
# Import OVITO modules.
from ovito.io import *
from ovito.modifiers import *
from ovito.data import *
from ovito.vis import *
from PyQt5 import QtCore

# This python script is meant to be run with Ovitos to create 4 movies 
# representing the 4 main windows of the ovito screen. The input data
# is the snapshots deform.*.dat files and perhaps the png snaphots of
# a stress-strain curve


pers = ArgumentParser()

pers.add_argument("fin", help="Input file", metavar="File-IN")
pers.add_argument("-r","--range", help="Range", nargs=3, type=int, metavar=("START","END","INCR"))
pers.add_argument("-o","--output", help="Output File Name base", metavar="FOUT.DIM.avi")
pers.add_argument("-c","--cna", help="CNA cutoff", type=float, metavar="CUTOFF")
pers.add_argument("-t","--time", help="Total length of video in seconds",type=int,metavar="TIME")
args = pers.parse_args()

# find wild-card in filename
iwild = args.fin.find("*")
if iwild == -1:
    print("Fatal Error: filename does not contain wild-card '*'")
    exit
print(args)
# Load file source 
node = import_file(args.fin)
print("... "+str(node.source.source_path)+ " loaded")
print("... "+str(node.source.num_frames)+" frames found")
# Handle warnings
if args.cna is None:                     #cna cutoff unspecified
    print("-->warning CNA cutoff undefined, using default 3.86")
    args.cna = 3.86
if args.range is None:
    args.range = [0,node.source.num_frames,1]
if args.range and args.range[1]>node.source.num_frames: #range is too high
    print("-->warning range specified too large")
    args.range[1] = node.source.num_frames
if args.output is None:
    args.output = "mov"
print("... range = %d:%d:%d"%(args.range[0],args.range[1],args.range[2]))

# MODIFIERS
#----------------------------------------------------------
print("... CNA cutoff = "+str(args.cna))
cna = CommonNeighborAnalysisModifier(mode = CommonNeighborAnalysisModifier.Mode.FixedCutoff, cutoff=args.cna)
node.modifiers.append(cna)
select_normal = SelectExpressionModifier(expression = 'StructureType>0')
node.modifiers.append(select_normal)
node.modifiers.append(DeleteSelectedParticlesModifier())
print("... modifiers added")
    
# COMPUTE
#----------------------------------------------------------
load = ["\\","|","/","-"]
N = (args.range[1]-args.range[0])/args.range[2]  #number of video frames
print("...%1.0d video frames\n    processing:"%(N))
for i in range(args.range[0],args.range[1],args.range[2]):
    node.compute(i)
    prcnt = float(i)/N*100
    print(load[i%4]+ "     %d %%"%(prcnt))
    sys.stdout.write("\033[F")


# Render
#-----------------------------------------------------------
node.add_to_scene()
cell = node.source.cell
cell.display.enabled = True
cell.display.rendering_color = (1,1,1)



vp = Viewport()
vp.type = Viewport.Type.PERSPECTIVE
vp.zoom_all()
vp.camera_pos = (-100, -150, 150)
vp.camera_dir = (2, 3, -3)
vp.fov = math.radians(20.0)

#Create an overlay.
tripod = CoordinateTripodOverlay()
tripod.size = 0.07
tripod.alignment = QtCore.Qt.AlignRight ^ QtCore.Qt.AlignBottom

# Attach overlay to the active viewport.
#viewport = ovito.dataset.viewports.active_vp
#viewport.overlays.append(tripod)
vp.overlays.append(tripod)

if args.time is not None:
    fps = N/args.time
else:
    fps = 10
ovito.dataset.anim.frames_per_second = fps

rs = RenderSettings(
    filename = args.output+".avi",
    size = (1024, 768),
    background_color=(0,0,0),
    custom_range = (args.range[0],args.range[1])
)
#rs.renderer = TachyonRenderer()
rs.range = RenderSettings.Range.CUSTOM_INTERVAL
print(rs.range)
rs.renderer.antialiasing = True
print("... rendering movie")
vp.render(rs)
