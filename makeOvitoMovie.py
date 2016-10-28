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
from PyQt5.QtGui import QPainter, QColor, QPen
import os

# This python script is meant to be run with Ovitos to create 4 movies 
# representing the 4 main windows of the ovito screen. The input data
# is the snapshots deform.*.dat files and perhaps the png snaphots of
# a stress-strain curve

pers = ArgumentParser()

pers.add_argument("fin", help="Input file", metavar="File-IN")
pers.add_argument("-r","--range", help="Range", nargs=3, type=int, metavar=("START","END","INCR"))
pers.add_argument("-o","--output",help="Output File Name base", metavar="FOUT.DIM.avi")
pers.add_argument("-c","--cna",  help="CNA cutoff", type=float, metavar="CUTOFF")
pers.add_argument("-t","--time", help="Total length of video in seconds",type=float,metavar="TIME")
pers.add_argument("-v", "--vpcount",
                  help="Number of Viewports to output (0=pers,1=pers+front,2=pers+front+top",
                  type=int, default=1)
pers.add_argument("-d", "--dxa", help="Perform dislocation analyis", action="store_true")
pers.add_argument("--tdump", help="simulation dump interval in ps",nargs=2, default=[1," "])
args = pers.parse_args()
args.tdump[0] = float(args.tdump[0])


def tsim_render_overlay(painter, **kargs):
    txt_color = QColor(255,255,255)
    painter.setPen(txt_color)
    tsim = args.tdump[0]*ovito.dataset.anim.current_frame
    painter.drawText(10,20,"%1.4f %s"%(tsim,args.tdump[1]))


# find wild-card in filename
iwild = args.fin.find("*")
if iwild == -1:
    print("Warning: ONLY ONE FILE SPECIFIED")
    args.range = [0,1,1]
    
print(args)
# Load file source 
node = import_file(args.fin)
print("... "+str(node.source.source_path)+ " loaded")
print("... "+str(node.source.num_frames)+" snapshots found")
# Handle warnings
if not args.dxa and args.cna is None:                     #cna cutoff unspecified
    print("--> warning CNA cutoff undefined, using default 3.86")
    args.cna = 3.86
if args.range is None:
    args.range = [0,node.source.num_frames,1]
if args.range[1]>node.source.num_frames: #range is too high
    print("-->warning range specified too large, respecifying...")
    args.range[1] = node.source.num_frames
if args.output is None:
    args.output = "mov"
print("... range = %d:%d:%d"%(args.range[0],args.range[1],args.range[2]))

# MODIFIERS
#----------------------------------------------------------
if not args.dxa:
    print("... CNA cutoff = "+str(args.cna))
    cna = CommonNeighborAnalysisModifier(mode   = CommonNeighborAnalysisModifier.Mode.FixedCutoff,
                                         cutoff = args.cna)
    node.modifiers.append(cna)
    select_normal = SelectExpressionModifier(
        expression = 'StructureType>0 || abs(Position.Z)>CellSize.Z/4')
    node.modifiers.append(select_normal)
    node.modifiers.append(DeleteSelectedParticlesModifier())
else:
    dxamod = DislocationAnalysisModifier()
    dxamod.input_crystal_type = DislocationAnalysisModifier.Lattice.HCP
    node.modifiers.append(dxamod)
    print("... Dislocation analysis for %s crystal"%(dxamod.input_crystal_type))
print("... modifiers added")
    
# COMPUTE
#----------------------------------------------------------
load = ["\\","|","/","-"]
N_video_frames = (args.range[1]-args.range[0])/args.range[2]  #number of video frames
print("... %1.0d video frames\n    processing modifiers:"%(N_video_frames))
for i in range(args.range[0],args.range[1],args.range[2]):
    node.compute(i)
    prcnt = float(i)/N_video_frames*100
    print(load[i%4]+ "     %d %%"%(prcnt))
    sys.stdout.write("\033[F")


# Render
#-----------------------------------------------------------
if args.dxa:
    print("... Found %i dislocation segments"%len(node.output.dislocations.segments))
    for segment in node.output.dislocations.segments:
        print("Segment %i: length=%f, Burgers vector=%s" %
              (segment.id, segment.length, segment.true_burgers_vector))
        print(segment.points)
    node.output.particle_properties.position.display = None
    node.output.dislocations.display.enabled = True
    
node.add_to_scene()
cell = node.source.cell
cell.display.enabled = True
cell.display.rendering_color = (1,1,1)

for ivp in range(args.vpcount):
    vp = Viewport()
    if ivp == 0:
        vp.type = Viewport.Type.PERSPECTIVE
        vp.camera_pos = (-70, -500, 100)
        vp.camera_dir = (0.2, 1, -0.25)
    elif ivp ==1:
        vp.type = Viewport.Type.FRONT
        vp.zoom_all()
    elif ivp ==2:
        vp.type = Viewport.Type.TOP
        vp.zoom_all()
    else:
        print("Error: Invalid number of Viewports")

    tripod = CoordinateTripodOverlay()
    tripod.size = 0.1
    tripod.alignment = QtCore.Qt.AlignRight ^ QtCore.Qt.AlignBottom

    tsim_overlay = PythonViewportOverlay(function = tsim_render_overlay)
    # Create an overlay.
    #tstep_text_overlay = TextLabelOverlay(
    #    text = str('timestep: [Timestep]'), 
    #alignment = QtCore.Qt.AlignLeft ^ QtCore.Qt.AlignBottom,
    #offset_y = 0.05,
    #font_size = 0.01,
    #text_color = (1,1,1))
    #tstep_text_overlay.source_node = node

# Attach overlay to the active viewport.
#viewport = ovito.dataset.viewports.active_vp
#viewport.overlays.append(overlay)
#vp.zoom_all()
#vp.fov = math.radians(20.0)
#Create an overlay.

# Attach overlay to the active viewport.
#viewport = ovito.dataset.viewports.active_vp
#viewport.overlays.append(tripod)
    vp.overlays.append(tripod)
    vp.overlays.append(tsim_overlay)
    #vp.overlays.append(tstep_text_overlay)
    if args.time is not None:
        fps = int(N_video_frames)/args.time
        if fps < 1:
            fps = 1
            print("--> warning cannot yet use fps lower than 1, resetting fps to 1.0")
            print("    new video length = %1.2f seconds"%(fps*N_video_frames))
    else:
        fps = 10
    ovito.dataset.anim.frames_per_second = fps
    print("... using %1.4d fps"%(fps))
    rs = RenderSettings(
        filename = args.output+"."+str(vp.type)+".avi",
        size = (1024, 768),
        background_color=(0,0,0),
        custom_range = (args.range[0],args.range[1])
        )
    #rs.renderer = TachyonRenderer()
    if args.range[1]>1:
        rs.range = RenderSettings.Range.CUSTOM_INTERVAL  #video
    else:
        rs.filename=args.output+"."+str(vp.type)+".png"                  #image

    rs.renderer.antialiasing = True
    print("*** rendering movie: %s, %d fps, %s ***"%(str(vp.type),fps, str(rs.range)))
    vp.render(rs)
    
    print("... done writing to "+os.getcwd()+"/"+rs.filename)
    
