"""
Python code to modify node coordinates in <graph> type files
Code capabilities: 
  * shift along X direction (with periodic wrapping, entagelement should be corrected manually)
  * set several node coordinates to their average value (create a straight dislocation)
  * merge two segments
"""
import sys
import os
import glob
import numpy as np
from argparse import ArgumentParser
import inspect
from xml.etree import ElementTree as ET

dir2ind = {"x": 0, "y":1, "z":2, "X":0, "Y":1, "Z":2}

class node(object):
    def __init__(self,tag,seg,pinned,x):
        self.tag = tag #node tag
        self.seg = seg #tag of segment it belongs to
        self.pinned = pinned
        self.x = [] #coordinates
        for i in range(3): self.x.append(x[i])

    def __str__(self):
        return "%i pinned=%i [%1.4f, %1.4f, %1.4f]"%(self.tag, self.pinned, self.x[0],
                                                     self.x[1], self.x[2])

        
class nodelist:
    def __init__(self):
        self.nodes = []
        self.box_lim = [[1,1,1],[-1,-1,-1]]
    def __getitem__(self,i):
        return self.nodes[i]

    def tag2ind(self,tag_list):
        isalist = isinstance(tag_list,list)
        if False == isalist: tag_list= [tag_list]
        i_list = []
        for t in tag_list:
            ind = -1
            for i in range(len(self.nodes)):
                if self.nodes[i].tag == t: ind = i; break
            i_list.append(ind)
        if False == isalist: i_list = ilist[0]
        return i_list
    
    def add_node(self,p_node):
        if (isinstance(p_node,node)):
            eps = 1e-5
            #check box limits
            for ix in range(3):
                if p_node.x[ix] > self.box_lim[1][ix]:
                    self.box_lim[1][ix] = p_node.x[ix]
                    #check negative limits
                if p_node.x[ix] < self.box_lim[0][ix]:
                    self.box_lim[0][ix] = p_node.x[ix]
            self.nodes.append(p_node)
        else:
            print("Error: cannot add object, is not instance of class: node")


    def wrap_nodes(self,dirks):
        """
        wrap the nodes along specific directions to inside the box limits.
        This could be required because of a previous "displace" command  
        The function is intended for internal command only but can be used externally
        However, dirks in this case cannot be a string
        dirks <> 0 to 2 (inclusive)
        """
        for inode in self.nodes:
            for ix in dirks:
                if inode.x[ix] < self.box_lim[0][ix]: #node outside to the left
                    tmp = self.box_lim[0][ix]-inode.x[ix]
                    tmp = tmp % (self.box_lim[1][ix]-self.box_lim[0][ix])
                    inode.x[ix] = self.box_lim[1][ix] - tmp
                    print(">>> wrapping node %i along +%i to %1.4f"%(
                        inode.tag,ix,inode.x[ix]))
                          
                elif inode.x[ix] > self.box_lim[1][ix]: #node outside to the right
                    tmp = inode.x[ix]-self.box_lim[1][ix]
                    tmp = tmp % (self.box_lim[1][ix]-self.box_lim[0][ix])
                    inode.x[ix] = self.box_lim[0][ix] + tmp
                    print(">>> wrapping node %i along +%i to %1.4f"%(
                        inode.tag,ix,inode.x[ix]))
            
    def displace(self,dirk,val):
        """
        displace the nodes along a specific direction (with rewrapping)
        dirk = x,y,z
        """
        idirk = dir2ind[dirk]
        print(">>> WARNING: displacing nodes along %s, limits = %1.5f-->%1.5f by %1.4f"%(
            dirk, self.box_lim[0][idirk],self.box_lim[1][idirk], val))
        
        for i in range(len(self.nodes)): self.nodes[i].x[idirk] += val
        self.wrap_nodes([idirk]) # wrap nodes
                
    
#########################################################################################   

pers = ArgumentParser()

pers.add_argument("fin", help="input file name")
pers.add_argument("--dx", help="displace nodes along X", type=float)
pers.add_argument("-o", "--fout", help="ouput file name")
pers.add_argument("--ave", help="set DIM coordinate to average of nodes", nargs="+")

args = pers.parse_args()


fin = ET.parse(args.fin)
print("... opening XML file "+args.fin)
nodes = nodelist()
num_nodes = 0

#### Read nodes from XML and create a node list
for child in fin.iter():
    if child.tag == "node" and len(child.attrib) > 3: #child represents node
        pinned = True
        x = []
        if child.attrib["pinned"] == 0: pinned = False
        tag = int(child.attrib["tag"])
        for s in ["x","y","z"]: x.append(float(child.attrib[s]))
        tmp = node(tag,pinned, x) #create temp. node-object
        nodes.add_node(tmp) #add to node list
        num_nodes += 1
print("... done reading %i nodes"%(num_nodes))
print("box limits: \n    %1.4f --> %1.4f\n    %1.4f --> %1.4f\n    %1.4f --> %1.4f"%(
    nodes.box_lim[0][0],nodes.box_lim[1][0],
    nodes.box_lim[0][1],nodes.box_lim[1][1],
    nodes.box_lim[0][2],nodes.box_lim[1][2]))
for i in range(num_nodes):
    print str(nodes[i])


##### Modify the node list with optional args
## Displace
if args.dx is not None:
    nodes.displace("x", args.dx) #displace the nodes
    print("... displaced nodes along x by %1.5f"%(args.dx))
    for i in range(num_nodes): print("     "+str(nodes[i]))

## Average
dim_2_i = {"x":0, "X":0, "y":1, "Y":1, "z":2, "Z":2}
if args.ave is not None:
    if len(args.ave) < 3:                                                       #! not enough nodes
        print(">>> ERROR: insufficient nodes to average")
        sys.exit()
    if args.ave[0] not in ["x","y","z","X","Y","Z"]:                            #! Invalid dimension
        print(">>> ERROR: invalid dimension specified for averaging command")
        sys.exit()
    av = 0                                                                      #  Initialize the average
    idim = dim_2_i[args.ave[0]]                                                 #  Find the index of the dim
    try:
        taglist = [int(tmp) for tmp in args.ave[1:]]                            #  change tags to ints
    except ValueError:
        print(">>> ERROR: Invalid node tag for averaging, must be an integer")  #! tag is not an int
        sys.exit()
    inode = nodes.tag2ind(taglist)                                              #  list indices of the tags
    if -1 in inode:                                                             #! tag does not exsist
        print(">>> ERROR: Node tag does not exist")
        print(args.ave); print(inode); sys.exit()     
    for i in inode: av += nodes[i].x[idim]                                      #  calc the average
    av /= len(inode)
    for i in inode: nodes[i].x[idim] = av                                       #  Modify the nodes
    print("... averaged nodes along dim-%i to %1.5f"%(idim,av))
    for i in inode: print("     "+str(nodes[i]))

    
##### Modify the XML
ichild = -1
for child in fin.iter():
    if child.tag == "node" and len(child.attrib) > 3: #find all "valid" node children i.e. with coords
        ichild += 1
        if int(child.attrib["tag"]) == nodes[ichild].tag: #tags match
            child.attrib["x"] = "%1.16f"%(nodes[ichild].x[0]) #modify the x-coordinate
            child.attrib["y"] = "%1.16f"%(nodes[ichild].x[1])
            child.attrib["z"] = "%1.16f"%(nodes[ichild].x[2])
        else:
            print("ERROR: inconsistent tags - XML or node list order modified")
            print("Original file not modified")
            exit
print("... original data has been modified but not written to disk")


    
if args.fout is not None:
    fin.write(args.fout)
    print("... done writing XML to file %s"%(args.fout))
else:
    print("... no XML file written")

