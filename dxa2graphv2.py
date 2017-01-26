"""
This program transforms a DXA file to a graph file for DD simulation
Assumptions:
 ** Working for HCP structure only
     Other structures are not considered
 ** Box is assumed to be orthogonal even for screw dislocation simulations
     This is to help in wrapping vertices across the periodic boundaries
 ** Main directions are set
     X=[11-20], Z=[-1100], Y=[0001]
"""
import sys
import os
import glob
import numpy as np
from argparse import ArgumentParser
import decimal
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.dom import minidom

class mb:
    """
    The class is meant to handle Miller-Bravais style vectors
    Its methods are still unimplemented
    To complete the code we need to finish all the arithmetic functions: add, divide, multiply
    and the inner, outer product functions
    For simplicity and efficiency these functions are not yet implemented
    """
    def __init__(self,data):
        for idat in range(4):
            if np.abs(data[idat]-int(data[idat]))<.1: data[idat] = int(data[idat])
            else: data[idat] = round(data[idat],1)
        self.coords = list(data)
        
    def __str__(self):
        strtmp = "[%i%i%i%i]"%(self.coords[0],self.coords[1],self.coords[2],self.coords[3])
        return strtmp

    def norm(self):
        nrm = 0.0
        for i in range(4): nrm += self.coords[i]**2
        nrm -= (self.coords[0]*self.coords[1]+
                self.coords[0]*self.coords[2]+
                self.coords[1]*self.coords[2])
        return nrm

    def __getitem__(self,key):
        return self.coords[key]
    def __setitem__(self,key,val):
        self.coords[key] = val

    def verify(self):
        """ return True if the instance has coordinates that are quotients of 6 i.e. 1/6 and
        i=-(h+k)
        """
        for i in self.coords:
            if np.abs(6*i-int(6*i))>0.1: return False
        if np.abs(self.coords[2]+self.coords[0]+self.coords[1]) > 0.1: return False
        return True

class node:
    """
    Representation of a segment node
    """
    def __init__(self,data):
        # data = [tag,x,y,z]
        self.coords = list(data[1:4])
        self.tag = data[0]
    def __getitem__(self,key):
        return self.coords[key]
    def __setitem__(self,key,val):
        self.coords[key] = val
    def __str__(self):
        return "tag:%i coords:%1.3f, %1.3f, %1.3f"%(self.tag,self.coords[0],self.coords[1],self.coords[2])
    def __add__(self,p2):
        return [x+y for x,y in zip(self.coords,p2.coords)]
    def __sub__(self,p2):
        return [x-y for x,y in zip(self.coords,p2.coords)]
        
class cluster:
    def __init__(self,index,stype,A):
        self.index = int(index)
        self.stype = stype
        self.A = np.array(A) #matrix is defined row by row
    def __str__(self):
        tmps = "cluster %i: %s type\n%1.3f %1.3f %1.3f\n%1.3f %1.3f %1.3f\n%1.3f %1.3f %1.3f"%(
            self.index, self.stype, self.A[0,0],self.A[0,1],self.A[0,2],
            self.A[1,0],self.A[1,1],self.A[1,2],
            self.A[2,0],self.A[2,1],self.A[2,2])
        return tmps

class segment:
    """
    class representation of a segment.
    Data members:
    * index = index from DXA file
    * localb = true Burgers vector calculated from input data
    * xi = normal to glide plane
    * vertices = a list containing links to the nodes that belong to the segment
    """
    def __init__(self,index,b_coords,xi_coords,vertices):
        self.index = index
        self.b = mb(b_coords)
        self.xi = mb(xi_coords)
        self.vertices = list(vertices)
        
    def __str__(self):
        tmps = "segment %i:\n  6*b= %1.3f %1.3f %1.3f %1.3f\n  6*xi= %1.3f %1.3f %1.3f %1.3f\n  vertices { "%(
            self.index,
            self.b[0], self.b[1], self.b[2], self.b[3],
            self.xi[0], self.xi[1],self.xi[2], self.xi[3],
            )
        for i in range(len(self.vertices)):
            tmps += "%i "%(self.vertices[i].tag)
        tmps += "}"
        return tmps
    
    def __len__(self):
        return len(self.vertices)


class configuration:
    def __init__(self,list_segments,list_clusters,topo):
        """
        function constructs a configuration based on input date from DXA file
        * LIST_VERTICES:  vertices data [[x0,y0,z0], ...]
        * LIST_SEGMENTS:  segment data  [id, [bx,by,bz], cluster, N_vertices, [x0,y0,z0], [x1,y1,z1], ... ]
        * LIST_CLUSTERS:  cluster data  [id, structure, [A00,A01,A02], [A10,A11,A12], [A20,A21,A22]]
        * LIST_JUNCTIONS: 
        """
        self.nodes = []
        self.repeated = []
        N_vert = 0
        self.count = 0  #number of distinct nodes
        #######################################################################
        # self.nodes contains references to node objects
        # list_vertices is a list containing all vertex listings in the DXA file
        # ordered by their appereance in the DISLOCATION section.
        # Repetition of a node is handeled by having the same object address
        # self.nodes will thus have the same length as list_vertices
        # self.nodes[i] cannot thus be identified to belong to which segment
        # exactly
        #######################################################################
        for d in list_segments:
            for vertex in d[4:]:
                data = [self.count,None,None,None]
                for ix in range(3): data[ix+1] = vertex[ix]
                        
                # check for repetition, a tolerance is added because some repeated nodes might get wrapped
                #  hence the coordinates are not exactly the same
                found = False
                if len(self.nodes)>0:
                    for i_node in self.nodes:
                        if np.sqrt((i_node[0]-vertex[0])**2+
                                   (i_node[1]-vertex[1])**2+(i_node[2]-vertex[2])**2 ) < 0.1: #same point
                            self.nodes.append(i_node) #add the same address to the list
                            found = True
                            self.repeated.append(1)
                            break
                    
                if found == False:
                    self.repeated.append(0)
                    self.nodes.append(node(data)) #new node
                    self.count+=1

        self.segments = []
        ############################################################################
        # Determine the true Burgers vector for each dislocation
        #  Calculate the Spatial Burgers Vector by multiplying the
        #  local Burgers vector by the matrix of the cluster
        #  divide the spatial-b by the lattice spacings along the 3 main directions
        #
        # *************************************************************************
        #
        # Determine the glide plane of the segment by taking the cross product
        #  between the dislocation line sense(first to last vertex) and the Burgers
        #  vector. This can be done in mb or real space and then transformed
        #  The sense vector is determined from the unwrapped coordinates of the
        #  segments because it will help for the case when the segment passes
        #  through a periodic boundary
        ############################################################################
        
        X = [1./3,1./3,-2./3,0.]; Z = [-1.,1.,0.,0.]; Y = [0.,0.,0.,1.]
        box_spc = [3.232,5.165,3.232*np.sqrt(3)]
        print("\n... principla directions and spacings along")
        print("     X = "+str(X)+"  "+str(box_spc[0])+" A")
        print("     Y = "+str(Y)+"  "+str(box_spc[1])+" A")
        print("     Z = "+str(Z)+"  "+str(box_spc[2])+" A")
        X = np.array(X); Y = np.array(Y); Z = np.array(Z)
        ivert = 0                                                 #index to first vertex of current segment
        
        for d in list_segments:
            tmp = np.array(d[1])                                 # spatial Burgers vector
            N_verts = d[3]                                       # number of vertices
            idcluster = d[2]-1                                   # cluster id

            # Burgers vector
            rotmat = np.array(list_clusters[idcluster][2:5])
            spcb = rotmat.dot(tmp)                                  # left-mult local Burgers by rotation matrix
            for i in range(3):  spcb[i] = spcb[i]/box_spc[i]       # scale vector (no units)
            trueb = 6.0*(X*spcb[0] + Y*spcb[1] + Z*spcb[2])         # transform to mb notation

            # Determine glide plane normal
            sens = np.diff(d[4:4+N_verts],axis=0)                      # calculate partial sense vectors
            sens = np.sum(sens,axis=0)                                 # sum
            sens = sens/N_verts                                        # average
            b2pr = (np.sqrt(2)*sens[1])/np.sqrt(sens[0]**2+sens[2]**2) # basal/prismatic ratio
            if b2pr < 1.0: sens[1] = 0
            normsens = np.linalg.norm(sens)
            for i in range(3): sens[i] = sens[i]/box_spc[i]    # scale vector (no units)
            #for i in range(3): sens[i] /= normsens            # normalize vector
            print("########## DEBUG ###########")
            print(sens)
            print(spcb)
            print(b2pr)
            truexi = np.cross(spcb,sens)
            print(truexi)
            truexi = 3.0*(X*truexi[0] + Y*truexi[1] + Z*truexi[2])      # transform xi to mb
            print(truexi)
            print("########## DEBUG ###########")

            # Vertex list
            truevertex = list(self.nodes[ivert:ivert+N_verts])
            self.segments.append(segment(d[0],trueb,truexi,truevertex))
            ivert = ivert+N_verts

        self.entangled_nodes = []
        ############################################################################
        # split segments running across a periodic boundary into two. There might be
        # some segments which pass through a periodic boundary. After wrapping the
        # nodes it is important to split these segments for NUMODIS. The segments
        # can be identified by a jump in the wrapped node coordinates of by exceeding
        # the box limits in the unwrapped coordinates
        #
        # ! For the moment it is assumed that segments only cross one dimension. If
        # the segment crosses 2 or more it is uncertain whether the code will work
        ############################################################################
        for d in self.segments:
            for ix in range(3):
                for iv in range(len(d.vertices)):
                    pnt = d.vertices[iv][ix]
                    pnt_out = ( pnt <= topo[1][ix]) and ( pnt >= topo[0][ix] )
                    if iv == 0 : #first point
                        first_out = pnt_out 
                        continue
                    else: #not first point
                        if pnt_out != first_out: #changed location (can be inside or outside)
                            self.segments.append(segment(d.index,d.b,d.xi,d.vertices[iv:]))
                            self.entangled_nodes.append([d.vertices[iv-1],d.vertices[iv]])
                            d.vertices = d.vertices[:iv]
                            print("   ****splitting segment %i at node %i"%(d.index,iv))
                            break

        #############################################################################
        # wrap nodes across periodic boundaries
        #############################################################################
        for vertex in self.nodes:
            if vertex[ix]<topo[0][ix]: vertex[ix] = vertex[ix]-topo[0][ix]+topo[1][ix] # x<origin
            elif vertex[ix]>topo[1][ix]: vertex[ix] = vertex[ix]-topo[1][ix]+topo[0][ix] # x>limit

        print("... Done building configuration")
        
        
    def build_xml_tree(self):
        
        xml_root = Element("root")
        xml_nodes = SubElement(xml_root,"nodes")
        for inode in range(len(self.nodes)):
            thisnode = self.nodes[inode]
            repeated = False
            if self.repeated[inode] == 1: continue

            tmp = SubElement(xml_nodes,"node",{
                      "pinned":"0",
                      "tag":str(thisnode.tag),
                      "x":"%1.12f"%(thisnode[0]),
                      "y":"%1.12f"%(thisnode[1]),
                      "z":"%1.12f"%(thisnode[2])})

        xml_lines = SubElement(xml_root,"lines")
        for seg in self.segments:
            xml_line = SubElement(xml_lines,"line",{
                           "bh":str(seg.b[0]),
                           "bk":str(seg.b[1]),
                           "bi":str(seg.b[2]),
                           "bl":str(seg.b[3]),
                           "tag":str(seg.index)})
            xml_constraint = SubElement(xml_line,"constraints")
            xml_plane = SubElement(xml_constraint,"plane",{
                           "h":str(seg.xi[0]),
                           "k":str(seg.xi[1]),
                           "i":str(seg.xi[2]),
                           "l":str(seg.xi[3])})
            xml_line_nodes = SubElement(xml_line,"nodes")
            for vert in seg.vertices:
                SubElement(xml_line_nodes,"node",{"tag":str(vert.tag)})

        xml_entangled = SubElement(xml_root,"entangled_nodes")
        for enode in self.entangled_nodes:
            tmp = SubElement(xml_entangled,"node",{
                    "tag1":str(enode[0]),
                    "tag2":str(enode[1])})

        return xml_root

            


#########################################################################
# Main program
#  1. Read data from in-file
#  2. Create configuration
#  3. Build XML tree
#  4. Ouput to file
#########################################################################

def prettify(elem):
    from xml.etree import ElementTree
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ElementTree.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")


pers = ArgumentParser()

pers.add_argument("fin", help="Name of input DXA file")
pers.add_argument("-o","--out", help="Name of output file",
                  default="dx2graphout")

args = pers.parse_args()

### 1. Reading DXA file
#-------------------------
fin = open(args.fin)
print("... Reading DXA file: "+args.fin)
seg_section = []
box_section = []
cluster_section = []

while True:
    line = fin.readline()
    if line == "": break   #EOF
    line = line.strip().split()

    if line[0] == "SIMULATION_CELL_ORIGIN":
        box_section.append([float(tmp) for tmp in line[1:4]])
        box_section.append([float(tmp) for tmp in line[1:4]])
        continue
    if line[0] == "SIMULATION_CELL_MATRIX":
        for i in range(3):
            line = fin.readline().strip().split()
            box_section[1][i] += float(line[i])
        continue

    if line[0] == "CLUSTER":
        clust = [int(line[1])]     #cluster id
        line = fin.readline().strip().split()
        clust.append(int(line[1])) #cluster structure
        line = fin.readline()
        for i in range(3):
            line = [float(x) for x in fin.readline().strip().split()]
            clust.append(line)     #rotation matrix
        cluster_section.append(clust)
        continue
    
    if line[0] == "DISLOCATIONS":
        total_num_segments = int(line[1])
        for i in range(total_num_segments):
            seg_section.append([])
            seg_section[-1].append(int(fin.readline().strip()))
            seg_section[-1].append([float(x) for x in fin.readline().strip().split()])
            seg_section[-1].append(int(fin.readline().strip()))
            seg_section[-1].append(int(fin.readline().strip()))
            for ivert in range(seg_section[-1][3]):
                seg_section[-1].append([float(x) for x in fin.readline().strip().split()])
        continue

print("... Done reading file")
print("    +++ simulation cell:")
for i in range(2):
    print("        "+str(box_section[i]))
print("    +++ %i segments, %i clusters"%(len(seg_section),len(cluster_section)))

cnfg = configuration(seg_section,cluster_section,box_section)
print("... Configuration built, %i distinct nodes"%(cnfg.count))

print("\n----------------NODE LIST-----------------")
for inode in range(len(cnfg.nodes)):
    if cnfg.repeated[inode] == 1: continue
    else: print(str(cnfg.nodes[inode]))
print("-------------END NODE LIST-----------------")

print("\n-------------SEGMENT LIST------------------")
for iseg in range(len(cnfg.segments)):
    print(str(cnfg.segments[iseg]))
print("-----------END SEGMENT LIST----------------")

xmldata = cnfg.build_xml_tree()

#print(prettify(xmldata))

