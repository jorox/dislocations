"""
This program transforms a DXA file to a graph file for DD simulation
Assumptions:
 ** Working for HCP structure only
     Other structures are not considered
 ** Box is assumed to be orthogonal even for screw dislocation simulations
     This is to help in wrapping vertices across the periodic boundaries
 ** Main directions are set
     X=[11-20], Z=[-1100], Y=[0001]
 *** Only one type of stacking fault is recognized I2 (basal)
 *** Need to deal with glide plane for screw dislocations
 *** Need to with hikl-->hkil (correct) in XML output
 *** Glide plane normals need to be integers or decimal free floats
 *** Partial burgers vectors are identified by significant decimal component in any of their spatial Burgers
     vector componenets
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
    def __init__(self,index,b_coords,xi_coords,vertices,ispart):
        self.index = index
        self.b = mb(b_coords)
        self.xi = mb(xi_coords)
        self.vertices = list(vertices)
        self.ispartial = ispart
        
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
        self.partials = []
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
        self.X = X; self.Y = Y; self.Z = Z
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
            ispart = False                                          # non-perfect Burgers
            rotmat = np.array(list_clusters[idcluster][2:5])        # Rotation matrix for Burgers
            spcb = rotmat.dot(tmp)                                  # left-mult local Burgers by rotation matrix
            for i in range(3):  spcb[i] = spcb[i]/box_spc[i]        # scale vector (no units)
            for ispcb in range(3):
                if np.abs(spcb[ispcb]-int(spcb[ispcb])) > 0.05: ispart = True
            trueb = 6.0*(X*spcb[0] + Y*spcb[1] + Z*spcb[2])         # transform to mb notation
            self.partials.append(ispart)

            # Determine glide plane normal
            sens = np.diff(d[4:4+N_verts],axis=0)                      # calculate partial sense vectors
            sens = np.sum(sens,axis=0)                                 # sum
            sens = sens/N_verts                                        # average
            b2pr = (np.sqrt(2)*sens[1])/np.sqrt(sens[0]**2+sens[2]**2) # basal/prismatic ratio
            if b2pr < 1.0: sens[1] = 0
            normsens = np.linalg.norm(sens)
            for i in range(3): sens[i] = sens[i]/box_spc[i]    # scale vector (no units)

            print("########## DEBUG ###########")
            print("sense = " + str(sens))
            print("spcb = " + str(spcb))
            print(b2pr)
            truexi = np.cross(spcb,sens)
            print("cross = " + str(truexi))
            truexi /= np.linalg.norm(truexi)
            print("norm corss = " + str(truexi))
            fac = (truexi[0]+truexi[2])/(truexi[0]-truexi[2])
            print(fac)
            truexi = 3.0*(X*truexi[0] + Y*truexi[1] + Z*truexi[2])      # transform xi to mb
            for itx in range(4): truexi[itx] = round(truexi[itx],0)
            if truexi[3]<0 : truexi[3] = -truexi[3]
            print(truexi)
            print("########## DEBUG ###########")

            # Vertex list
            truevertex = list(self.nodes[ivert:ivert+N_verts])
            self.segments.append(segment(d[0],trueb,truexi,truevertex,ispart))
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
                            self.segments.append(segment(len(self.segments),d.b,d.xi,d.vertices[iv:],d.ispartial))
                            self.partials.append(d.ispartial)
                            self.entangled_nodes.append([d.vertices[iv-1],d.vertices[iv]])
                            print("   ****splitting segment %i at node #%i-tag:%i"%(
                                d.index,iv,d.vertices[iv].tag))
                            d.vertices = d.vertices[:iv]
                            break

        #############################################################################
        # wrap nodes across periodic boundaries
        #############################################################################
        for ix in range(3):
            for vertex in self.nodes:
                if vertex[ix]<topo[0][ix]: vertex[ix] = vertex[ix]-topo[0][ix]+topo[1][ix] # x<origin
                elif vertex[ix]>topo[1][ix]: vertex[ix] = vertex[ix]-topo[1][ix]+topo[0][ix] # x>limit

        print("... Done building configuration")
        
    def fix_nodes_plane(self):
        """
        The function fixes the coordinates of the nodes belonging to each segment to all be
        in the same glide plane.
        In MD the nodes can be slightly outside the same glide plane (in parallel planes).
        Numodis will block in such a case
        The function will calculate the equation of the plane using the normal to the plane
        i.e. the xi[] vector and a point. The point chosen will be the average of the nodes 
        belonging to the segment. 
        if xi = [a,b,c] and p=[x0,y0,z0] then plane s: 0=ax+by+cz+d where d = -ax0-by0-cz0
        Then say we want the projection of point A=[x1,y1,z1] to A'=[x2,y2,z2]
        we can say that the line (AA') has the parametric equation: (x=x1+mt, y=y1+nt, z=z1+pt)
        but the direction (AA') is parallel to xi. Hence, the coordinates of A' are
        (x2=x1+at, y2=y1+bt, z2=z1+ct).
        To find "t" we note that the coordinates of A' satisfy the equation of the plane:
        t = - (d+ax1+by1+cz1)/(a**2+b**2+c**2) = - (a(x1-x0)+b(y1-y0)+c(z1-z0))/(a**2+b**2+c**2)
        """
        # get the a1,a2,a3,a4 in terms of XYZ to transform xi to XYZ rep.
        a1 = np.array([2,-1,-1,0]); a2=np.array([-1,2,-1,0]); a3=np.array([-1,-1,2,0]); a4=np.array([0,0,0,1])
        repmat = np.array([self.X,self.Y,self.Z])
        repmat = np.transpose(repmat)
        repmat = np.linalg.inv(repmat) # [X|Y|Z]a = [ ]
        a1xyz = np.dot(repmat,a1); a2xyz=np.dot(repmat,a2); a3xyz=np.dot(repmat,a3); a4xyz=np.dot(repmat,a4)
        print("########### DEBUG ###########")
        print("a1 = "+str(a1xyz))
        print("a2 = "+str(a2xyz))
        print("a3 = "+str(a3xyz))
        print("a4 = "+str(a4xyz))
        print("##############################")
        for seg in self.segments:
            # find xi in XYZ rep
            n = a1xyz*xi[0] + a2xyz*xi[1] + a3xyz*xi[2] + a4xyz*xi[3]
            p0 = np.array([0.,0.,0.])
            for tmp_node in seg.vertices:
                p0 += np.array(tmp_node.coords)
            p0 /= len(seg.vertices)
            for iv in range(len(seg.vertices)):
                p1 = np.array(seg.vertices[iv].coords)
                tmp = p1-p0
                t = -(xi[0]*tmp[0]+xi[1]*tmp[1]+xi[2]*tmp[2])/(np.sum(xi**2))
                seg.vertices[iv].coords = list(p1+xi*t)
                print("############ DEBUG - fix_nodes_plane #############")
                print("old = "+str(p1))
                print("new = "+str(seg.vertices[iv].coords))
                print("##################################################")
        
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
                           "bh":str(int(seg.b[0])),
                           "bk":str(int(seg.b[1])),
                           "bi":str(int(seg.b[2])),
                           "bl":str(int(seg.b[3])),
                           "tag":str(seg.index)})
            xml_constraint = SubElement(xml_line,"constraints")
            
            tmp = list(seg.xi)
            for i in range(4): tmp[i] = round(tmp[i],0)   # round to nearest integer
            tmp[2] = -(tmp[0]+tmp[1])                     # Miller-Bravais
            
            xml_plane = SubElement(xml_constraint,"plane",{
                           "h":str(tmp[0]),
                           "k":str(tmp[1]),
                           "i":str(tmp[2]),
                           "l":str(tmp[3])})
            
            if seg.ispartial:
                #print("PARTIAL FOUND")
                xml_sf = SubElement(xml_line, "stackingfaults")
                SubElement(xml_sf,"sf",{
                    "h":str(tmp[0]), "k":str(tmp[1]), "i":str(tmp[2]), "l":str(tmp[3]),
                     "direction":str(tmp[3]), "type":"I2"})
                
            xml_line_nodes = SubElement(xml_line,"nodes")
            for vert in seg.vertices:
                SubElement(xml_line_nodes,"node",{"tag":str(vert.tag)})

        xml_entangled = SubElement(xml_root,"entangled_nodes")
        for enode in self.entangled_nodes:
            tmp = SubElement(xml_entangled,"node",{
                    "tag1":str(enode[0].tag),
                    "tag2":str(enode[1].tag)})

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
                  default="dx2graphout.xml")

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
fout = open(args.out,"w")
fout.write(prettify(xmldata))
print("... done writing output to %s"%(fout.name))
print("    %i nodes, %i segments of which %i partials"%(len(cnfg.nodes),len(cnfg.segments),np.sum(cnfg.partials)))
#print(prettify(xmldata))

