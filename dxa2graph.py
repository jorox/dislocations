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

REC_STRUCTURE_TYPES = {"HCP":2}
class mb:
    """
    The class is meant to handle Miller-Bravais style vectors
    Its methods are still unimplemented
    To complete the code we need to finish all the arithmetic functions: add, divide, multiply
    and the inner, outer product functions
    For simplicity and efficiency these functions are not yet implemented
    """
    def __init__(self,u,v,w,t):
        self.coords = [u,v,w,t]
    def __str__(self):
        strtmp = "[%i%i%i%i]"%(self.coords[0],self.coords[1],self.coords[2],self.coords[3])
        return strtmp

    def norm(self):
        nrm = 0.0
        for i in range(4): nrm += self.coords[i]**2
        nrm -= 0.5*(self.coords[0]*self.coords[1]+
                    self.coords[0]*self.coords[2]+
                    self.coords[1]*self.coords[2])

    def __add__(self,other):
        res = self.scale*np.array(self.coords)+other.scale*np.array(other.coords)
    
        
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
    def __init__(self,index,localb,cluster_id,vertices):
        self.index = index
        self.localb = [x for x in localb]
        self.cluster_id = cluster_id
        self.vertex = list(vertices)
    def __str__(self):
        tmps = "segment %i:\n%1.3f %1.3f %1.3f\ncluster = %i\n%i vertices"%(
            self.index, self.localb[0], self.localb[1], self.localb[2], self.cluster_id, len(self.vertex))
        for i in range(len(self.vertex)):
            tmps += "\n%1.3f %1.3f %1.3f"%(self.vertex[i][0],self.vertex[i][1],self.vertex[i][2])
        return tmps
    def __len__(self):
        return len(self.vertex)

    def wrap_vert(self,origin,scvec):
        res = 0
        for vert in self.vertex:
            for i in range(3):
                uplimit = scvec[i][i]+origin[i]
                if vert[i] < origin[i]:
                    vert[i]= uplimit - (origin[i]-vert[i])
                    res = -1
                elif vert[i] > uplimit:
                    vert[i] = origin[i] + (uplimit-vert[i])
                    res = 1
        return res


pers = ArgumentParser()

pers.add_argument("fin", help="Name of input DXA file")
pers.add_argument("-o","--out", help="Name of output file")

args = pers.parse_args()

### Reading DXA file
#----------------------------------------------------------------
fin = open(args.fin)
print("... Reading DXA file: "+args.fin)
clusters = []
list_disloc = []
list_junc = []
sc_origin = []
sc_matrix = []

while True:
    line = fin.readline()
    if line == "": break #EOF
    line = line.strip().split()
    
    # Simulation cell dimensions
    if line[0] == "SIMULATION_CELL_ORIGIN":
        sc_origin = [float(tmp) for tmp in line[1:4]]
        print("+++Simulation cell origin")
        print("%1.4f %1.4f %1.4f"%(sc_origin[0], sc_origin[1], sc_origin[2]))
        continue
    if line[0] == "SIMULATION_CELL_MATRIX":
        for i in range(3):
            line = fin.readline().strip().split()
            sc_matrix.append([float(x) for x in line])
        print("+++Simulation cell matrix")
        print(sc_matrix)
        continue
            

    # Clusters
    if line[0] == "CLUSTER":
        cluster_mat = []
        cluster_index = int(line[1]) # cluster index
        line = fin.readline().strip().split() #cluster structure
        cluster_stype = int(line[1])
        line = fin.readline() #CLUSTER_ORIENTATION
        for i in range(3):
            line = fin.readline().strip().split()
            cluster_mat.append([float(tmp) for tmp in line])
        if cluster_stype == REC_STRUCTURE_TYPES["HCP"]:
            clusters.append(cluster(cluster_index,"HCP",cluster_mat))
            print("    \n+++cluster added")
            print(str(clusters[-1]))
        continue

    # Dislocations
    total_num_vertices = 0
    if line[0] == "DISLOCATIONS":
        total_num_dislocs = int(line[1])
        for i in range(total_num_dislocs):
            dis_index = int(fin.readline().strip())   #index
            line = fin.readline().strip().split()
            dis_lb = [float(tmp) for tmp in line]   #local Burgers vector
            dis_cluster = int(fin.readline().strip())         #cluster
            dis_size = int(fin.readline().strip())            #size = number of vertices
            dis_vertices = []                                 #vertices
            for j in range(dis_size):
                line = fin.readline().strip().split()
                dis_vertices.append([float(tmp) for tmp in line])
            list_disloc.append(segment(dis_index,dis_lb,dis_cluster,dis_vertices))
            
        print("    \n+++loaded %i dislocations"%(len(list_disloc)))
        for d in list_disloc:
            print(str(d))
            total_num_vertices += len(d)
        print("-------------- end dislocation list ------------------")
        
        total_num_dislocs = len(list_disloc)
        print("    +++ %i vertices for %i segments"%(total_num_vertices,total_num_dislocs))
        continue
    
    # Dislocation Junctions
    if line[0]=="DISLOCATION_JUNCTIONS":
        print("    +++ junctions:")
        for i in range(2*total_num_dislocs):
            line = fin.readline().strip().split()
            list_junc.append([int(tmp) for tmp in line])
            print("     "+str(line))
        print("-------------- end junctions list -------------------")
        continue
        
### Create Distinct Node List
list_node = []
for i_dis in range(len(list_disloc)):
    d = list_disloc[i_dis]
    for i_vert in range(len(d.vertex)):
        skip_vert = False
        if (i_dis>0 and i_vert==0) or (i_dis>0 and i_vert==1):
            for i_junc in range(2*i_dis):
                if list_junc[i_junc] == [i_vert,i_dis]: skip_vert = True
        if skip_vert: continue
        list_node.append(list(d.vertex[i_vert]))

print("\n... %i distinct nodes for %i segments"%(len(list_node),len(list_disloc)))

### Determine the true Burgers vector for each dislocation
##  Calculate the Spatial Burgers Vector by multiplying the
##  local Burgers vector by the matrix of the cluster
##  divide the spatial-b by the lattice spacings along the 3 main directions
X = [1./3,1./3,-2./3,0.]; Z = [-1.,1.,0.,0.]; Y = [0.,0.,0.,1.]
box_spc = [3.232,5.165,3.232*np.sqrt(3)]
print("\n... principla directions and spacings along")
print("     X = "+str(X)+"  "+str(box_spc[0])+" A")
print("     Y = "+str(Y)+"  "+str(box_spc[1])+" A")
print("     Z = "+str(Z)+"  "+str(box_spc[2])+" A")
X = np.array(X); Y = np.array(Y); Z = np.array(Z)
list_trueb = []
list_spceb = []
print("\n... True Burgers vectors:")
for d in list_disloc:
    found_cluster = False
    for clstr in clusters:
        if d.cluster_id == clstr.index:
            found_cluster = True
            spceb = clstr.A.dot(d.localb)              # left-multiply local Burgers by rotation matrix
            for i in range(3):  spceb[i] /= box_spc[i] # get scaled-spatial Burgers vector
            trueb = X*spceb[0] + Y*spceb[1] + Z*spceb[2] # true Burgers vector 
            print("     %i) -- 1/3[%1.3f,%1.3f,%1.3f,%1.3f]"%(
                d.index, 3*trueb[0],3*trueb[1],3*trueb[2],3*trueb[3]))
            #print("     "+str(spcb))
            list_trueb.append(trueb)
            list_spceb.append(spceb)
            break #cluster found - stop searching
        #if cluster has not been found
    if found_cluster == False:
            print("ERROR: Could not find cluster %i for segment %i"%(d.cluster_id,d.index))
                
print("--------------- end True Burgers section ------------------")


### Determine glide plane for each dislocation
##  Currently only prismatic and basal glide planes are defined
##  The sense vector is defined by calculating the vector that joins the first and
##  final points
print("\n... Glide plane directions")
list_gplane = []
i = 0
for d,spcb in zip(list_disloc,list_spceb):
    tmp = np.diff(d.vertex,axis=0)
    tmp = np.sum(tmp,axis=0)
    tmp = tmp/len(d.vertex)
    for i in range(3): tmp[i] /= box_spc[i]
    tmp = np.cross(tmp,spcb)
    tmp = 3.*X*tmp[0] + 3.*Y*tmp[1] + 3.*Z*tmp[2]
    tmp = tmp/np.max(np.abs(tmp))
    print("    %i) -- "+ str(tmp))
    i +=1

### Wrap vertices for dislocations passing through a periodic boundary
##
print("\n... Wrapping dislocations: +/-1 wrapped positive/negative, 0 no wrapping")
for d in list_disloc:
    res = d.wrap_vert(sc_origin,sc_matrix)
    print("    %i) -- %i"%(i,res))
    i+=1
                
            
        
