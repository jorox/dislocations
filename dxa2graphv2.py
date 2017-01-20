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

class mb:
    """
    The class is meant to handle Miller-Bravais style vectors
    Its methods are still unimplemented
    To complete the code we need to finish all the arithmetic functions: add, divide, multiply
    and the inner, outer product functions
    For simplicity and efficiency these functions are not yet implemented
    """
    def __init__(self,data):
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
        return "%i %1.3f, %1.3f, %1.3f"%(self.tag,self.coords[0],self.coords[1],self.coords[3])
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
        self.b = mb(c_coords)
        self.xi = mb(xi_coords)
        self.vertices = list(vertices)
        
    def __str__(self):
        tmps = "segment %i:\n%1.3f %1.3f %1.3f\ncluster = %i\n%i vertices"%(
            self.index, self.localb[0], self.localb[1], self.localb[2], self.cluster_id, len(self.vertex))
        for i in range(len(self.vertex)):
            tmps += "\n%1.3f %1.3f %1.3f"%(self.vertex[i][0],self.vertex[i][1],self.vertex[i][2])
        return tmps
    
    def __len__(self):
        return len(self.vertex)


class configuration:
    def __init__(self,list_vertices,list_segments,list_junctions,topo):
        self.nodes = []
        #######################################################################
        # wrap vertices across boundaries and create the list of nodes
        # self.nodes contains references to node objects
        # The list is ordered by appearance of nodes in segment information
        # Repetition is handeled by having the same object address
        #######################################################################
        for vertex in list_vertices:
            # wrap vertices if needed across the boundaries
            data = [vertex[0],None,None,None]
            for ix in range(3):
                if vertex[ix+1]<topo[0][ix]: data[ix+1] = vertex[ix+1]-topo[0][ix]+topo[1][ix] # x<origin
                elif vertex[ix+1]>topo[1][ix]: data[ix+1] = vertex[ix+1]-topo[1][ix]+topo[0][ix] # x>limit
                else: data[ix+1] = vertex[ix+1]

            # check for repetition, a tolerance is added because some repeated nodes might get wrapped
            #  hence the coordinates are not exactly the same
            found = False
            if len(self.nodes)>0:
                for i_node in self.nodes:
                    if np.sum(np.abs(i_node-node(data))) < 0.1: #same point
                        self.nodes.append(i_node) #add the same address to the list
                        found = True
                        break
                    
            if found == False: self.nodes.append(node(coords)) #new node

        ######################################################################
        # Determine the true Burgers vector for each dislocation
        #  Calculate the Spatial Burgers Vector by multiplying the
        #  local Burgers vector by the matrix of the cluster
        #  divide the spatial-b by the lattice spacings along the 3 main directions
        ######################################################################
        self.segments = []
        X = [1./3,1./3,-2./3,0.]; Z = [-1.,1.,0.,0.]; Y = [0.,0.,0.,1.]
        box_spc = [3.232,5.165,3.232*np.sqrt(3)]
        print("\n... principla directions and spacings along")
        print("     X = "+str(X)+"  "+str(box_spc[0])+" A")
        print("     Y = "+str(Y)+"  "+str(box_spc[1])+" A")
        print("     Z = "+str(Z)+"  "+str(box_spc[2])+" A")
        X = np.array(X); Y = np.array(Y); Z = np.array(Z)
        
        for d in list_segments:
            tmp = np.array(d[1])
            tmp = np.array(list_clusters[d[2]][2:5]).dot(tmp) # left-multiply local Burgers by rotation matrix
            for i in range(3):  tmp[i] /= box_spc[i] # get scaled-spatial Burgers vector
            trueb = X*tmp[0] + Y*tmp[1] + Z*tmp[2] # transform to mb notation

            
            d.realb = [bi for bi in trueb]
            break #cluster found - stop searching
        #if cluster has not been found
    if found_cluster == False:
            print("ERROR: Could not find cluster %i for segment %i"%(d.cluster_id,d.index))
                
print("--------------- end True Burgers section ------------------")
        
