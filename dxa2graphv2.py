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
    def __init__(self,list_nodes,list_segments,list_junctions,topo):

        #
