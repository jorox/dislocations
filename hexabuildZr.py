import numpy as np
import sys

ver = "01"

class mbvec:
    def __init__(self,v):
        self.u=v[0]
        self.v=v[1]
        self.w=v[2]
        self.x=v[3]

    def norm(self):
        # calculates the norm of a Miller-Bravais vector
        tmp = np.sqrt(self.u**2+self.v**2+self.w**2+self.x**2
                      -self.u*self.v-self.u*self.w-self.v*self.w)
        return tmp
    
    def mb2m(self):
        # returns a Numpy array with the Miller representation of the vector
        tmp = [self.u-self.w, self.v-self.w, self.x]
        return np.array(tmp)

    def coord(self):
        # returns a list with Miller-Bravais coordinates
        # CONSIDER CHANGING TO __getitem__
        return [self.u,self.v,self.w,self.x]
    
    def printme(self):
        # a function to return a human-readable version of the vector
        print("%1.4f %1.4f %1.4f %1.4f"%(self.u,self.v,self.w,self.x))
    
    def dot(self,other):
        # calculates the dot product of vector with another MB vector
        A = other.coord()
        tmp=(self.u*A[0]+self.v*A[1]+self.w*A[2]+self.x*A[3]
              -0.5*(self.u*A[1]+self.u*A[2])
              -0.5*(self.v*A[0]+self.v*A[2])
              -0.5*(self.w*A[0]+self.w*A[1]))
        return tmp
    
    def length(self,a,c):
        # returns the length in whatever units a and c are in
        tmp = np.sqrt(a**2*(self.u**2+self.v**2+self.w**2 
                         -self.u*self.v-self.u*self.w-self.v*self.w)
                     +c**2*self.x**2)
        return tmp

    def __add__(self,other):
        # adds two MB vectors
        a = other.coord()
        tmp = mbvec([self.u+a[0],self.v+a[1],self.w+a[2],self.x+a[3]])
        return tmp
    def __sub__(self,other):
        # subtracts two MB vectors
        a = other.coord()
        tmp = mbvec([self.u-a[0],self.v-a[1],self.w-a[2],self.x-a[3]])
        return tmp
    def __mul__(self,other):
        # multiplication by a scalar
        tmp = mbvec([other*self.u,other*self.v,other*self.w,other*self.x])
        return tmp
    def  __rmul__(self,other):
        # Right-multiplication by a scalar
        tmp = mbvec([other*self.u,other*self.v,other*self.w,other*self.x])
        return tmp

class crystalbasis:
    motif = [np.array([0,0,0]),np.array([2./3,1./3,1./2])]
    def __init__(self):
        #The main directions of the crystal are:
        #  a1 = 1/3[2-1-10], a2 = 1/3[-12-10], a3 = [0001]
        # ==> motif = = [0,0,0]; [2/3,1/3,0.5]
        #The initializer sets te X, Y, Z to
        # X = [11-20], Y=[-21-10], Z
        X = mbvec(1./3*np.array([1,1,-2,0]))
        Y = mbvec(1./3*np.array([-2,1,-1,0]))
        Z = mbvec
        

def banner(text, ch='=', length=78):
    spaced_text = ' %s ' % text
    banner = spaced_text.center(length, ch)
    return banner

def main():
    print("\n"+banner("Zr-builder ver."+ver,length=100)+"\n")
    fin = open(sys.argv[1])
    flog = open("build.log",'w')
    print("... reading input file " + sys.argv[1])
    flog.write("... reading input file " + sys.argv[1]+'\n')
    
# determine the direction of XYZ
    for line in fin:
        line = line.strip().split()
        # if line[0]=="X" or line[0]=="x":
        #     print("     >>found new X direction ")
        #     flog.write("     >>found new X direction \n")
        #     X=mbvec(np.array(line[1:],dtype='float'))
        # if line[0]=="Y" or line[0]=="y":
        #     print ("     >>found new Y direction ")
        #     Y=mbvec(np.array(line[1:],dtype='float'))
        # if line[0]=="Z" or line[0]=="z":
        #     print("     >>found new Z direction ")
        #     Z=mbvec(np.array(line[1:],dtype='float'))
        if line[0]=="nx":
            print("      ++found x-limits")
            nx = np.array(line[1:3],dtype='int')
        if line[0]=="ny":
            print("      ++found y-limits")
            ny = np.array(line[1:3],dtype='int')
        if line[0]=="nz":
            print("      ++found z-limits")
            nz = np.array(line[1:3],dtype='int')
        if line[0]=="edge":
            print("      ++edge dislocation")
            createedge = True
        if line[0] == "screw":
            print("      ++screw dislocation")
            createscrew = True
        if line[0]=="loop":
            createloop=True
            loop_center = int(line[1])
            n0001 = [0, 1]
            n1100 = [0, 1]
            print("      ++loop, centered "+str(loop_center)+", unit size")
        if line[0]=="n0001":
            print("      ++loop n0001 direction changed")
            n0001 = np.array(line[1:3],dtype='int')
        if line[0]=="n1100":
            print("      ++loop n1100 direction changed")
            n1100 = np.array(line[1:3],dtype='int')
        if line[0]=="Temp" or line[0]=="temp" :
            print("      ++found temperature")
            T = int(line[1])
    
    if T==10:
        a = 3.23412 
        c = 5.16798 #calculated at 10K and 300K by minimizing the potential energy
    if T == 300:
        a = 3.23289
        c = 5.17302

            
# !Need to check that directions are orthogonal!
# !Need to check that limits do not coincide!

    X = 

    print("\n... orienting crystal");flog.write("... orienting crystal\n")
    print("     X/Y/Z directions: "); flog.write("     X/Y/Z directions: \n")
    X.printme();Y.printme();Z.printme()
