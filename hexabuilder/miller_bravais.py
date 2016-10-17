import numpy as np

class Point(object):
    """A class representing a point in Miller-Bravais notation. 
    Points have the following properties:

    Attributes:
    coords: A numpy array of 4 numbers
    """
    def __init__(self,u,v,w,t):
        self.coords = np.array([u,v,w,t])
    def __init__(self):
        self.coords = np.array([0,0,0,0])
    def __init__(self, a):
        if len(a) != 4:
            raise ValueError('Input array must have exactly 4 elements')
        self.coords = np.array(a)
    
    def __str__(self):
        return "[%1.3f, %1.3f, %1.3f, %1.3f]"%(self.u,self.v,self.w,self.t)

class Vector(Point):
    """A class representing a vector or direction in Miller-Bravais notation.
    DVector has the following properties:

    Attributes:
    coords: A numpy array of 4 numbers
    """
    def __init__(self,u,v,w,t):
        self.coords = np.array([u,v,w,t])
    def __init__(self):
        self.coords = np.array([0,0,0,0])

    def __str__(self):
        return "[%1.3f, %1.3f, %1.3f, %1.3f]"%(self.u,self.v,self.w,self.t)

    def __getitem__(self,key):
        return self.coords[key]
    def __setitem__(self,key,value):
        self.coords[key] = value
    def __iter__(self):
        self.n = 0
        return self.coords[0]
    def __next__(self):
        if self.n<4:
            self.n += 1
            return self.coords[n]
        else:
            raise StopIteration
    def __add__(self,other):
        return Vector(self.coords + other.coords)
    def __subtract__(self,other):
        return Vector(self.coords-other.coords)
    def __mul__(self,a):
        if isinstance(a, (int,long,float)):
            return Vector(a*self.coords)
        elif isinstance(a, Vector):
            return self.dot(a)
        else:
            raise ValueError("Error: Multiplication cannot be performed")
    def __rmul__(self,a):
        if isinstance(a, (int,long,float)):
            return Vector(a*self.coords)
        elif isinstance(a, Vector):
            return self.dot(a)
        else:
            raise ValueError("Error: Multiplication by a non-number")
    
    def dot(self, A):
        """Returns the dot product between two Miller-Bravais vectors
        """
        if isinstance(A,Vector):
            tmp=(self.u*A[0]+self.v*A[1]+self.w*A[2]+self.x*A[3]
              -0.5*(self.u*A[1]+self.u*A[2])
              -0.5*(self.v*A[0]+self.v*A[2])
              -0.5*(self.w*A[0]+self.w*A[1]))
            return tmp
        else:
            raise ValueError(
                "Error: cannot perform inner-product on non Miller-Bravais 
                vector")
    
    def norm(self):
        """This function returns the norm of the vector calculated according
        to the standard non-orthogonal Hexagonal basis
        """
        return self.dot(Vector(self))**0.5

        
    
