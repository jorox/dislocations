'''
=================
3D wireframe plot of a PBC dislocation dipole
=================

A very basic demonstration of a wireframe plot.
'''
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def u_dip(x,y,z,b,Rb,Rd):
    """
    return the displacement at (x,y,z) due to a dislocation dipole located at 
    Rb and Rd with burgers b,-b respectively
    """
    ux = b/2.0/3.1415*(np.arctan2(z-Rb[2],y-Rb[1])
                      -np.arctan2(z-Rd[2],y-Rd[1]))
                      
    return np.array([ux,0.,0.])

def u_sum(x,y,z,b,Rb,Rd,cx,cy,cz,n):
    """
    return the displacement at (x,y,z) due to a dislocation dipole and 
    its images.
    cx, cy, cz are the image translation vectors.
    n is a vector specifying the number of images to include along each
    image translation vector
    """
    Rb = np.array(Rb)
    Rd = np.array(Rd)
    cx = np.array(cx)
    cy = np.array(cy)
    cz = np.array(cz)
    tmp = np.array([0.,0.,0.])
    
    n[np.where(np.isnan(Rb))[0][0]] = 0 #np.where returns a tuple of arrays
    
    for ix in range(-n[0],n[0]+1,1):
        for iy in range(-n[1],n[1]+1,1):
            for iz in range(-n[2],n[2]+1,1):
                #print(ix,iy,iz)
                tmp += u_dip(x,y,z,b,Rb+ix*cx+iy*cy+iz*cz,Rd+ix*cx+iy*cy+iz*cz)

    return tmp


def u_pbc(x,y,z,b,Rb,Rd,cx,cy,cz,n,origin,u0,d=None):
    """
    return the displacement at (x,y,z) due to a dislocation dipole and its 
    periodic images with correction to restore periodic boundary conditions
    origin is coordinates of the central image.
    
    u_pbc = u_sum - d.r -u0
    0 = u_sum(r0) - d.r0 - u0 where r0 is a point on the central-image limit   
    
    d is a 3x3 constant matrix
    
    To find d we need to solve 3 linear systems
    | x1 y1 z1 | |di1|    |u(r1)i|
    | x2 y2 z2 | |di2| =  |u(r2)i|   where i=1,2,3 
    | x3 y3 z3 | |di3|    |u(r3)i| 

    """
    # some prep
    Rb = np.array(Rb)
    Rd = np.array(Rd)
    cx = np.array(cx)
    cy = np.array(cy)
    cz = np.array(cz)
    origin = np.array(origin)
    u0 = np.array(u0)
    
    # get the 4 special points
    if (d==None):
        print("uerr pNonearameters not given, calculating...")
        rsp = []
        rsp.append(origin + cy)
        rsp.append(origin + cz)
        rsp.append(origin + cy + cz)
        #rsp.append(origin + cx + cy)
        #rsp = np.insert(rsp,3,1,axis=1) # insert a column of '1's to the end
        print(rsp)
        print(" determinant of mat = %1.5f"%(np.linalg.det(rsp)))
        
        # caluclate the values of u_sum at those points
        u_sum_sp = []
        for i in range(3):
            u_sum_sp.append(u_sum(rsp[i][0],rsp[i][1],rsp[i][2],b,Rb,Rd,cx,cy,cz,n))
    
        # calculate the inverse matrix 
        inv_rsp = np.linalg.inv(rsp)
        d = []
        #u0 = []
        
        for i in range(3):
            tmp = [u_sum_sp[j][i] for j in range(3)]
            tmp = np.array(tmp)
            du_row = np.dot(inv_rsp,tmp) 
            d.append(du_row[:3])
            #u0.append(du_row[3])
            # determine u_0 knowing that d.r+u_0 = [0,0,z] at central image 
            # limits
        d = np.array(d)
        
    tmp = u_sum(x,y,z,b,Rb,Rd,cx,cy,cz,n)
    tmp = tmp - np.dot(d,np.array([x,y,z])) -u0
    
    return tmp,d
    
def main():
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # 1. Define the full region (zero-th image + periodic)
    x_range = np.arange(-1,1.1,0.05); N = len(x_range)
    y_range = np.arange(-0.5,0.51,0.05); M = len(y_range)
    
  

    X = []
    Y = []
    for i in range(M): X.append(x_range) # X = MxN
    for i in range(N): Y.append(y_range) # Y = NxM
    X = np.array(X)
    Y = np.transpose(np.array(Y)) # Y = MxN
    Z = np.zeros((M,N))

    Rb = np.array([0.5,0,np.nan])  # location of first disloc
    Rd = np.array([-0.5,0,np.nan]) # location of last disloc
    derr = None
    u0err = np.array([0.,0.,0.])
    for i in range(M):
        for j in range(N):
            #shift = u_dip(X[i,j],Y[i,j],Z[i,j],1,Rb,Rd)
            
            #shift = u_sum(X[i,j],Y[i,j],Z[i,j],1,Rb,Rd,
            #            [2.,0.,0.],[0.,1.,0.],[0.,0.,1.],[5,10,0])
            
            shift, derr = u_pbc(X[i,j], Y[i,j], Z[i,j], 1, Rb, Rd,
                        [2.,0.,0.], [0.,1.,0.], [0.,0.,3.], [5,10,3],
                        [-1.,-0.5,0.0001], u0err,derr) 
            X[i,j] += shift[0]; Y[i,j] += shift[1]; Z[i,j] += shift[2]
            
            
    # Plot a basic wireframe.
    ax.set_xlim3d([-1 , 1])
    ax.set_ylim3d([-0.5,0.5])
    ax.set_zlim3d([-1, 1])
    ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

    plt.show()
    
    
if __name__=="__main__":
    main()


