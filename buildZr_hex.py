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
#    motif = np.array([[0,0,0], [1./3,2./3,0.5], [0,1,0], [-2./3,2./3,0.5]]) # config2
    motif = np.array([[0,0,0], [-1./3,1./3,0.5], [0,1,0], [-1./3,4./3,0.5]]) # config1
#    motif = np.array([[-1,1,0], [-1./3,1./3,0.5], [0,1,0], [-1./3,4./3,0.5]]) #reverse C along [11-20]
    def __init__(self,X,Y,Z):
        e1 = mbvec([1,0,0,0])
        e2 = mbvec([0,1,0,0])
        e3 = mbvec([0,0,1,0])
        e4 = mbvec([0,0,0,1])
        # determine the coordinates of a1,2,3,4 in the XYZ basis
        self.a = [ [X.dot(e1)/X.norm(), Y.dot(e1)/Y.norm(), Z.dot(e1)/Z.norm()], 
                   [X.dot(e2)/X.norm(), Y.dot(e2)/Y.norm(), Z.dot(e2)/Z.norm()], 
                   [X.dot(e3)/X.norm(), Y.dot(e3)/Y.norm(), Z.dot(e3)/Z.norm()],
                   [X.dot(e4)/X.norm(), Y.dot(e4)/Y.norm(), Z.dot(e4)/Z.norm()] ]
        self.a = np.array(self.a)

        print("\n         !!  Debug Info - a  !!\n!!----------------------------!!")
        print(self.a)
        print("!!----------------------------!!\n")
    

        # determine the coordinates of the motif atoms in the XYZ basis
        self.motif_XYZ = []
        for i in range(len(self.motif)):
            tmp = (self.a[0]*self.motif[i,0]+
                   self.a[1]*self.motif[i,1]+
                   self.a[3]*self.motif[i,2])
            tmp /= [X.norm(), Y.norm(), Z.norm()]
            self.motif_XYZ.append(tmp)
        print("\n         !!  Debug Info - motif_xyz  !!\n!!----------------------------!!")
        print(self.motif_XYZ)
        print("!!----------------------------!!\n")
    
    def __getitem__(self,index):
        return self.a[index]
    def get_xyz_motif(self,i):
        return self.motif_XYZ[i%len(self.motif_XYZ)] #returns an np array

def banner(text, ch='=', length=78):
    spaced_text = ' %s ' % text
    banner = spaced_text.center(length, ch)
    return banner

def main():
    print("\n"+banner("Zr-HEXA-builder ver."+ver,length=100)+"\n")

# define common constants
    T = 10


    # default directions for the lab axes 
    X = mbvec(1./3*np.array([1,1,-2,0]))
    Z = mbvec(1./2*np.array([-1,1,0,0]))
    Y = mbvec([0,0,0,1]) #default directions
    nx = [-1, 1]
    ny = [-1, 1]
    nz = [-1, 1] #default sizes
    createedge = False
    createloop = False
    createscrew = False
    load=["\\","|","/","-"]


    # input file for lattice construction
    fin = open(sys.argv[1])
    print("... reading input file " + sys.argv[1])

# determine the direction of XYZ
    for line in fin:
        line = line.strip().split()
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

# determine the crystal orientation
    print("\n... orienting crystal")
    print("     X/Y/Z directions: ")
    X.printme();Y.printme();Z.printme()
    basis = [np.array([0,0,0]),np.array([0.0,0.5,2./3])]

# determine the spacings and number of box cells
    print("\n... setting up simulation box")
    xspacing = X.length(a,c);
    yspacing = Y.length(a,c);
    zspacing = Z.length(a,c);
    print("     spacings (A): %1.5f, %1.5f, %1.5f"%(xspacing,yspacing,zspacing))
    print("     number of unit cells in XYZ directions: %1.0d %1.0d %1.0d"
          %(np.diff(nx), np.diff(ny), np.diff(nz)))
    
# determine the box limits
    SQRT3 = np.sqrt(3)
    nx[0] = 0; ny[0] = 0; nz[0]=0
    box = [[nx[0]*xspacing,nx[1]*xspacing],
           [ny[0]*yspacing,ny[1]*yspacing],
           [nz[0]*zspacing,nz[1]*zspacing]]
    if nz[1] > nx[1]: print("warning box too skewed")
    box_tilt = np.array([0, -nz[1]/2*xspacing, 0])
    box = np.array(box)
    
    
    print("\n... box dimensions (A):")
    print("     X %1.5f --> %1.5f  %3d -> %3d = %1.5fnm\n     Y %1.5f --> %1.5f    %3d -> %3d = %1.4fnm\n     Z %1.5f --> %1.5f   %3d -> %3d = %1.5fnm\n     tilt %1.5f %1.5f %1.5f" %(
                box[0,0],box[0,1],nx[0],nx[1],np.diff(box[0]/10),
                box[1,0],box[1,1],ny[0],ny[1],np.diff(box[1]/10),
                box[2,0],box[2,1],nz[0],nz[1],np.diff(box[2]/10), 
                box_tilt[0], box_tilt[1], box_tilt[2] ))
# create atoms
    natoms = (nx[1]-nx[0])*(ny[1]-ny[0])*(nz[1]-nz[0])*len(basis) #4 basis atoms
    nghost = 0

    print("... expecting %1.0d atoms"%natoms)
    atoms = []
    ia = -1 # atoms id (0-based)
    for iz in range(nz[0],nz[1]):
        for iy in range(ny[0],ny[1]):
            for ix in range(nx[0],nx[1]):
                for ib in range(len(basis)): # 2 atoms per basis poin
                    tmp = np.array([ix,iy,iz])+basis[ib]+np.array([-iz*0.5,0,0]) 
                    #tmp[2] += 0.5
                    tmp[0] *= xspacing 
                    #flog.write(str(tmp)+"\n")
                    # Change to cartesian coordinates
                    tmp[1] *= yspacing
                    tmp[2] *= zspacing
                    #flog.write(str(tmp))
                    ia += 1
                    atoms.append([ia+1, tmp[0],tmp[1],tmp[2],[ix,iy,iz,ib]])
                    prcnt = float(ia)/natoms*100
                    print(load[ia%4]+ "     %d %%"%(prcnt))
                    sys.stdout.write("\033[F")
     
# create the edge dislocation
    if createedge:
        print("\n /////////////// Building Edge dislocation \\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        nghost = buildedge1120(atoms,nx,ny,nz,a,box,box_tilt)
        natoms -= nghost
        print("... %d atoms, %d ghost atoms"%(natoms,nghost))
        print("/\\\\\\\\\\/////////\\\\\\\\\\//////////\\\\\\\\\\/////////////")

# create the screw dislocation
    if createscrew:
        print("\n /////////////// Building Screw dislocation \\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        buildscrew(atoms,nx,ny,nz,a)
        fixids(atoms)
        box = [[nx[0]*xspacing,nx[1]*xspacing],
               [ny[0]*yspacing,ny[1]*yspacing],
               [nz[0]*zspacing,nz[1]*zspacing]]
        box = np.array(box)
        print("/\\\\\\\\\\/////////\\\\\\\\\\//////////\\\\\\\\\\/////////////")
    
# create the loop 
    if createloop:
        print("\n/////////////// Buidling SIA Loop \\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        buildloop(atoms,basis,a,c,nx,ny,nz,n0001,n1100,loop_center)
        natoms -= nghost
        print("/\\\\\\\\\\/////////\\\\\\\\\\//////////\\\\\\\\\\/////////////")
    

# write LAMMPS input file
    foutname = "zr_hex.data"
    flammps = open(foutname,"w")
    flammps.write("#ver. "+ver+" HCP-Zr [11-20],[0001],[-1100] + dislocation and loop, metal units, "+str(T)+"K")
    flammps.write("\n    %d atoms\n    %d atom types\n"%(natoms,1)) #N atoms, 1 type
    flammps.write("\n%1.6f %1.6f xlo xhi\n%1.6f %1.6f ylo yhi\n%1.6f %1.6f zlo zhi\n%1.6f %1.6f %1.6f xy xz yz\n"%
                  (box[0,0],box[0,1], box[1,0],box[1,1], box[2,0],box[2,1],
                   box_tilt[0],box_tilt[1],box_tilt[2]))
    flammps.write("\nMasses\n")
    flammps.write("\n    1 91.224\n")
    flammps.write("\nAtoms\n")
    ia = 0
    for i in range(len(atoms)):
       # print("DEBUG-atoms",i, atoms[i][4], "DEBUG-atoms")
        if atoms[i][0] == -1: continue #ghost atom
        ia += 1
        flammps.write("\n\t%d\t%d\t%1.6f %1.6f %1.6f"%(
                ia,1, atoms[i][1], atoms[i][2], atoms[i][3]))
    flammps.close()


    print("\n... done writing %1.0d atoms to file %s and lammps file"%(natoms,foutname))    
    
# Write a LAMMPS snapshot for OVITO visualzation
#
#         ITEM: TIMESTEP
#         0
#         ITEM: NUMBER OF ATOMS
#         3800
#         ITEM: BOX BOUNDS pp pp ss
#         -16.16 12.928
#         -25.825 25.825
#         -27.9955 26.1295
#         ITEM: ATOMS id type x y z 

    fsnap = open("zr_hex.ov"+ver,"w")
    fsnap.write("ITEM: TIMESTEP\n0")
    fsnap.write("\nITEM: NUMBER OF ATOMS\n%d"%(natoms))
    fsnap.write("\nITEM: BOX BOUNDS xy xz yz pp pp ss\n%1.5f %1.5f %1.5f\n%1.5f %1.5f %1.5f\n%1.5f %1.5f %1.5f"%
                      (box[0,0],box[0,1], box_tilt[0], 
                       box[1,0],box[1,1], box_tilt[1], 
                       box[2,0],box[2,1], box_tilt[2]) )
    fsnap.write("\nITEM: ATOMS id type x y z")
    ia = 0
    for i in range(len(atoms)):
        if atoms[i][0] == -1: continue #ghost atom
        ia +=1
        fsnap.write("\n%d 1 %1.4f %1.4f %1.4f"%(ia,atoms[i][1],atoms[i][2],atoms[i][3]))
    fsnap.close()
    print("... done writing to snapshot file zr_hex.ov"+ver)

                    
def buildedge1120(atoms,nx,ny,nz,a,box,box_tilt):
    # creates an edge dislocation
    #   b=1/3[11-20], n=1/3[-1100], y=0

    # First determine where to cut along the Z-axis
    zcut = np.sum(nz)/2 # example: nz = [-1 1]-->0 OR [-2 1]-->0 OR [-2 7]-->2
    xcut = nx[1]-1  # example: nx = [-1 1], xcut=0
    nremoved = 0

    

    #Remove Atoms:
    #--------------------------------------------------
    # start from iz = zcut
    #  in the hexagonal picture i.e the primitive lattice there are two atoms to be removed
    #  in the cubic picture there is a wrapping proceduce that will take place
    for iy in range(ny[0],ny[1]):
        ix = xcut
        for iz in range(zcut-1,nz[0]-1,-1):
            for ib in range(2):
                aid = returnid([ix,iy,iz,ib],nx,ny,nz,[0,2])
                print("DEBUG",'\n',atoms[aid],'-',aid,'\n','DEBUG')
                atoms[aid][0] = -1; atoms[aid][1] = -9999999; atoms[aid][2] = -9999999;atoms[aid][3] = -9999999;
                nremoved += 1
    


    #Fix box boundaries
    #---------------------------------------------------
    N = nx[1]-nx[0]
    exx2 = 0.5/(N-1)
    exx1 = -0.5/N
    xspacing2 = a*(1.+0.5/(N-1))
    xspacing1 = a*(1.-0.5/N)

    box[0,0] = xspacing1*nx[0]
    box[0,1] = xspacing1*nx[1]
    box_tilt[1] = 
    
    L = box[0][1]-box[0][0]

    #Shift and wrap
    #---------------------------------------------------
    for i in range(len(atoms)):
        if atoms[i][0] == -1: #ghost atom
            continue
        if atoms[i][4][2] < zcut:
            atoms[i][1] += atoms[i][1]*exx2
         #   atoms[i][1] += nx[0]*(xspacing1-xspacing2)
            #if atoms[i][1] > box[0][1]: # some atoms need to be wrapped; these correspond to atoms on the 1/6(11-20) plane
            #    atoms[i][1] = atoms[i][1]-box[0][1] + box[0][0]
                #!!!! NEED TO FIX IX HERE !!!
            #elif atoms[i][1] < box[0][0]: # some atoms here might also need to be wrapped, but I'm not sure
            #    atoms[i][1] = box[0][0]-atoms[i][1]
                #!!!! NEED TO FIX IX HERE !!!

            #atoms[i][1] += nx[0]*(xspacing1-xspacing2) #put first planes into coincidence
        else:
            atoms[i][1] += atoms[i][1]*exx1 # no need to wrap these atoms they are in their proper place


                    
    print("... removed %1.0d atoms and rescaled the rest and the box"%(nremoved))
    print("      new box dimension along X (A): %1.4f %1.4f"%(box[0][0], box[0][1]))
    return nremoved

    #strain in top region is e_top=-1/(2N), in bottom region e_bot=1/2(N-1)

def buildscrew(atoms,nx,ny,nz,a):
    # The function creates a screw dislocation in the crystal defined by atoms
    # The method used is that taken from Serra Bacon 2013 and Khater Bacon 2010
    # The displacement field is given by 
    #    ux = (b/2/pi) (arctan(z/y)-pi) = n/2/pi (theta-pi) theta = 0 to 2pi
    # For y=0, theta can be either pi/2 or 3pi/2 or 0,
    #  Intro to Dislocations p.66 that for x=0 check the drawing
    for i in range(len(atoms)):
        if atoms[i][2] == 0: # y=0 
            if atoms[i][3] > 0 : #z= positive ==> theta = pi/2
                ux = a/2/np.pi*(np.pi/2-np.pi)
            else: 
                if atoms[i][3] == 0: ux = 0
                else: ux = a/2/np.pi*(3.0/2*np.pi-np.pi) #z=negative ==> theta=3pi/2
 
        else: ux = a/2/np.pi*(np.arctan(atoms[i][3]/atoms[i][2])-np.pi)
        
        atoms[i][1] += ux

    # To restablish periodicity across the Y-direction a displacement of +/-b/2 
    # was added to the x-coordinate of atoms on the +/- y boundaries
    for ix in range(nx[0],nx[1]):
        for iz in range(nz[0],nz[1]):
            for ib in range(4): # 4-atoms per box cell
                aid = returnid([ix,ny[0],iz,ib],nx,ny,nz, [0,4])
                atoms[aid][1] -= a/2 # -y boundary
                aid = returnid([ix,ny[1]-1,iz,ib], nx, ny, nz, [0,4])
                atoms[aid][1] += a/2 # +y boundary

def buildloop(atoms,basis,a,c,nx,ny,nz,n0001,n1100,xloop):
    # This function creates a dislocaiton loop with "nsia" ISAs
    # The loop has a Burgers vector [-12-10]
    # It follows the procedure as described by Bacon and Osetsky and de Diego in the papers:
    #    "Structure and Properties of vacancy and interstitial clusteris in a-Zirconium"
    #      J. Nucl. Mat. 374 p.87 (2008)
    #    " Mobility of interstitial clusters in a-zirconium"
    #      Metallurgical and Mater. Trans. A 33A p.783 (2002)
    # n0001 and n1100 are the number of stacks along the 0001 and -1100 directions
    # xloop : coordinates ix,iy,iz of the loop center
    #
    # Bugs:
    # -------------------------------
    # The function has a problem when creating loops with a center near the edge. Mainly
    #  when xloop is equal to n[x]-2. Perioding image wrapping is necessary
    
    loop_atoms = []

# Find the atoms which belong to loop-center
# Remember that atoms have been displaced in the x-direction above and below the zcut plane
# 3 of these atoms belong to the center loop cell
# They are the second, third, and fourth box-unit basis atoms, check basis.motif
# The fourth sia is the 1st motif atoms in the next x-cell ie ix+1,iy,iz
# Hence the 4 base atoms are
#  ca0=[ix,iy,iz,0]   ca1=[ix,iy,iz,1]   ca2=[ix, iy, iz, 2] ca3=[ix,iy,iz+1,3]
#  Then for one stack along the +/-[0001] direction we get 4 more such that
#    cai[1]+/-1 ...
# Then for one stack along the +/-[-1100] direction Z-direction
#    cai[2]+/-1
# The base atoms for the loop are the same as those in the paper. If one chooses a different motif
#   such as [0,0,0, ][1/3,2/3,0.5] the loop reverses it's Burgers vector
# The shift factor seems to affect the Burgers vector of the loop detected by DXA. For 0.1<eps<0.3
#  the Burgers vector direction is opposite to the disloc for eps>=0.3 it is the same  
    xloop = wrap(xloop,nx)

    for iy in range(n0001[0],n0001[1]):
        for iz in range(n1100[0],n1100[1]):
            iy = wrap(iy,ny)
            iz = wrap(iz,nz)
            catoms = [ [xloop,iy,iz,0], [xloop,iy,iz,1], 
                       [xloop,iy,iz,2], [xloop,iy,iz,3]]
            #print("DEBUG-catoms\n ------------- \n",catoms,"\n ------------\n DEBUG-catoms")
            for ic in catoms:
                loop_atoms.append(returnid(ic,nx,ny,nz,[0,4]))
                #print("DEBUG-loop \n",atoms[loop_atoms[-1]][4],"\n DEBUG-loop")
                if atoms[loop_atoms[-1]][0] == -1:
                    print("Warning: trying to modify ghost atom")

    Nloop = len(loop_atoms)
    last_index = atoms[-1][0]
    print("    +found %d atoms for the loop"%(Nloop))
    #print("DEBUG",loop_atoms,"DEBUG")

# Now that we've identified the atoms we need to add N with the same coordinates and then shift half of them 
#  left and half right along the [11-20] direction. I.E add +/-eps*(a1+a2) to their cartesian coordinates
#!! IMPORTANT Some crowdion atoms might change their box unit association because of the shift

    eps = 0.3 # shift magnitude
    shift =  np.array(eps*a*(basis[0]+basis[1]))
    print("DEBUG-shift",shift,"DEBUG-shift")

    for i in range(Nloop):
        last_index +=1
        oldatom = atoms[loop_atoms[i]]
     #   print("DEBUG-old",oldatom,"DEBUG-old")
        atoms.append([oldatom[0],oldatom[1],oldatom[2],oldatom[3],
                      [oldatom[4][0],oldatom[4][1],oldatom[4][2],oldatom[4][3]]])  #store old position
        for j in range(1,4):
            atoms[loop_atoms[i]][j] +=  shift[j-1] # add eps*a2 to atom 
            atoms[-1][j] -= shift[j-1] # add -eps*a2 to last atom
        atoms[-1][0] = last_index # fix index of new atoms
     #   print("DEBUG-new",atoms[loop_atoms[i]],atoms[-1],"DEBUG-new")

#    atoms[loop_atoms[2]][4][2] -= 1
#    atoms[loop_atoms[2]][4][0] -= 1
    print("    +doubled the atoms, shifted, and corrected, %d atoms now"%(len(atoms)))
    print("    +created dislocation loop with %d SIAs"%(Nloop))
    
    return 0 
        

def returnid(a,nx,ny,nz,nb):
    # returns the first id of the atom with ix = [ix,iy,iz,ib]
    # !!! THis function has been optimized and tested in another Python script!!!
    lx = nx[1]-nx[0]
    ly = ny[1]-ny[0]
    lz = nz[1]-nz[0]
    lb = nb[1]-nb[0]
    ix = a[0]; iy = a[1]; iz = a[2]; ib = a[3]

    if ix>=nx[1] or iy>=ny[1] or iz>=nz[1] or ib>=nb[1]:
        ia = -1
        print("error: invalid atim indices!")
        return ia
    # determine lbz
    tmp = iz-nz[0]
    lbz = tmp*lb*lx*ly
    
    # determine lby
    tmp = iy-ny[0]
    lby = lbz + tmp*lb*lx
    
    # determine lbx
    tmp = ix-nx[0]
    lbx = lby + tmp*lb
    
    ia = lbx+ib

    return ia
            

def wrap(i,a):
    # Taken and tested from 
    #  http://stackoverflow.com/questions/707370/clean-efficient-algorithm-for-wrapping-integers-in-c
    # Just need to decrease a[1] to account for the fact that the last index is the limit
    a[1] -=1
    r = a[1]-a[0]+1
    if i<a[0]:
        i += r*((a[0]-i)/r+1)
    tmp = a[0]+(i-a[0])%r
    a[1]+=1 # fix a because it is passed by reference
    return tmp

if __name__=="__main__":
    main()
