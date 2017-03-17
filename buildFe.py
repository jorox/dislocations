import numpy as np
import sys

ver = "01"

class crystalbasis:
    motif = np.array([[0.0,0.0,0.0], [0.5,0.5,0.5], [0.5,-0.5,0.5],
                      [0.0,0.0,1.0], [0.5,0.5,1.5], [0.5,-0.5,1.5],[1.0,0.0,1.0],
                      [0.0,0.0,2.0],                [0.5,-0.5,2.5],[0.0,-1.0,2.0],
                      [-0.5,-0.5, 1.5], [1,0,2]])
    conv_motif = np.array([[0.0,0.0,0.0],[0.5,0.5,0.5]])

    def find_motif(self,X,Y,Z):
        e1 = np.array([1,0,0])
        e2 = np.array([0,1,0])
        e3 = np.array([0,0,1])
        Xnorm = np.linalg.norm(X)
        Ynorm = np.linalg.norm(Y)
        Znorm = np.linalg.norm(Z)
        # determine the coordinates of a1,2,3 in the XYZ basis
        rotmat = [ [X.dot(e1)/Xnorm, Y.dot(e1)/Ynorm, Z.dot(e1)/Znorm], 
                   [X.dot(e2)/Xnorm, Y.dot(e2)/Ynorm, Z.dot(e2)/Znorm], 
                   [X.dot(e3)/Xnorm, Y.dot(e3)/Ynorm, Z.dot(e3)/Znorm] ]
        rotmat = np.array(rotmat) #rotational matrix for the change of basis

        minh = int(min([X[0],Y[0],Z[0]]))-2
        mink = int(min([X[1],Y[1],Z[1]]))-2
        minl = int(min([X[2],Y[2],Z[2]]))-2
        maxh = int(max([X[0],Y[0],Z[0]]))+2
        maxk = int(max([X[1],Y[1],Z[1]]))+2
        maxl = int(max([X[2],Y[2],Z[2]]))+2
        self.motif = []
        new_motif = []
        self.motif.append([0,0,0])
        new_motif.append([0,0,0])
        print("########### DEBUG motif-find ############")
        print(minh,maxh,mink,maxk,minl,maxl)
        print(" ")
        for l in range(minl,maxl,1):
            for k in range(mink,maxk,1):
                for h in range(minh,maxh,1):
                    for convbasis in self.conv_motif:
                        newpnt = True
                        tmp1 = np.array([h,k,l])+convbasis
                        tmp2 = (rotmat[0]*tmp1[0]+
                               rotmat[1]*tmp1[1]+
                               rotmat[2]*tmp1[2])
                        tmp2 /= [Xnorm, Ynorm, Znorm]
                        tmp2 = [round(x,4) for x in tmp2]
                        for i in range(3):
                            if np.abs(tmp2[i]-1.0)<0.1: tmp2[i] = 0.0
                            
                        if (tmp2[0]<1 and tmp2[1]<1 and tmp2[2]<1 and
                        tmp2[0]>=0 and tmp2[1]>=0 and tmp2[2]>=0):
                            for imot in range(len(new_motif)):
                                if np.sqrt(np.sum((np.array(new_motif[imot])-tmp2)**2))<0.01:
                                    newpnt = False
                                    break
                            if newpnt:
                                print("%s <-- %s"%(tmp2,tmp1))
                                self.motif.append(list(tmp1))
                                new_motif.append(tmp2)
                            
        
        self.motif = np.array(self.motif)
        print("#########################################")
        
    def __init__(self,X,Y,Z):
        self.find_motif(X,Y,Z)
        e1 = np.array([1,0,0])
        e2 = np.array([0,1,0])
        e3 = np.array([0,0,1])
        Xnorm = np.linalg.norm(X)
        Ynorm = np.linalg.norm(Y)
        Znorm = np.linalg.norm(Z)
        # determine the coordinates of a1,2,3 in the XYZ basis
        self.a = [ [X.dot(e1)/Xnorm, Y.dot(e1)/Ynorm, Z.dot(e1)/Znorm], 
                   [X.dot(e2)/Xnorm, Y.dot(e2)/Ynorm, Z.dot(e2)/Znorm], 
                   [X.dot(e3)/Xnorm, Y.dot(e3)/Ynorm, Z.dot(e3)/Znorm] ]
        self.a = np.array(self.a) #rotational matrix for the change of basis

        print("\n         !!  Debug Info - a  !!\n!!----------------------------!!")
        print(self.a)
        print("!!----------------------------!!\n")
    

        # determine the coordinates of the motif atoms in the XYZ basis
        self.motif_XYZ = []
        for i in range(len(self.motif)):
            tmp = (self.a[0]*self.motif[i,0]+
                   self.a[1]*self.motif[i,1]+
                   self.a[2]*self.motif[i,2])
            tmp /= [Xnorm, Ynorm, Znorm]
            for j in range(3):
                if np.abs(tmp[j]-1.0)<0.11: tmp[j]=0
            self.motif_XYZ.append(tmp)
        print("\n         !!  Debug Info - motif_xyz  !!\n!!----------------------------!!")
        for tmp in self.motif_XYZ: print(tmp)
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
    print("\n"+banner("Fe-builder ver."+ver,length=100)+"\n")

    # define main variables for Zr
    # define common constants
    T = 10


    # default directions for the lab axes 
    X = np.array([1.,1.,1.])
    Z = np.array([1.,-1.,0.])
    Y = np.array([-1.,-1.,2.]) #default directions
    nx = [-1, 1]
    ny = [-1, 1]
    nz = [-1, 1] #default sizes
    createedge    = False
    createloop111 = False
    createscrew   = False
    
    load=["\\","|","/","-"]


    # input file for lattice construction
    fin = open(sys.argv[1])
    flog = open("build.log",'w')
    print("... reading input file " + sys.argv[1])
    flog.write("... reading input file " + sys.argv[1]+'\n')
    
# determine the direction of XYZ
    for line in fin:
        line = line.strip().split()
        if line[0]=="X" or line[0]=="x":
            print("     >>found new X direction ")
            flog.write("     >>found new X direction \n")
            X=mbvec(np.array(line[1:],dtype='float'))
        if line[0]=="Y" or line[0]=="y":
            print ("     >>found new Y direction ")
            Y=mbvec(np.array(line[1:],dtype='float'))
        if line[0]=="Z" or line[0]=="z":
            print("     >>found new Z direction ")
            Z=mbvec(np.array(line[1:],dtype='float'))
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
            
        if line[0]=="loop111":
            createloop111=True
            loop_center = int(line[1])
            n112 = [0, 1]
            n110 = [0, 1]
            print("      ++loop, centered "+str(loop_center)+", unit size")
        
        if line[0]=="n112":
            print("      ++loop n0001 direction changed")
            n112 = np.array(line[1:3],dtype='int')
        if line[0]=="n110":
            print("      ++loop n1100 direction changed")
            n110 = np.array(line[1:3],dtype='int')
        
        if line[0]=="Temp" or line[0]=="temp" :
            print("      ++found temperature")
            T = int(line[1])
    
    if T==10:
        a = 2.856
    if T == 300:
        a = 2.856

            
# !Need to check that directions are orthogonal!
# !Need to check that limits do not coincide!        

# determine the crystal orientation
    print("\n... orienting crystal");flog.write("... orienting crystal\n")
    print("     X/Y/Z directions: "); flog.write("     X/Y/Z directions: \n")
    print("X = %s"%(str(X)))
    print("Y = %s"%(str(Y)))
    print("Z = %s"%(str(Z)))
    basis = crystalbasis(X,Y,Z)

# determine the spacings and number of box cells
    print("\n... setting up simulation box")
    xspacing = np.linalg.norm(X)*a;
    yspacing = np.linalg.norm(Y)*a;
    zspacing = np.linalg.norm(Z)*a;

    print("     spacings (A): %1.5f, %1.5f, %1.5f"%(xspacing,yspacing,zspacing))
    flog.write("     spacings (A): %1.5f, %1.5f, %1.5f\n"%(xspacing,yspacing,zspacing))
    print("     number of unit cells in XYZ directions: %1.0d %1.0d %1.0d"
          %(np.diff(nx), np.diff(ny), np.diff(nz)))
    flog.write("     number of unit cells in XYZ directions: %1.0d %1.0d %1.0d\n"
          %(np.diff(nx), np.diff(ny), np.diff(nz)))
    
# determine the box limits
    box = [[nx[0]*xspacing,nx[1]*xspacing],
           [ny[0]*yspacing,ny[1]*yspacing],
           [nz[0]*zspacing,nz[1]*zspacing]]
    box = np.array(box)
    print("\n... box dimensions (A):"); flog.write("... box dimensions (A):")
    print("     X %1.5f --> %1.5f  %3d -> %3d = %1.5fnm\n     Y %1.5f --> %1.5f    %3d -> %3d = %1.4fnm\n     Z %1.5f --> %1.5f   %3d -> %3d = %1.5fnm" %(box[0,0],box[0,1],nx[0],nx[1],np.diff(box[0]/10),
                box[1,0],box[1,1],ny[0],ny[1],np.diff(box[1]/10),
                box[2,0],box[2,1],nz[0],nz[1],np.diff(box[2])/10 ))
    flog.write("\n     X %1.5f --> %1.5f \n     Y %1.5f --> %1.5f \n     Z %1.5f --> %1.5f\n"
               %(box[0,0],box[0,1],box[1,0],box[1,1],box[2,0],box[2,1]))

# create atoms
    natoms = (nx[1]-nx[0])*(ny[1]-ny[0])*(nz[1]-nz[0])*len(basis.motif) #4 basis atoms
    nghost = 0

    print("... expecting %1.0d atoms"%natoms)
    atoms = []
    ia = -1 # atoms id (0-based)
    for iz in range(nz[0],nz[1]):
        for iy in range(ny[0],ny[1]):
            for ix in range(nx[0],nx[1]):
                for ib in range(len(basis.motif)): # 2 atoms per basis poin
                    tmp = np.array([ix,iy,iz])+basis.get_xyz_motif(ib)    # transform to Miller notation
                    #flog.write("\n\n")
                    #flog.write(str(tmp)+"\n") 
                    # Now add the motif to the box point. 
                    # NOTE: one could change the motif coordinates to a 'partial' 
                    #  Miller-Bravais vector. However, this using Miller notation makes 
                    #  reading the coordinates of the atoms easier to read
                    # IMPORTANT: the motif is specific ONLY to the default directions 
                    #  specified. For a different direction CHANGE THE MOTIF!!
                    
                    tmp[0] *= xspacing
                    tmp[1] *= yspacing
                    tmp[2] *= zspacing
                    ia += 1
                    atoms.append([ia+1, tmp[0],tmp[1],tmp[2],[ix,iy,iz,ib]])
                    prcnt = float(ia)/natoms*100
                    print(load[ia%4]+ "     %d %%"%(prcnt))
                    sys.stdout.write("\033[F")

    print("... %d atoms, %d ghost atoms"%(natoms,nghost))
    
# create the edge dislocation
    if createedge:
        print("\n /////////////// Building 111 Edge dislocation \\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        nghost = buildedge111(atoms,nx,ny,nz,[xspacing,yspacing,zspacing],box)
        natoms -= nghost
        print("... new box dimensions along X: %1.4f %1.4f"%(box[0,0],box[0,1]))
        print("... %1.0f atoms, %1.0f ghost atoms"%(natoms,nghost))
        print("/\\\\\\\\\\/////////\\\\\\\\\\//////////\\\\\\\\\\/////////////")

# create the screw dislocation
    if createscrew:
        print("\n /////////////// Building Screw dislocation \\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        box = [[nx[0]*xspacing,nx[1]*xspacing],
               [ny[0]*yspacing,ny[1]*yspacing],
               [nz[0]*zspacing,nz[1]*zspacing]]
        box = np.array(box)
        buildscrew(atoms,nx,ny,nz,xspacing/2.0,box)
        fixids(atoms)
        
        print("/\\\\\\\\\\/////////\\\\\\\\\\//////////\\\\\\\\\\/////////////")
    
# create the loop 
    if createloop111:
        print("\n/////////////// Buidling SIA 111 Loop \\\\\\\\\\\\\\\\\\\\\\\\\\\\")
        nloop = buildloop_111(atoms,basis,a,nx,ny,nz,n112,n110,loop_center)
        natoms += nloop
        print("/\\\\\\\\\\/////////\\\\\\\\\\//////////\\\\\\\\\\/////////////")

# write the XYZ file
    foutname = "Fe.xyz"
    fout = open(foutname, 'w')
    fout.write("%d atoms"%(len(atoms)))
    for i in range(len(atoms)):
        if atoms[i][0] == -1: continue #ghost atom
        fout.write("\nZr %1.6f %1.6f %1.6f"%(
                             atoms[i][1], atoms[i][2],atoms[i][3]))
    fout.close()
    correct = True

# write LAMMPS input file
    flammps = open("fe.data","w")
    flammps.write("#ver. "+ver+" BCC-Fe [111],[-1-12],[1-10] + dislocation and loop, metal units, "+str(T)+"K")
    flammps.write("\n    %d atoms\n    %d atom types\n"%(natoms,1)) #N atoms, 1 type
    if correct and createscrew:
        print("--> warning: shearing box for screw dislocation")
        Ly = box[1,1]-box[1,0]
        yhi = np.sqrt(Ly**2-xspacing**2/4/4)+box[1,0]
        flammps.write("\n%1.6f %1.6f xlo xhi\n%1.6f %1.6f ylo yhi\n%1.6f %1.6f zlo zhi"%
                  (box[0,0],box[0,1], box[1,0],yhi, box[2,0],box[2,1]) )
        flammps.write("\n%1.4f %1.4f %1.4f xy xz yz\n"%(-xspacing/4,0.0,0.0))
    else:
        flammps.write("\n%1.6f %1.6f xlo xhi\n%1.6f %1.6f ylo yhi\n%1.6f %1.6f zlo zhi\n"%
                  (box[0,0],box[0,1], box[1,0],box[1,1], box[2,0],box[2,1]) )
    flammps.write("\nMasses\n")
    flammps.write("\n    1 55.845\n")
    flammps.write("\nAtoms\n")
    ia = 0
    for i in range(len(atoms)):
       # print("DEBUG-atoms",i, atoms[i][4], "DEBUG-atoms")
        if atoms[i][0] == -1: continue #ghost atom
        ia +=1
        flammps.write("\n\t%d\t%d\t%1.6f %1.6f %1.6f"%(
                ia,1, atoms[i][1], atoms[i][2], atoms[i][3]))
    flammps.close()


    print("\n... done writing %1.0d atoms to file %s and lammps file"%(natoms,foutname))    
    flog.write("\n... done writing %1.0d atoms to file %s\n"%(natoms,foutname))
    
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

    fsnap = open("fe.ov"+ver,"w")
    fsnap.write("ITEM: TIMESTEP\n0")
    fsnap.write("\nITEM: NUMBER OF ATOMS\n%d"%(natoms))
    if correct and createscrew:
        print("--> warning: shearing box for screw dislocation")
        fsnap.write("\nITEM: BOX BOUNDS xy xz yz pp pp ss\n%1.5f %1.5f %1.5f\n%1.5f %1.5f %1.5f\n%1.5f %1.5f %1.5f"%
                      (box[0,0],box[0,1],-xspacing/4, box[1,0],yhi,0.0, box[2,0],box[2,1],0.0) )
    else:
        fsnap.write("\nITEM: BOX BOUNDS pp pp ss\n%1.5f %1.5f\n%1.5f %1.5f\n%1.5f %1.5f"%
                      (box[0,0],box[0,1], box[1,0],box[1,1], box[2,0],box[2,1]) )
    fsnap.write("\nITEM: ATOMS id type x y z")
    ia = 0
    for i in range(len(atoms)):
        if atoms[i][0] == -1: continue #ghost atom
        ia += 1
        fsnap.write("\n%d 1 %1.4f %1.4f %1.4f"%(ia,atoms[i][1],atoms[i][2],atoms[i][3]))
    fsnap.close()
    if ia != natoms:
        print("WARNING: incosistent number of atoms in OVITO file")
    print("... done writing to snapshot file fe.ov"+ver)

#################################################################################################

def buildedge111(atoms,nx,ny,nz,spacing,box):
    zcut = (nz[1]+nz[0])/2.0
    N = 2*(nx[1]-nx[0])
    xcut = nx[1]-1
    xmin = nx[0]*spacing[0]
    nghost = 0
    exx1 = -0.5/N
    exx2 = 0.5/(N-1)

    for ia in range(len(atoms)):
        iz=atoms[ia][4][2]; iy = atoms[ia][4][1]; ix=atoms[ia][4][0]; ib = atoms[ia][4][3]
        xunit =(atoms[ia][1]-ix*spacing[0])/spacing[0] #the relative position of the atom in its unit cell
        zunit =(atoms[ia][3]-iz*spacing[2])/spacing[2]
        
        if iz>zcut or (iz==zcut and zunit>0.45) :
            atoms[ia][1] += (atoms[ia][1]-xmin)*exx1
        else:
            if ix==xcut and xunit>0.49:
                atoms[ia][0] = -1
                nghost +=1
                continue
            else:
                atoms[ia][1] += (atoms[ia][1]-xmin)*exx2
        

    box[0][1] -= 0.25*spacing[0]
    return nghost


def buildscrew(atoms,nx,ny,nz,a,box):
    # The function creates a screw [11-20] dislocation in the basal plane
    # The method used is that taken from Serra Bacon 2013 and Khater Bacon 2010
    # The method follows the technique in Bacon "Intro to Dislocations" eq. (4.11) with z->x, y->z, x->y
    # The displacement field is given by 
    #    ux = (b/2/pi) (arctan(z/y)-pi) = n/2/pi (theta-pi) theta = 0 to 2pi
    # For y=0, theta can be either pi/2 or 3pi/2 or 0,
    #  Intro to Dislocations p.66 that for x=0 check the drawing
    cntr_pnt = np.sum(box,1)/2.0
    print cntr_pnt
    fout = open("dirty.dat","w")
    fout.write("#z y theta")
    Lx = np.diff(box[0],1)[0]
    for i in range(len(atoms)):
        y = atoms[i][2]-cntr_pnt[1]
        z = atoms[i][3]-cntr_pnt[2]
        #if y == 0: # y=0 
        #    if z > 0 :  ux = a/4.0 #theta = pi/2 
        #    elif z < 0: ux = -a/4.0 #theta=3pi/2 
        #    else: ux = 0 #core (theta = pi) 
        #elif z == 0:
        #    if y > 0: ux = a # theta = 2pi
        #    elif y < 0: ux = a/2.0 #theta = pi
        #    else: ux = a/2.0
        #else:
        ux = a/2/np.pi*(myatan(z,y)) #default
        fout.write("\n%1.0f %1.3f %1.3f %1.3f"%(i,z,y,ux/a))
        atoms[i][1] += ux
        atoms[i][1] = wrap_box(atoms[i][1],box[0])
        #if atoms[i][1] < box[0,0] or atoms[i][1]>= box[0,1]: print(
        #        "ERROR_buildscrew: atom %i %1.4f"%(i,atoms[i][1]))

    # To restablish periodicity across the Y-direction a displacement of +/-b/2 
    # was added to the x-coordinate of atoms on the +/- y boundaries
    #for ix in range(nx[0],nx[1]):
    #    for iz in range(nz[0],nz[1]):
    #        for ib in range(4): # 4-atoms per box cell
    #            aid = returnid([ix,ny[0],iz,ib],nx,ny,nz, [0,4])
    #            atoms[aid][1] -= a/2 # -y boundary
    #            aid = returnid([ix,ny[1]-1,iz,ib], nx, ny, nz, [0,4])
    #            atoms[aid][1] += a/2 # +y boundary

def fixids(atoms):
    for i in range(len(atoms)):
        atoms[i][0] = i+1
    return 0

def buildloop_111(atoms,basis,a,nx,ny,nz,n112,n110,xloop):
    """
    Adds interstitial atoms to create a loop with Burgers vector [111]
    """
    loop_atoms = [] #holds ids of loop atoms = 4xseed_list
    catoms = [] #holds ix,iy,iz,ib of loop atoms = 4xseed_list
    seed_list = get_seeds_111(n112,n110,xloop)
    print("DEBUG seed_list",seed_list,"#DEBUG")
    iseed = -1
    for seed in seed_list:
        iseed +=1
        #if seed[3]%2 == 0: print("ERROR: Wrong seed for loop60 constrcution") 
        catoms.append(list(seed)) #append seed atom
        for itemp in range(len(basis.motif)-1):
            seed[3] += 1
            catoms.append(list(seed)) #append all the basis atoms in a cell (we are assuming they are 12)
            
    # Now get the ids if the loop atoms
    for ic in catoms:
        tmp = returnid(ic,nx,ny,nz,[0,len(basis.motif)])
        if atoms[tmp][1]-a*np.sqrt(3)*atoms[tmp][4][0] < 0.5*a*np.sqrt(3): #ensure atom is in the first half of the unit cell along X -- we are assuming that X = [111]
            loop_atoms.append(tmp)
            #print("DEBUG-loop \n",atoms[loop_atoms[-1]][4],"\n DEBUG-loop")
            if atoms[loop_atoms[-1]][0] == -1:
                print("Warning: trying to modify ghost atom")
            
    
    shift = [0.05*a*0.5, 0.0, 0.0]
    add_sias(atoms,loop_atoms,shift)
    #print("DEBUG",loop_atoms,"DEBUG")
    return len(loop_atoms)

def get_seeds_111(n110,n112,xloop):
    """
    Get a list of seed atoms for the construction of a SIA loop with Burgers = [111].
    n110 and n112 specify the max-min stacking indices along those two directions
    n110/n112 is a vector of the kind [i,j] where i and j are the indices of the unit cell in the crystal
    For example:
          0      1  ...
     --|---- -|------------
       |      |
 1     |  X   |   X
     --|----- |------------
       |  X   |   X
 0     |      |   
    - -|----- |------------
    n110 = [0 1] and n112= [0 1] 
    xloop is the location of the top most segment
    """
    seed_list = []
    for i112 in range(n112[1],n112[0]-1,-1):
        seed = [xloop,i112,n110[1],0] # THIS IS THE FIRST BASIS ATOM
        seed_list.append(list(seed))  # it is the list of seed atoms to construct the loop in each unit cell
        for i110 in range(n110[1]-1,n110[0]-1,-1):
            seed = [xloop,i112,i110,0]
            seed_list.append(list(seed))
            
    return seed_list

def add_sias(atoms, loop_atoms, shift):
    """
    Add Self-interstitial atoms to the list specified in loop_atoms.
    The shift is a vector in XYZ space defined in the variable shift
    """
    last_index = atoms[-1][0]
    Nloop = len(loop_atoms)
    print("    +found %d atoms for the loop"%(Nloop))
    print("DEBUG-shift",shift,"DEBUG-shift")
    
    for i in range(Nloop):
        last_index +=1
        oldatom = atoms[loop_atoms[i]]
        print("DEBUG-old",oldatom,"DEBUG-old")
        atoms.append([oldatom[0],oldatom[1],oldatom[2],oldatom[3],
                      [oldatom[4][0],oldatom[4][1],oldatom[4][2],oldatom[4][3]]])  #store old position
    
        for ix in range(1,4):
            atoms[loop_atoms[i]][ix] +=  shift[ix-1]
            atoms[-1][ix] -= shift[ix-1]
        atoms[-1][0] = last_index # fix index of new atoms
        print("DEBUG-new",atoms[loop_atoms[i]],atoms[-1],"DEBUG-new")

#    atoms[loop_atoms[2]][4][2] -= 1
#    atoms[loop_atoms[2]][4][0] -= 1
    print("    +doubled the atoms, shifted, and corrected, %d atoms now"%(len(atoms)))
    print("    +created dislocation loop with %d SIAs"%(Nloop))
    
    return Nloop 

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
    a[1]+=1 # fix a, passed by reference
    return tmp

def wrap_box(x,a):
    # Wrap a float x into the range specified by a such that a[0]<=x_wrapped<a[1]
    # tested and verified with a simple example:
    # >>> for x in np.arange(-4,5,0.5):
    #...     print x, wrap_box(x,[-2.5,3])
    # -4.0 1.5
    # -3.5 2.0
    # -3.0 2.5
    # -2.5 -2.5
    # -2.0 -2.0
    # -1.5 -1.5
    # -1.0 -1.0
    # -0.5 -0.5
    # 0.0 0.0
    # 0.5 0.5
    # 1.0 1.0
    # 1.5 1.5
    # 2.0 2.0
    # 2.5 2.5
    # 3.0 -2.5
    # 3.5 -2.0
    # 4.0 -1.5
    # 4.5 -1.0

    if x<a[0]: x = x%a[0] + a[1]
    elif x>=a[1]: x = x%a[1] + a[0]
    if x >= a[1] or x<a[0]: print("ERROR: wrap_box x=%1.2f"%(x))
    return x

def myatan(x,y):
    PI = np.pi
    if x>0 and y>=0: #gives 0
        return np.arctan(y/x)
    elif x<0 and y>=0: #gives PI
        return PI-np.arctan(-y/x)
    elif x<0 and y<0:
        return PI+np.arctan(y/x)
    elif x>0 and y<0:
        return 2*PI-np.arctan(-y/x)
    else: #x==0:
        if y>0: return PI/2
        elif y<0: return 3*PI/2
        else: return 0
def get_loop_size(loop_atoms):
    """ returns an array of 3 elements size of loop in x,y,z
    """

if __name__=="__main__":
    main()
