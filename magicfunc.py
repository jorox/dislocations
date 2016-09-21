import sys
import os

def main():
    if os.path.isfile("dummy.txt"):
        print "... dummy.txt file already exists"
    else:
        atoms = []
        ia = -1
        for iz in range(-3,3):
            for iy in range(-3,3):
                for ix in range(-5,5):
                    for ib in range(4):
                        ia +=1
                        atoms.append([ia,[ix,iy,iz,ib]])

        fout = open("dummy.txt","w")
        fout.write("#dummy file to test the magic func")
        for i in range(len(atoms)):
            fout.write("\n %d | %d %d %d %d"%(
                    atoms[i][0], 
                    atoms[i][1][0], atoms[i][1][1], atoms[i][1][2],atoms[i][1][3]))
        fout.close()
        print "... done writing file dummy.txt"
    
    tmp = [int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4])]

    print("... id = "+ str(getid(tmp,[-5,5],[-3,3],[-3,3],[0,4])) )

    
def getid(a,nx,ny,nz,nb):
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
    
if __name__=="__main__":
    main()
