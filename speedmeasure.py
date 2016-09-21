# This python measures the speed of a dislocation segment from a LAMMPS dump file
# The LAMMPS dump file must contain all of just the atoms in the dislocations
# segment
# The header for the file should look like:
#   id x y z c_sym

import sys
import numpy as np

def main():
    print("... opening "+sys.argv[1]+" to read")
    fin_details = get_file_info(sys.argv[1])
    print("    %d timesteps, %d atoms"%(fin_details[0],fin_details[1]))
    print("    "+str(fin_details[2]))
    scaled = False
    for i in range(len(fin_details[2])):
        if fin_details[2][i][0] == "c_csym":
            i_csym = fin_details[2][i][1]
        if fin_details[2][i][0] == "x": ix = fin_details[2][i][1]
        if fin_details[2][i][0] == "xs": ix =fin_details[2][i][1]; scaled = True
        if fin_details[2][i][0] == "y": iy =fin_details[2][i][1]
        if fin_details[2][i][0] == "ys": iy =fin_details[2][i][1]; scaled = True
        if fin_details[2][i][0] == "z": iz =fin_details[2][i][1]
        if fin_details[2][i][0] == "zs": iz =fin_details[2][i][1]; scaled = True

    #########################################################################
    # Formulas for navigation of a LAMMPS file with constant number of atoms
    # in line-numbers starting from 0 and python countings
    #   timestep   = x+1
    #   # of atoms = x+3
    #   box        = [x+5,x+8]
    #   1st atom  = x+9
    #==> nth atom  = x+9+(n-1)
    #==> final atom= x+9+(N-1)
    # ---------------------------
    # ==> Total number of lines per ssection = 9+N_atoms
    #
    # First line of 
    #  1st step : x=0
    #  2nd step : x=9+N_atoms
    # ==>  nth step : x = (9+N_atoms)*(n-1)
    #
    #==> Total lines per file = (9+N_atoms)
    print("... calculating speed of dislocation")
    nl = -1
    nstep = 1
    fin = open(sys.argv[1])
    box = [[0,0],[0,0],[0,0]]
    com_trajectory = []
    
    # loading bar
    load_bar = ["\\","|","/","-"]
    ndis = 0 # number of atoms in dislocation
    for line in fin:
        nl +=1
        if nl==0: continue
        x = (9+fin_details[1])*(nstep-1)
        line = line.strip().split()
        s = nl-x # remove the x portion
        if s<5: continue

        if s==5: # box section
            box[0] = [float(line[0]),float(line[1])]
        if s==6: # box section
            box[1] = [float(line[0]),float(line[1])]
        if s==7: # box section
            box[2] = [float(line[0]),float(line[1])]
            com = np.array([0.0,0.0,0.0])
            lz = box[2][1]-box[2][0]
            lz_inner = lz/2-12 #this assumes the box bounds in the z-direction are -z and +z
            

        if s>=9: #atoms section
            # check if the c_csym is not 2 and is not boundary
            if line[i_csym] != "2":
                if scaled: tmp =float(line[iz])*lz+box[2][0] #unscale 
                else: tmp = float(line[iz])
                if np.abs(tmp)<lz_inner: # +/-z boundaries!!
                    ndis += 1
                    tmp = np.array([float(line[ix]),
                                float(line[iy]),
                                float(line[iz])])
                    com = com + tmp

            # rerached the last atom
            if s==9+fin_details[1]-1: 
                com = com/ndis
                if scaled:
                    for j in range(3):
                        com[j] = com[j]*(box[j][1]-box[j][0])+box[j][0]
                com_trajectory.append(com)
                ndis = 0
                nstep += 1
        
        # progress report
        prcnt = float(nstep)/fin_details[0]*100
        print(load_bar[nstep%4]+"    %d %%"%(prcnt))
        sys.stdout.write("\033[F")
        ##############################################
    
    outnm = sys.argv[2]
    fout = open(outnm,"w")
    fout.write("#Trajectory of dislocation from file "+sys.argv[1])
    fout.write("\n#step x(A) y(A) z(A)") 
    for i in range(len(com_trajectory)):
        fout.write("\n%d %1.5f %1.5f %1.5f"%(i,com_trajectory[i][0],
                                           com_trajectory[i][1],
                                           com_trajectory[i][2]))
    print("... done writing to "+outnm)
        

    




def get_file_info(fnm):
    fin = open(fnm)
    # Rerturns the number of atoms and timesteps in a file
    nStep = 0
    nATOMS = 0
    next_line_is_natoms = False
    nLINES = 0
    header = []

    for line in fin:
        line = line.strip()

        if next_line_is_natoms: #natoms 
            nATOMS = int(line)
            next_line_is_natoms = False 
            continue

        if line=="ITEM: NUMBER OF ATOMS": next_line_is_natoms=True # natoms
        if "ITEM: ATOMS" in line: # header line
            line = line.split()
            for i in range(2,len(line)): 
                    header.append((line[i],i-2))
            break

    fin.close()
    
    fin = open(fnm)
    for line in fin:
        nLINES +=1
    fin.close()
    
    if nLINES%(9+nATOMS)!=0:
        print("ERROR ATOMS COUNT NOT CONSTANT")
    else:
        nStep = nLINES/(9+nATOMS)

    return (nStep,nATOMS,header)


if __name__=="__main__":
    main()
