import sys
import os
import glob
import numpy as np
from argparse import ArgumentParser

class point:
    def __init__(self,tag,dof,x,v):
        self.x = [s for s in x]
        self.v = [s for s in v]
        self.tag = tag
        self.dof = dof
        self.vld = True

    def isPhysical(self):
        if self.dof == 1:
            return True
        else:
            return False
    def isValid(self):
        return self.vld
    def notValid(self):
        self.vld = False
    def sep(self,other,w):
        return np.sum( (np.array(self.x)-np.array(other.x))**2*np.array(w) )**0.5
    def __lt__(self,other):
        return NotImplemented
    def __gt__(self,other):
        return NotImplemented
    def __le__(self,other):
        return NotImplemented
    def __ge__(self,other):
        return NotImplemented
    def __eq__(self,other):
        return self.sep(other)<1e-5
    def __ne__(self,other):
        return self.sep(other)>1e-5
    def __str__(self):
        return "%i, %i, [%1.2f,%1.2f,%1.2f], [%1.2f,%1.2f,%1.2f]"%(self.tag,self.dof,
                                                                       self.x[0],self.x[1],self.x[2],
                                                                       self.v[0],self.v[1],self.v[2])


class pnodelist:
    def __init__(self):
        self.pnts = []
        self.current = 0

    def append(self,newpnt):
        self.pnts.append(newpnt)
        
    def __len__(self):
        return len(self.pnts)
    def __str__(self):
        s = ""
        for i in range(len(self.pnts)):
            s += "%s\n"%(self.pnts[i])
        s = s[:len(s)-1]
        return s
    def __iter__(self):
        self.current = -1
        return self
    
    def __next__(self):
        if self.current < len(self.pnts)-1:
            self.current+=1
            return self.pnts[self.current]
        else:
            raise StopIteration

    def __getitem__(self,i):
        return self.pnts[i]

    def sep(self,other,wghts):
        if len(self)!=len(other):
            N = min(len(self),len(other))
            print("WARNING: trying to calculate difference between unequal lists. Using minimum length")
        else: N = len(self)
        tmp = [-1]*N
        for i in range(N): tmp[i] = self.pnts[i].sep(other.pnts[i],wghts)
        return tmp

    def sort(self, other):
        # calculate distance from every point in other
        if len(self)==len(other):
            old = list(self.pnts)
            N = len(self)
            keys = range(N)
            for j in range(N-1): #last element does not need sorting
                sep_list = []
                for i in range(j,N):
                    #find separation between other[j] and points self[key[i]]
                    sep_list.append([other[j].sep(self[keys[i]]),i])
                sep_list = sorted(sep_list, key=lambda s: s[0])
                swap_i = sep_list[0][1]
                tmp = keys[swap_i]
                keys[swap_i] = keys[j]
                keys[j] = tmp
            old = list(self.pnts)
            for i in range(len(keys)):
                self.pnts[i] = old[keys[i]]
            if np.sum(np.abs(self.sep(other))) > 20:
                print("WARNING: sorting failed")
        else:
            print("WARNING: sorting unequal lists")

    def better_sort(self,other,wghts,debug=None):
        """
        Sort self according to other even if number of elements in self and other are different
        """
        if len(self)==1 and len(other)==1: return #1 physical node -- unchanged -- do nothing
        if debug is None: debug = False
        #calculate diff matrix
        old = list(self.pnts) #old points
        R_diff = [] #matrix of distances between all nodes (all permutations, unsorted)
        split_event = False
        merge_event = False
        if len(self)>len(other):
            split_event = True
            N_new_nodes = len(self)-len(other)
        if len(self)<len(other):
            merge_event = True
            N_merged_nodes = len(other)-len(self)
            
        for snode in self.pnts:
            R_diff.append([])
            for onode in other.pnts:
                R_diff[-1].append(snode.sep(onode,wghts))
            
        R_diff = np.array(R_diff)
        #return the minimum along the columns
        pos = np.argmin(R_diff,axis=1) #find min for each row=sorted position, can have repetitions and omitions

        if split_event:
            #find "N_new_nodes" nodes with greatest distance from all "other" nodes
            print("+++ better_sort(): Split event: %i new nodes"%(N_new_nodes))
            min_sep = np.amin(R_diff,axis=1) #min separation each node
            max_ind = np.argsort(min_sep)[-N_new_nodes:] #indices of new nodes with largest separation
            
            #set new nodes to final positions
            for i in range(N_new_nodes):
                pos[max_ind[i]] = len(other)+i
        elif merge_event:
            ### In case of a merge event invalidate point and keep last position
            # Find unused positions
            tmp = 0
            print("+++ better_sort(): Merge event: %i nodes lost"%(N_merged_nodes))
            for i in range(len(other)):
                if i not in pos:
                    tmp +=1
                    self.append(point(other[i].tag,other[i].dof,other[i].x,other[i].v))
                    pos = np.append(pos,i)
            if tmp != N_merged_nodes:
                print(">>> ERROR: better_sort() unhandeled event during merge event"); sys.exit()
        if debug:
            print("################## DEBUG ####################")
            print(pos)
            for i in range(len(other)): print(other.pnts[i])
            print("#-------------NEW--------------")
            for i in range(len(self)): print(self.pnts[i])
            print("#------------R-diff------------")
            for i in range(len(R_diff)): print(R_diff[i].tolist())
            print("#############################################")
            
        #Now sort
        if np.sum(pos) != np.sum(range(len(self))):
            print(">>> ERROR: sort_better(): Position-List not unique"); sys.exit()
        for i in range(len(old)):
            self.pnts[pos[i]] = old[i]

            
#### Main program
#---------------------------------------------------------------------------------------
        
pers = ArgumentParser()

pers.add_argument("fin", help="input file name pattern containing *", metavar="file-pattern")
pers.add_argument("-o","--out", help="name of output filename", default="get_v_vtk.dat")
pers.add_argument("-n","--nfiles", help="number of files to use",type=int,default=-1)
pers.add_argument("-s","--sigeps",help="sigeps file", default="SIGEPS")
pers.add_argument("-w","--weights", help="modify the contribution of xyz when calculating separation", nargs=3,
                  type = float, default=(1.,1.,1.))
pers.add_argument("-d","--debug", help="output information helpful for debugging",
                  default=False, action="store_true")
pers.add_argument("-c","--com", help="Do not sort the list of nodes, just output the center of mass coordinates",
                  action="store_true", default=False)

args = pers.parse_args()

print("... opening file: %s"%(args.sigeps))
fsigeps = open(args.sigeps)
time = [0]
for line in fsigeps:
    if line[0].isalpha(): continue
    time.append(float(line.strip().split()[1]))
print("      +++found %i time-step"%(len(time)))
    
iwild = args.fin.find("*")
if iwild == -1:
    print("ERROR: no wildcard found in file pattern")
    exit

flist = glob.glob(args.fin) #list of filenames
print("... %1.0d files found"%(len(flist)))
masterlist = [] #list to hold lists of physical nodes
vout = open(args.out+".vel", "w")
#tagout = open(args.out+".tag","w")
posout = open(args.out+".pos","w")
step = -1
load = ["\\","|","/","-"]

flist = sorted(flist,key=lambda s: int(s[s.rfind("_")+1:s.rfind(".")]))
flist = flist[:args.nfiles]
print("... Using %i files"%(len(flist)))
#############DEBUGGING###################
sepout = open("listsep.dat","w")
sepout.write("#iterative separation between p-node lists")
sepout.write("\n#step s1 s2 s3 s4")
#########################################
for ff in flist:
    step += 1
    fvtk = open(ff)
    line = fvtk.readline()
    load = ["\\","|","/","-"]
    while line:
    	line = line.strip()
        if line == "":
            line = fvtk.readline()
            continue

        #save positions of points
        if line.split()[0] == "POINTS":
            npoints = int(line.split()[1])
            pos = []
            for i in range(npoints):
                line = fvtk.readline().strip().split()
                pos.append([float(s) for s in line])

        #save tags
        if line == "SCALARS tags int":
            line = fvtk.readline()
            line = fvtk.readline().strip().split()
            tags = [int(s) for s in line]
            
        #save DOF    
        if line == "SCALARS dof int":
	   line = fvtk.readline()
	   line = fvtk.readline().strip().split()
	   dofs = [int(s) for s in line]

        #save velocities
        if line == "VECTORS NodalVelocities float":
            vel = []
            for i in range(npoints):
                line = fvtk.readline().strip().split()
                vel.append([float(s) for s in line])

        line = fvtk.readline()

    # create list for a file
    nlist = pnodelist()
    if args.com == False:
        for i in range(npoints):
            if dofs[i] != 1: continue #avoid non-physical nodes
            if tags[i] in tags[:i]: continue #avoid repetitions by checking tags!
            nlist.append(point(tags[i],dofs[i],pos[i],vel[i]))
    else:
        comx= np.array([0.0,0.0,0.0])
        comv = np.array([0.0,0.0,0.0])
        for i in range(npoints):
            comx += np.array(pos[i])
            comv += np.array(vel[i])
        comx /= npoints
        comv /= npoints
        nlist.append(point(-1,1,comx,comv))
        
    print("step %i: %i nodes"%(step,len(nlist)))
    
    if step>0 and args.com==False:
        nlist.better_sort(masterlist[step-1],args.weights,args.debug)
        #############DEBUGGING###################
        tmp = masterlist[step-1].sep(nlist,args.weights)
        if np.sum(np.abs(tmp))>20*np.sqrt(np.sum(args.weights)):
            print(">>> WARNING: Possible invalid separation after sorting list step %i"%(step))
        sepout.write("\n%1.2d %1.4f %1.4f %1.4f %1.4f"%(step,tmp[0],tmp[1],tmp[2], tmp[3]))
        #########################################
    masterlist.append(nlist)
    
    if step == 0:
        print(nlist)
        print("... initial physical nodes %i"%(len(nlist)))


    
#### Output
#-----------------------------------------------------------------------
maxpnodes = 0
biggest_list = -1
for i in range(len(masterlist)):
    if len(masterlist[i]) > maxpnodes:
        maxpnodes = len(masterlist[i])
        biggest_list = i 
print("... max p-nodes %i in file %i"%(maxpnodes,biggest_list))
#print("   %s"%(masterlist[biggest_list]))
print("... creating output file with %i columns"%(3*maxpnodes))

vout.write("#velocity date from vtk files for %i nodes"%(maxpnodes))
posout.write("#velocity date from vtk files for %i nodes"%(maxpnodes))
vout.write("\n#step time(ns)")
posout.write("\n#step time(ns)")
tmp = ["x","y","z"] #i/maxpnodes = 0="x" or 1="y" or 2="z"
for i in range(3*maxpnodes):
    vout.write(" v%s_%i($%i)"%(tmp[i/maxpnodes],i%maxpnodes+1,i+3)) #x.y.z, 1-->mpn, 2-->3*mpn+1 
    posout.write(" %s_%i($%i)"%(tmp[i/maxpnodes],i%maxpnodes+1,i+3))
    
step = 0
dout = open("dout.dat","w")
dout.write("#debug file, indices")
i_t = -1
for ilist in masterlist:
    i_t +=1
    vout.write("\n%1.2d %1.4f"%(step, time[i_t]))
    posout.write("\n%1.2d %1.4f"%(step, time[i_t]))
    dout.write("\n%1.2d"%(step))
    step += 1
    N = len(ilist) #number of nodes in list
    for i in range(3*maxpnodes):
        dout.write(" (%i,%i)"%(i%maxpnodes,i/maxpnodes))
        if i%maxpnodes < N: # index exists
            vout.write(" %1.4f"%(ilist[i%maxpnodes].v[i/maxpnodes]))
            posout.write(" %1.4f"%(ilist[i%maxpnodes].x[i/maxpnodes]))
        else: #index doesnot exist (yet)
            vout.write(" NaN")
            posout.write(" NaN")

if args.com: print("... calculated COM data")
print("... done writing to "+vout.name)

            
    
                
	   

