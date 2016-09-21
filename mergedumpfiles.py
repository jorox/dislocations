import sys

start = int(sys.argv[2])
end = int(sys.argv[3])
base = sys.argv[1]
fout = open(base+".full",'w')

for i in range(start,end+1):
    istimestep = False
    nlines = 0
    fin = open(base+"."+str(i))

    for line in fin:
        nlines +=1
        if istimestep:
            print(base+"."+str(i)+" --> TIMESTEP: "+line+"   ")
            istimestep = False
        if line == "ITEM: TIMESTEP\n":
            istimestep = True

        fout.write(line)
    
    print(nlines)

fout.close()
