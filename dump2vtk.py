from argparse import ArgumentParser

pers = ArgumentParser()

pers.add_argument("fin", help="Input file", metavar="File-IN motif")
pers.add_argument("-o", "--output", metavar="File-OUT motif", default="batchdxa")

args = pers.parse_args()
print(args)
d = dump(args.fin)
t = d.time()

for it in time:
    natoms = len(d.vecs(it,"id"))
    tmp = [1]*natoms
    print("... %i step --> %i atoms"%(it,natoms))
    d.tselect.one(it)
    
    d.setv("type",tmp)
    v = vtk(d)
    v.one(args.output+"."+str(it))
    

