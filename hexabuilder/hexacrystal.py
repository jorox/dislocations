import miller_bravais as mb

class Atom(Vector):
    """This class represents an atom in a hexagonal crystal. Atom has the
    following properties:

    Attributes:
    index: a 4 integer array that specifies the unit cell to which it belongs,
    and it's basis index.
    xyz: a 3 float numpy array that specifies the coordinates in 'real' space
    """
    def __init__(self,index,pos):
        self.index = index
        self.xyz = pos

    def __get
