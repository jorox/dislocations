import miller_bravais as mb

class Atom(Vector):
    """This class represents an atom in a hexagonal crystal. Atom inherits Vector and adds
    to it some "real" space features such as XYZ coordinates, unit cell indices.
    These two new properties are **calculated** by specifying the Miller-Bravais directions 
    along which the main X, Y, and Z axes lie or in other words by giving the Vector
    representations of the î,ĵ,k vectors 

    Attributes:
    index: a 4 integer array that specifies the unit cell to which it belongs,
    and it's basis index.
    xyz: a 3 float numpy array that specifies the coordinates in 'real' space
    """
    xyz = []
    
