'''
Created on 16-01-2011

@author: lisu
'''

class Atom(object):
    '''
    Class holds atom objects with its coordinates and symbol
    '''
    def __init__(self, symbol, x = 0., y = 0. , z = 0.):
        self._x = x
        self._y = y
        self._z = z
        self._symbol = symbol
    @property    
    def x(self):
        return self._x
    @x.setter    
    def x(self, value):
        """X coordinate of atom in the 3d space"""
        self._x = value
    
    @property
    def y(self):
        return self._y
    @y.setter    
    def y(self, value):
        """Y coordinate of atom in the 3d space"""
        self._y = value
    
    @property
    def z(self):
        return self._z
    @z.setter    
    def z(self, value):
        """Z coordinate of atom in the 3d space"""
        self._z = value
    
    @property
    def symbol(self):
        return self._symbol
    @symbol.setter
    def symbol(self, value):
        """Symbol of the atom"""
        self._symbol = value
    @property
    def x0(self):
        """List access to position coordinates"""
        return [self.x, self.y, self.z]
    @x0.setter
    def x0(self, x0):
        self.x, self.y, self.z = x0
        
class Frame(object):
    
    def __init__(self):
        self._atomsList = []
        self._comment = ""
    @property
    def atoms(self):
        """Atoms in the frame"""
        return self._atomsList
    @property
    def comment(self):
        """Comment String"""
        return self._comment
    
    @comment.setter
    def comment(self, commentStr):
        self._comment = commentStr
  
    
    @property
    def boxSize(self):
        return NotImplemented