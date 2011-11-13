"""
This classes contain lattices for atoms
"""
__author__ =  'Mateusz Lis'
__version__=  '0.1'


from math import sin, cos, pi

class Position(object):
    """Class holds coordinates of position"""
    def __init__(self, x0 = [0.,0.]):
        self._x0 = x0
    @property
    def x0(self):
        """list of coordinates"""
        return self._x0
    @x0.setter
    def x0(self, x0):
        self._x0 = x0
    def __repr__(self):
	return str(self.x0)
       
        
class HexLattice(object):
    #FIXME: size!!!!
    def __init__(self, n = 8, x0 = (0,0), dx = (6.5,6.5), 
                 pbcSites = [[0,0,0], [-1,0,0], [-1,-1,0], [0,-1,0],[0,1,0],[-1,1,0], [1,0,0], [1,1,0], [1,-1,0]]):
        """
        Creates instance of hexagonal lattice class.
        n - denotes number of particles in one line of lattice (it must be equilateral)
        lattice contains n x n particles
        x0 - denotes top left corner coordinates of the hexagonal lattice
        pbcSites - denotes sites of square frame that have to be copied using pbc
        """
        self._x0 = x0
        self._n = n
        self._positions = []
	self._angle = pi/3. #FIXME: Fixed constant
        self._dx = dx
        self._setPositions()
        self._pbcSites = pbcSites
         
    @property    
    def x0(self):
        """Top left corner of lattice"""
        return self._x0
    @x0.setter
    def x0(self, x0):
        self._movePos(x0)
        self._x0 = x0
        
    def _movePos(self, x0):
        """Translates all positions from self.x0 to x0"""
        diff = [x0[i] - self._x0[i] for i in range(len(self._x0))]
        for pos in self.positions:
            pos.x0 = [pos.x0[i] + diff[i] for i in range(len(pos.x0))]
    @property
    def n(self):
        """Number of particles in one line of lattice"""
        return self._n
    
    @property
    def dx(self):
        """vector of lengths of edges"""
        return self._dx
    @dx.setter
    def dx(self, dx):
        return NotImplemented
    @property
    def boxVectors(self):
        v1 = [self.dx[0], 0, 0]
        
        v2 = [ - cos(self._angle) * self.dx[1], sin(self._angle) * self.dx[1],0]
        
        v3 = [0, 0, 0]
        return [v1, v2, v3]
        
    @property
    def positions(self):
        """List of positions in lattice"""
        return self._positions
    
    @property
    def pbcSites(self):
        return self._pbcSites
    
    
    def _setPositions(self):
        """create position field for each particle in lattice and set its position"""
        xLine0 = self.x0[0] #+ 0.5 * (self.dx[0] / self.n - cos(self._angle) * self.dx[1] / self.n) 
        yLine0 = self.x0[1] #+ 0.5 * self.dx[1] / self.n * sin(self._angle)
        print sin(self._angle), "sin"
	print cos(self._angle), "cos"
        for _ in range(self.n):
            for x in range(self.n):
                self.positions.append(Position([xLine0 + x * self.dx[0] / self.n, yLine0]))
            xLine0 -= cos(self._angle) * self.dx[1] / self.n
            yLine0 += sin(self._angle) * self.dx[1] / self.n
	print self.positions[0], self.positions[1], "Positions"
    
