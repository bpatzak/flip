class Geometry:
    def __init__(self):
        pass

class PolyLine(Geometry):
    def __init__(self, x,y,z, fmt=""):
        Geometry.__init__(self)
        self.x=x
        self.y=y
        self.z=z
        self.fmt = fmt

class TriSurf(Geometry):
    def __init__(self, nodes, cells):
        Geometry.__init__(self)
        self.nodes=nodes
        self.cells=cells
        


