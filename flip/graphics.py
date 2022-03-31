from flip import geometry
import matplotlib as mpl
import matplotlib.tri as tri
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from flip import domain

class Graphics:
    def __init__(self):
        pass
    def add(self, obj):
        if (isinstance(obj, geometry.PolyLine)):
            self.addPolyline(obj)
        elif (isinstance(obj, geometry.TriSurf)):
            self.addTriangle(obj)
    def addPolyline(self, polyline):
        pass
    def addTriSurf(self, obj):
        pass
    def plot(self):
        pass



class MatplotlibGraphis2D(Graphics):
    def __init__(self):
        Graphics.__init__(self)
        plt.axes().set_aspect('equal')
        #self.ax = self.fig.gca(projection='3d')
        self.showcolorbar=False
    def addPolyline(self, pl):
        plt.plot(pl.x, pl.z, pl.fmt)
    def addTriSurf(self, nx, ny, cells, nodeVars=None, cellVars=None):
        fmt='ko-'
        t = tri.Triangulation(nx, ny, cells)
        if (cellVars is None) and (nodeVars is None):
            plt.triplot(t, 'ko-')
        elif (nodeVars is None):
            tcf=plt.tripcolor(t, facecolors = cellVars)
            plt.colorbar(tcf)
        #tcf=plt.tricontourf(t, tr.vals)
        #self.showcolorbar=True
        #else:
        #print(plt.triplot(tr.x, tr.y, [(0,1,2)], fmt, vmin=-100, vmax=100))
    def plot(self):
        if self.showcolorbar:
            plt.colorbar()
        plt.show()
