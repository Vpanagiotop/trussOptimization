import numpy as np
class Element:
    def __init__(self, tag, E, A, startNode, endNode):
        self.tag = tag
        self.startNode = startNode
        self.endNode = endNode
        self.E = E
        self.A = A
        self.nodeTag = [startNode.tag, endNode.tag]
        self.dof = np.array([startNode.dof, endNode.dof]).reshape(1,4)[0]
        self.length, self.cos, self.sin = getElementInfo(startNode, endNode)
        self.nominalStiffness = getNominalStiffness(self.length, self.cos, self.sin, E, A)
        self.stiffness = self.nominalStiffness.copy()

def getElementInfo(startNode, endNode):
    length = np.sqrt((startNode.x - endNode.x)**2 + (startNode.y - endNode.y)**2)
    Lx = endNode.x - startNode.x
    Ly = endNode.y - startNode.y
    cos = Lx/length
    sin = Ly/length
    return length, cos, sin

def getNominalStiffness(length, c, s, E, A):
    return E * A / length * np.array([
        [c**2,  c*s,   -c**2, -c*s ],
	    [c*s,   s**2,  -c*s,  -s**2],
        [-c**2, -c*s,  c**2,  c*s ],
	    [-c*s,  -s**2, c*s,  s**2]
    ])