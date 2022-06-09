import numpy as np
from scipy.sparse import coo_matrix, linalg
from Node import Node
from Element import Element
class Generator:
    def __init__(self, length, height, rmin: int, mesh={'x': float, 'y': float}, Amin:float = 1.e-9, Amax:float = 1, pen=3.0):
        self.length = length
        self.height = height
        self.mesh = mesh
        self.Amin = Amin
        self.Amax = Amax
        self.pen = pen
        points = []
        xCoordinates = np.linspace(0, length, int(length/mesh['x'] + 1), dtype=float)
        yCoordinates = np.linspace(0, height, int(height/mesh['y'] + 1), dtype=float)
        for x in xCoordinates:
            for y in yCoordinates:
                points.append([x,y])
        self.node = assignNodes(points)
        self.dofs = assigndegreesOfFreedom(self.node)
        self.element = assignElements(self, yCoordinates)
        self.elementDofs = assignElementDofs(self)
        self.Kel, self.Kmatrix = getElasticStiffness(self)
        
    def format(self):
        self.fixedDofs, self.freeDofs = setBoundaryConditions(self)
        self.F = setNodalLoads(self)

    def modify(self, xPhys):
        getElementStiffness(self, xPhys)
        self.Kel, self.Kmatrix = getElasticStiffness(self)

    def solveFE(self):
        Kff = self.Kmatrix[self.freeDofs,:][:,self.freeDofs]
        u = np.zeros((len(self.dofs),1))
        u[self.freeDofs,0] = linalg.spsolve(Kff,self.F[self.freeDofs,0])
        return u

def assignNodes(points):
    node = []
    for i in range(len(points)):
        tag = i
        x, y = points[i][0], points[i][1]
        node.append(Node(tag, x, y)) 
    return node

def assigndegreesOfFreedom(node):
    dof = np.linspace(0,len(node)*2 - 1, len(node)*2, dtype=int)
    for i in range(len(node)):
        node[i].dof = [dof[2 * i], dof[2 * i + 1]]
    return dof

def assignElements(self, yCoordinates):
    element = []
    tag = 0
    for node in self.node[:-1]:
        startNode = node
        if startNode.x == self.length:
            endNodeTags = np.array([node.tag + 1])
        else:
            if startNode.y == self.height:
                endNodeTags = np.array([node.tag + len(yCoordinates), node.tag + len(yCoordinates) -1])
            elif startNode.y == 0:
                endNodeTags = np.array([node.tag+1, node.tag + len(yCoordinates), node.tag + len(yCoordinates) +1])
            else:
                endNodeTags = np.array([node.tag+1, node.tag + len(yCoordinates), node.tag + len(yCoordinates) +1, node.tag + len(yCoordinates) -1])
        endNodes = np.array(self.node)[endNodeTags]
        for endNode in endNodes:
            element.append(Element(tag, 1, 1, startNode, endNode))
            tag = tag+1
    return element

def assignElementDofs(self):
    elementDofs = np.zeros((len(self.element),4), dtype =int)
    for ele in self.element:
        elementDofs[ele.tag] = np.array([ele.dof])
    return elementDofs

def getElasticStiffness(self):
    Kel = np.zeros((self.dofs.size, self.dofs.size), float)
    for ele in self.element:
        for iEle in range(len(ele.dof)):
            i = int(ele.dof[iEle])
            for jEle in range(len(ele.dof)):
                j = int(ele.dof[jEle])
                Kel[i,j] = Kel[i,j] + float(ele.stiffness[iEle, jEle])
    Kmatrix = coo_matrix(Kel).tocsc()
    return Kel, Kmatrix

def setBoundaryConditions(self):
    fixedDofs = np.empty((0,0), dtype=int)
    for node in self.node:
        if node.xFixed == True:
           fixedDofs =  np.append(fixedDofs,node.dof[0])
        if node.yFixed == True:
           fixedDofs =  np.append(fixedDofs,node.dof[1])
    freeDofs = np.setdiff1d(self.dofs,fixedDofs)
    return fixedDofs, freeDofs

def setNodalLoads(self):
    F = np.zeros((len(self.dofs),1))
    for node in self.node:
        if node.Px != 0:
            F[node.dof[0]] = node.Px
        if node.Py!= 0:
            F[node.dof[1]] = node.Py
    return F

def getElementStiffness(self, xPhys):
    for ele in self.element:
        ele.stiffness = ele.nominalStiffness * (self.Amin+(xPhys[ele.tag])**self.pen*(self.Amax-self.Amin))