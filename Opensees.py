
import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import openseespy.postprocessing as opsv

def Analysis(Ai, model, plotmodel = False):
    E = 210000000
    ops.wipe()
    ops.model('basic', '-ndm', 2, '-ndf', 2)

    for i in range(len(model.node)):
        node = model.node[i]
        ops.node(node.tag, node.x, node.y)
        if node.xFixed == True: x = 1
        else: x = 0
        if node.yFixed == True: y = 1
        else: y = 0
        ops.fix(node.tag, x, y)
    # ops.uniaxialMaterial("Steel01", 1, 250000, E, 0)
    ops.uniaxialMaterial('Elastic', 1, E)

    for i in range (len(model.element)):
        element = model.element[i]
        ops.element("Truss", element.tag, element.startNode.tag, element.endNode.tag, Ai[i], 1)

    ops.timeSeries('Linear', 1)
    ops.pattern("Plain", 1, 1)

    for i in range(len(model.node)):
        if model.node[i].load.Py!=0:
            ops.load(model.node[i].tag, 0, -100)

    ops.constraints('Plain')
    ops.numberer('RCM')
    ops.system('BandSPD')
    ops.test('NormDispIncr', 1.0e-6, 6, 2)
    ops.algorithm('Linear')
    ops.integrator('LoadControl', 1)
    ops.analysis('Static')
    ops.analyze(1)

    if plotmodel == True:
        opsv.plot_defo(sfac=0)
        plt.grid()
        opsv.plot_model()
        plt.grid()
    Ni = np.array([])
    for i in range(1,len(model.element)+1):
        Ni = np.append(Ni,ops.basicForce(i)[0])

    u = np.array([])
    v = np.array([])
    for i in range(len(model.node)):
        nodeTag = model.node[i].tag
        model.node[i].u = ops.nodeDisp(nodeTag,1)
        model.node[i].v = ops.nodeDisp(nodeTag,2)
        u = np.append(u,ops.nodeDisp(nodeTag,1))
        v = np.append(v,ops.nodeDisp(nodeTag,2))
    model.u = u
    model.v = v
    strainEnergy = np.array([])
    for i in range(len(model.element)):
        ele = model.element[i]
        ele.Ni = ops.basicForce(ele.tag)[0]
        ui = ele.startNode.u
        uj = ele.endNode.u
        vi = ele.startNode.v
        vj = ele.endNode.v
        ele.elogation = ele.cos * (uj-ui) + ele.sin * (vj-vi)
        ele.strainEnergy = np.round(ele.elogation * ele.Ni * 0.5, 6)
        strainEnergy = np.append(strainEnergy, ele.strainEnergy)
    model.strainEnergy = strainEnergy
    return Ni
