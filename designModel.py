import numpy as np
import copy

def designModel(model):
    model0 = copy.deepcopy(model)
    for ele in model0.element:
        ele.A = model0.xPhys[ele.tag]
        
    deleteElements =  np.argwhere(np.array(model0.xPhys)==0)
    model0.element = np.delete(model0.element, deleteElements, 0)

    model0.A = []
    deleteNodes = []
    for ele in model0.element:
        ele.startNode.elements.append(ele)
        ele.endNode.elements.append(ele)
        model0.A.append(model0.xPhys[ele.tag])
    for i in range(len(model0.node)):
        if len(model0.node[i].elements) == 0:
            deleteNodes.append(i)
    model0.node = np.delete(model0.node, deleteNodes, 0)
    print(len(model0.node), len(model0.element))
    return model0
    
