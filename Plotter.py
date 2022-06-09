import matplotlib.pyplot as plt

def plotNodes(model):
    fig2, ax2 = plt.subplots()
    for node in model.node:
        ax2.scatter(node.x, node.y, color='b', s=10)

def plotElement(model):
    plotELe = []
    fig1, ax1 = plt.subplots()
    for node in model.node:
        ax1.scatter(node.x, node.y, color='b', s=10)
    for element in model.element:
        plotELe.append( ax1.plot([element.startNode.x, element.endNode.x], [element.startNode.y, element.endNode.y], 'b'))
        # ax1.plot([element.startNode.x, element.endNode.x], [element.startNode.y, element.endNode.y], 'b')
    return plotELe

