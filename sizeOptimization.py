import numpy as np
from scipy.optimize import minimize
from designModel import designModel
from scipy.optimize import minimize
# from Opensees import Analysis
def sizeOptimization(model):
    model0 = designModel(model)
    L = np.array([ele.length for ele in model.element])
    A = np.array([1.e-5 for i in model.element], dtype=float)

    def stressDistribution(x):
        Ni = Analysis(x, model)
        return Ni/x
    
    def uperLimit(x):
        stress = stressDistribution(x)
        return -stress + 235000

    def lowerLimit(x):
        stress = stressDistribution(x)
        return stress + 220000

    def maxDeflection(x):
        Ni= Analysis(x, model)
        vi = model.v
        return vi

    def DeflectionLimit(x):
        vi = maxDeflection(x)
        return -max(abs(vi)) + 0.005
    
    def Deflectionbnd(x):
        Ni= Analysis(x, model)
        ui = model.u
        return ui[6]

    def objectiveFunction(x):
        volume = L * x
        return sum(volume)

    cons = ({'type': 'ineq', 'fun': uperLimit},
            {'type': 'ineq', 'fun': lowerLimit},
            {'type': 'ineq', 'fun': DeflectionLimit},
            {'type': 'eq', 'fun': Deflectionbnd},
    )
    
    bnds = [(ele.A,0.9) for ele in model.element]
    bnds = tuple(bnds)
    
    res = minimize(objectiveFunction, model.A, method='SLSQP', bounds=bnds, constraints=cons, options={'maxiter': 1000, 'gtol': 1e-6, 'disp': True})
    A = res.x
    print(res)
    Ni= Analysis(A, model, plotmodel=False)
    stress = Ni/A
    for i in range (len(model.element)):
        model.element[i].A = res.x[i]
    model.modify()
    print('----------------------------- Size Optimization -----------------------------')
    print(res)
    np.set_printoptions(suppress = True, precision=5)
    print('----------------------------- Stress (MPa)-----------------------------')
    stress = Ni/A
    print(stress / 1000)
    print('----------------------------- max Deflection (m)-----------------------------')
    print(max(abs(model.u)))
    return (res, stress/1000)