import numpy as np


def RungeKutta4(inVec, t, dt, computeDByDt, cookie):
    vecLen = len(inVec)
    outVec = np.zeros(vecLen)
    k1 = np.zeros(vecLen)
    k2 = np.zeros(vecLen)
    k3 = np.zeros(vecLen)
    k4 = np.zeros(vecLen)

    k1 = dt * computeDByDt(inVec, t, cookie)
    k2 = dt * computeDByDt(inVec + (k1/2.), t + (dt/2.), cookie)
    k3 = dt * computeDByDt(inVec + (k2/2.), t + (dt/2.), cookie)
    k4 = dt * computeDByDt(inVec + k3, t + dt, cookie)

    outVec = inVec + (k1 + 2*k2 + 2*k3 + k4) / 6.
    return outVec


def ODESolve(method, dt, inVec, computeDByDt, cookie, params):
    [nEquations, tEnd, maxSteps, batchSteps, filePrefix] = params
    # TODO: add batched handling, flushing data to file
    data = np.zeros((maxSteps, nEquations))

    if (method == 'RungeKutta4'):
        ODESolver = RungeKutta4
    else:
        print("Unknown ODESolver %s" % (method))
        return (-1., data[:n])

    t = 0.
    vec = inVec

    nSteps = 0
    data[0] = vec
    # TODO: start timer
    while (nSteps < maxSteps - 1 and t < tEnd):
        vec = RungeKutta4(vec, t, dt, computeDByDt, cookie)
        nSteps += 1
        data[nSteps] = vec
        t += dt
    # TODO: end timer

    # TODO: compute timer and add to return data
    return (t, data[:nSteps], 0.)
