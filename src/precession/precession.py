import math
import numpy as np
import matplotlib.pyplot as plt

from ode_solve import *


class Config:
  def __init__(self, tEnd, maxSteps, batchSteps, filePrefix, onlyMercury, decoupledMercury, grMercury, enableConvergenceTest):
    self.tEnd = tEnd
    self.maxSteps = maxSteps
    self.batchSteps = batchSteps
    self.filePrefix = filePrefix
    self.onlyMercury = onlyMercury
    self.decoupledMercury = decoupledMercury
    self.grMercury = grMercury
    self.enableConvergenceTest = enableConvergenceTest

class Planet:
  def __init__(self, config, name, M, T, R):
    self.name = name
    self.M = M/(2.e+30)
    self.T = T
    self.R = R
    self.V = (2 * math.pi * self.R) / self.T
    self.L = self.R * self.V
    self.GMS = self.R * self.V**2
    if (config.decoupledMercury):
      self.GMP = 0.
    else:
      self.GMP = self.GMS * self.M

class Mercury:
  def __init__(self, config, GMS):
    self.GMS = GMS
    self.M = (2.4e+23)/(2.e+30)
    self.T = 0.24
    self.e = 0.206
    self.a = 0.387
    self.RMin = self.a * (1 - self.e)
    self.RMax = self.a * (1 + self.e)
    self.vMax = math.sqrt((((1 + self.e) * (1 + self.M)) / self.RMin) * self.GMS)
    self.L = self.a * (1 - self.e) * self.vMax
    if (config.decoupledMercury):
      self.GMM = 0.
    else:
      self.GMM = self.GMS * self.M
    if (config.grMercury):
      # Constant for General Relativistic correction on Mercury due to Sun
      # This is the value that goes into the force equation term as:
      #  +3*GMS*alpha/r^4
      # And update to potential term in Lagrangian will be something like:
      #  -GMS*alpha/r^3
      # Where 3*alpha value taken from Giordano & Nakanishi is: 1.1e-8
      #self.GMR = self.GMS * (1.1e-8)
      self.GMR = self.GMS / (6.25e+4)**2
    else:
      self.GMR = 0.
    print("rmin = %.7f, vmax = %.7f, L = %.7f, GMS = %.7f, GMM = %.12f, GMR = %.12f" %(self.RMin, self.vMax, self.L, self.GMS, self.GMM, self.GMR))


def InitialValues(config, planet, mercury):
  # initial values are arranged as:
  # [rP(0), thetaP(0), vP(0), omegaP(0), rM(0), thetaM(0), vM(0), omegaM(0)]
  # note that v = dr_dt and omega = dtheta_dt
  # place the planet and mercury with phase difference of pi
  if (config.onlyMercury):
    return [0., 0., 0., 0., mercury.RMin, 0., 0., (mercury.vMax / mercury.RMin)]
  else:
    return [planet.R, math.pi, 0., (planet.L / planet.R**2), mercury.RMin, 0., 0., (mercury.vMax / mercury.RMin)]

def ComputeDByDt(inVec, t, cookie):
  [config, planet, mercury] = cookie

  # incoming values are arranged as:
  # [rP(t-), thetaP(t-), vP(t-), omegaP(t-), rM(t-), thetaM(t-), vM(t-), omegaM(t-)]
  # note that v = dr_dt and omega = dtheta_dt
  [rP, thetaP, vP, omegaP, rM, thetaM, vM, omegaM] = inVec

  thetaMP = thetaM - thetaP
  rMP = math.sqrt(rM**2 + rP**2 - 2 * rM * rP * math.cos(thetaMP))
  rMP3 = rMP**3
  alphaP = planet.GMP / rMP3
  alphaM = mercury.GMM / rMP3

  if (config.onlyMercury):
    drP_dt = 0.
    dthetaP_dt = 0.
    dvP_dt = 0.
    domegaP_dt = 0.
  else:
    drP_dt = vP
    dthetaP_dt = omegaP
    dvP_dt = (rP * omegaP**2) - (planet.GMS / rP**2) - (alphaM * (rP - rM * math.cos(thetaMP)))
    domegaP_dt = ((-2. * vP * omegaP) / rP) + (alphaM * (rM / rP) * math.sin(thetaMP))

  drM_dt = vM
  dthetaM_dt = omegaM
  dvM_dt = (rM * omegaM**2) - (mercury.GMS / rM**2)
  domegaM_dt = ((-2. * vM * omegaM) / rM)
  if (config.grMercury):
    #dvM_dt -= (mercury.GMR / rM**4)
    rMRed2 = rM - (2 * mercury.GMR)
    rMRed3 = rM - (3 * mercury.GMR)
    dvM_dt -= (2 * mercury.GMR * omegaM**2)
    dvM_dt += ((2 * mercury.GMS * mercury.GMR) / rM**3)
    dvM_dt += ((3 * mercury.GMR * vM**2) / (rM * rMRed2))
    domegaM_dt *= (rMRed3 / rMRed2)
  dvM_dt -= (alphaP * (rM - rP * math.cos(thetaMP)))
  domegaM_dt -= (alphaP * (rP / rM) * math.sin(thetaMP))

  # return values are arranged as:
  # [drP_dt(t), dthetaP_dt(t), dvP_dt(t), domegaP_dt(t), drM_dt(t), dthetaM_dt(t), dvM_dt(t), domegaM_dt(t)]
  return np.array([drP_dt, dthetaP_dt, dvP_dt, domegaP_dt, drM_dt, dthetaM_dt, dvM_dt, domegaM_dt])


def Plot(config, targetTolerance, tEnd, dt, data, planet, mercury):
  # incoming values are arranged as:
  # [rP(t), thetaP(t), vP(t), omegaP(t), rM(t), thetaM(t), vM(t), omegaM(t)]
  # note that v = dr_dt and omega = dtheta_dt
  nRows = len(data)

  if (not(config.onlyMercury)):
    plt.figure("%s orbit, data for %d points" %(planet.name, nRows))
    plt.grid()
    plt.axes().set_aspect('equal')
    plt.xlabel('X')
    plt.ylabel('Y')
    xValues = np.zeros(nRows)
    yValues = np.zeros(nRows)
    rPMin = planet.R - targetTolerance
    rPMax = planet.R + targetTolerance
    t = 0.
    print("%s out of tolerance orbit data for last 100 yrs (tEnd %.7f, rP %.9f, rPMin %.9f, rPMax %.9f): count, t, rP, thetaP, phiP" %(planet.name, tEnd, planet.R, rPMin, rPMax))
    for i in range(nRows):
      rP = (data[i])[0]
      thetaP = (data[i])[1]
      xValues[i] = rP * math.cos(thetaP)
      yValues[i] = rP * math.sin(thetaP)
      if (t >= (tEnd - 100.) and not(rP >= rPMin or rP <= rPMax)):
        phiP = thetaP % (2 * math.pi)
        print("%d, %.9f, %.9f, %.9f, %.9f" %(i, t, rP, thetaP, phiP))
      t += dt
    print("=========")
    plt.plot(xValues, yValues)
    plt.show()

  plt.figure("Mercury orbit, GR effects %s, data for %d points" %(["disabled", "enabled"][config.grMercury], nRows))
  plt.grid()
  plt.axes().set_aspect('equal')
  plt.xlabel('X')
  plt.ylabel('Y')
  xValues = np.zeros(nRows)
  yValues = np.zeros(nRows)
  rMMaxMin = mercury.RMax - targetTolerance
  rMMaxMax = mercury.RMax + targetTolerance
  t = 0.
  rMMax1 = 0.
  rMMax2 = 0.
  print("Mercury perihelion data for first and last 50 yrs (tEnd %.7f, rMMax %.9f, rMMaxMin %.9f, rMMaxMax %.9f): count, t, rMMax, thetaM, phiM, phiMDelta" %(tEnd, mercury.RMax, rMMaxMin, rMMaxMax))
  for i in range(nRows):
    rM = (data[i])[4]
    thetaM = (data[i])[5]
    xValues[i] = rM * math.cos(thetaM)
    yValues[i] = rM * math.sin(thetaM)
    if ((t <= 50.) and (rMMax1 < rM)):
      rMMax1 = rM
      phiM = thetaM % (2 * math.pi)
      print("%d, %.9f, %.9f, %.9f, %.9f, %.9f" %(i, t, rM, thetaM, phiM, phiM - math.pi))
    if ((t >= (tEnd - 50.)) and (rMMax2 < rM)):
      rMMax2 = rM
      phiM = thetaM % (2 * math.pi)
      print("%d, %.9f, %.9f, %.9f, %.9f, %.9f" %(i, t, rM, thetaM, phiM, phiM - math.pi))
    t += dt
  print("=========")
  plt.plot(xValues, yValues)
  plt.show()


if __name__ == "__main__":
  # TODO: read config values from config file
  # confg variables --
  #   tEnd, maxSteps, batchSteps, filePrefix,
  #   onlyMercury, decoupledMercury, grMercury
  #   enableConvergenceTest
  config = Config(450, 10000000, 10000000, 'precession', False, False, True, True)

  # TODO: include dt, dtFactor, targetTolerance into config later
  targetTolerance = 1.0e-7
  #dt = 0.001
  dt = 0.00005
  dtFactor = 2

  nEquations = 8
  params = [nEquations, config.tEnd, config.maxSteps, config.batchSteps, config.filePrefix]

  #planet = Planet(config, 'Venus', 4.87e+24, 0.616, 0.7)
  #planet = Planet(config, 'Earth', 5.97e+24, 1., 1.)
  planet = Planet(config, 'Jupiter', 1.9e+27, 12., 5.203)
  mercury = Mercury(config, planet.GMS)
  cookie = [config, planet, mercury]

  dt1 = dt
  print("\n\n=========")
  print("calling ODESolve for dt1 %.7f ..." %(dt1))
  (tEnd1, data1, timer1) = ODESolve('RungeKutta4', dt1, InitialValues(config, planet, mercury), ComputeDByDt, cookie, params)
  print("completed ODESolve for dt1, total-run-time %d, nSteps %d, tEnd %0.7f, final state %s" %(timer1, len(data1), tEnd1, data1[-1]))
  # TODO: print statistics on timimg
  Plot(config, targetTolerance, tEnd1, dt1, data1, planet, mercury)

  if (config.enableConvergenceTest):
    dt2 = dt1 * dtFactor
    print("\n\n=========")
    print("calling ODESolve for dt2 %.7f ..." %(dt2))
    (tEnd2, data2, timer2) = ODESolve('RungeKutta4', dt2, InitialValues(config, planet, mercury), ComputeDByDt, cookie, params)
    print("completed ODESolve for dt2, total-run-time %d, nSteps %d, tEnd %0.7f, final state %s" %(timer2, len(data2), tEnd2, data2[-1]))
    # TODO: print statistics on timimg
    Plot(config, targetTolerance, tEnd2, dt2, data2, planet, mercury)

    dt3 = dt2 * dtFactor
    print("\n\n=========")
    print("calling ODESolve for dt3 %.7f ..." %(dt3))
    (tEnd3, data3, timer3) = ODESolve('RungeKutta4', dt3, InitialValues(config, planet, mercury), ComputeDByDt, cookie, params)
    print("completed ODESolve for dt3, total-run-time %d, nSteps %d, tEnd %0.7f, final state %s" %(timer3, len(data3), tEnd3, data3[-1]))
    # TODO: print statistics on timimg
    Plot(config, targetTolerance, tEnd3, dt3, data3, planet, mercury)

    diff12 = (data1[::2*dtFactor])[:len(data3)] - (data2[::dtFactor])[:len(data3)]
    diff23 = (data2[::dtFactor])[:len(data3)] - data3

    # TODO: print detailed data on convergence etc.
    print("========= diff12 =========")
    print(diff12)
    print("========= diff23 =========")
    print(diff23)
