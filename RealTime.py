import matplotlib.pyplot as plt
import math
# Defining the constants
M = 100
b = 0.8
K = 200

def derivative(t, v, theta, disturbanceForce, M = M, K = K, b = b):
    p = (50-theta)/50
    return (K*p -v + disturbanceForce)*b/M

def disturbanceForce(t):
  return  0

def rungeKutta(disturbanceForce, theta, P1, t2, dx = 1):
    ''' theta - gas pedal angle
        P1 - tuple of x, y coordinate for initial point
        t2 - x coordinate for final point
        dx - step value(default - 1e-3)'''
    
    # taking the initial points
    t, v = P1

    # Arrays for storing the values corresponding to coordinates
    V = [v]
    T = [t]
    i = 0
    while(t + dx <= t2):
        # determining the position of the theta value
        if (i < len(theta[0]) - 1):
          if t+dx >= theta[1][i + 1]:
            i += 1
        
        # defining the K's
        k1 = derivative(t, v, theta[0][i], disturbanceForce(t))
        k2 = derivative(t + dx/2, v + (k1 * dx)/2, theta[0][i], disturbanceForce(t))
        k3 = derivative(t + dx/2, v + (k2*dx)/2, theta[0][i], disturbanceForce(t))
        k4 = derivative(t + dx, v + k3 * dx, theta[0][i], disturbanceForce(t))

        # Updating the new v
        vn = v + ((k1 + 2 * k2 + 2 * k3 + k4)*dx)/6

        # Updating the old values
        v = vn
        t += dx

        # Appending to the coordinate arrays
        V.append(v)
        T.append(t)
    return (T, V)

theta = [[1], [0]]


T, V = rungeKutta(disturbanceForce, theta, (0,0), 600)
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.5)
l, = plt.plot(T, V)

# Global n
n = 1000

def updateTheta(text):
  C = tuple(text.split())
  global theta, T, V
  theta = [[], []]
  for c in C:
      a, b = tuple(c.split(':'))
      a = float(a)
      b = float(b)
      theta[0].append(b)
      theta[1].append(a)
  print(theta)
  T, V = rungeKutta(disturbanceForce, theta, (0,0), n)
  l.set_ydata(V)
  ax.set_ylim(min(V), max(V))
  l.set_xdata(T)
  ax.set_xlim(min(T), max(T))
  axbox4 = plt.axes([0.1, 0.275, 0.8, 0.075])
  TextBox(axbox4, "", initial = f'Max velocity: {max(V)}')
  axbox3 = plt.axes([0.1, 0.2, 0.8, 0.075])
  saturationVelocity, t = getSaturationVelocity(T, V)
  if (saturationVelocity == None):
      TextBox(axbox3, "", initial = f'No Saturation in the given range')
  else:
      TextBox(axbox3, "", initial = f'Saturation velocity: {saturationVelocity}, Ts = {t}')
  plt.draw()

def updateN(text):
  global n, T, V
  # Updating the global n
  n = float(text)
  T, V = rungeKutta(disturbanceForce, theta, (0,0), n)
  l.set_ydata(V)
  ax.set_ylim(min(V), max(V))
  l.set_xdata(T)
  ax.set_xlim(min(T), max(T))
  axbox4 = plt.axes([0.1, 0.275, 0.8, 0.075])
  TextBox(axbox4, "", initial = f'Max velocity: {max(V)}')
  axbox3 = plt.axes([0.1, 0.2, 0.8, 0.075])
  saturationVelocity, t = getSaturationVelocity(T, V)
  if (saturationVelocity == None):
      TextBox(axbox3, "", initial = f'No Saturation in the given range')
  else:
      TextBox(axbox3, "", initial = f'Saturation velocity: {saturationVelocity}, Ts = {t}')
  plt.draw()

def updateDisturbance(text):
    global disturbanceForce
    def disturbanceForce(t):
        return eval(text)
    T, V = rungeKutta(disturbanceForce, theta, (0,0), n)
    l.set_ydata(V)
    ax.set_ylim(min(V), max(V))
    l.set_xdata(T)
    ax.set_xlim(min(T), max(T))
    axbox4 = plt.axes([0.1, 0.275, 0.8, 0.075])
    TextBox(axbox4, "", initial = f'Max velocity: {max(V)}')
    axbox3 = plt.axes([0.1, 0.2, 0.8, 0.075])
    saturationVelocity, t = getSaturationVelocity(T, V)
    if (saturationVelocity == None):
        TextBox(axbox3, "", initial = f'No Saturation in the given range')
    else:
        TextBox(axbox3, "", initial = f'Saturation velocity: {saturationVelocity}, Ts = {t}')
    plt.draw()
    

def getSaturationVelocity(T, V, tol = 1e-6):
    for i in range(len(T) - 1):
        if abs((V[i+1] - V[i])/V[i+1]) < tol:
            return V[i+1], T[i+1]
    return None, None

from matplotlib.widgets import TextBox
initial_text = "0:0"
plt.title("Velocity (vs) time")
plt.xlabel("time(t)")
plt.ylabel("velocity(Km/hr)")
axbox4 = plt.axes([0.1, 0.275, 0.8, 0.075])
TextBox(axbox4, "", initial = f'Max velocity: {max(V)}')
axbox3 = plt.axes([0.1, 0.2, 0.8, 0.075])
saturationVelocity, t = getSaturationVelocity(T, V)
if (saturationVelocity == None):
    TextBox(axbox3, "", initial = f'No Saturation in the given range')
else:
    TextBox(axbox3, "", initial = f'Saturation velocity: {saturationVelocity}, Ts = {t}')
axbox5 = plt.axes([0.1, 0.35, 0.8, 0.075])
axbox1 = plt.axes([0.1, 0.05, 0.8, 0.075])
axbox2 = plt.axes([0.1, 0.125, 0.5, 0.075])
thetaText = TextBox(axbox1, 'Control', initial=initial_text)
distForceText = TextBox(axbox5, 'Disturbance', initial='0')
distForceText.on_submit(updateDisturbance)
thetaText.on_submit(updateTheta)
thetaN = TextBox(axbox2, 'End Time', initial=str(n))
thetaN.on_submit(updateN)

plt.show()
