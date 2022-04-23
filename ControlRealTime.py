from matplotlib.widgets import TextBox, Slider
import matplotlib.pyplot as plt
import math  # Do not change the imports, as some can be used in real time

# Defining the constants
M = 100         # Mass
b = 20          # Resistive Coefficient
Ke = 0.6 * 100  # Engine Gain ( x 100 So that the 'p' represents the fraction)
Kc = 1e-3
P = 50          # Maximum angle that can be pressed
vo = 0        # Initial velocity 

################ Defining the equations and forces ####################


def equationNoControl(t, v, theta, disturbanceForce, M=M, K=Ke, b=b):
    # Theta denotes the angle pressed
    p = (theta)/50  # p denotes the fraction of the amount that can be pressed
    return (K*p - v + disturbanceForce(t)/b)*b/M

def vd(t):
    '''Rquired velocity'''
    return 100

def disturbanceForce(t):
    return 0

def equationWithControl(t, v, disturbanceForce, M=M, Ke=Ke, Kc = Kc, b=b):
    # Theta denotes the angle pressed
    return (Kc*Ke/M) * vd(t) - v * (b + Kc*Ke)/M + disturbanceForce(t)/M

###################################################################


################################ Utility functions ########################

flag = False # Just to prevent warnings 
# Range Kutta Method


def rungeKutta(disturbanceForce, theta, P1, t2, dx=1, flag = flag, Kc = Kc):
    ''' theta - gas pedal angle
        P1 - tuple of x, y coordinate for initial point
        t2 - x coordinate for final point
        dx - step value(default - 1e-3)'''

    # taking the initial points
    t, v = P1
    
    # Arrays for storing the values corresponding to coordinates
    V = [v * 18/5]
    T = [t]
    i = 0
    while(t + dx <= t2):
        # determining the position of the theta value
        if (i < len(theta[0]) - 1):
            if t+dx >= theta[1][i + 1]:
                i += 1

        # defining the K's
        if (not flag):
            derivative = equationNoControl
            k1 = derivative(t, v, theta[0][i], disturbanceForce)
            k2 = derivative(t + dx/2, v + (k1 * dx)/2,
                                theta[0][i], disturbanceForce)
            k3 = derivative(t + dx/2, v + (k2*dx)/2,
                                theta[0][i], disturbanceForce)
            k4 = derivative(t + dx, v + k3 * dx,
                               theta[0][i], disturbanceForce)
        else:
            derivative = equationWithControl
            k1 = derivative(t, v, disturbanceForce, Kc = Kc)
            k2 = derivative(t + dx/2, v + (k1 * dx)/2,
                                disturbanceForce, Kc = Kc)
            k3 = derivative(t + dx/2, v + (k2*dx)/2,
                                disturbanceForce, Kc = Kc)
            k4 = derivative(t + dx, v + k3 * dx,
                               disturbanceForce, Kc = Kc)

        # Updating the new v
        vn = v + ((k1 + 2 * k2 + 2 * k3 + k4)*dx)/6

        # Updating the old values
        v = vn
        t += dx

        # Appending to the coordinate arrays
        V.append(v * 18/5)
        T.append(t)
    return (T, V)


# Utility for updating the vlaues of theta


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
    T, V = rungeKutta(disturbanceForce, theta, (0, vo), n, flag = flag)
    l.set_ydata(V)
    ax.set_ylim(min(V), max(V))
    l.set_xdata(T)
    ax.set_xlim(min(T), max(T))
    axbox4 = plt.axes([0.1, 0.2, 0.8, 0.075])
    axbox3 = plt.axes([0.1, 0.125, 0.8, 0.075])
    TextBox(axbox4, "", initial=f'Max velocity: {max(V)}')
    saturationVelocity, t = getSaturationVelocity(T, V)
    if (saturationVelocity == None):
        TextBox(axbox3, "", initial=f'No Saturation in the given range')
    else:
        TextBox(
            axbox3, "", initial=f'Saturation velocity: {saturationVelocity}, Ts = {t}')
    plt.draw()

# Utility for updating the value of n


def updateN(text):
    global n, T, V
    # Updating the global n
    n = float(text)
    T, V = rungeKutta(disturbanceForce, theta, (0, vo), n, flag = flag)
    l.set_ydata(V)
    ax.set_ylim(min(V), max(V))
    l.set_xdata(T)
    ax.set_xlim(min(T), max(T))
    axbox4 = plt.axes([0.1, 0.2, 0.8, 0.075])
    axbox3 = plt.axes([0.1, 0.125, 0.8, 0.075])
    TextBox(axbox4, "", initial=f'Max velocity: {max(V)}')
    saturationVelocity, t = getSaturationVelocity(T, V)
    if (saturationVelocity == None):
        TextBox(axbox3, "", initial=f'No Saturation in the given range')
    else:
        TextBox(
            axbox3, "", initial=f'Saturation velocity: {saturationVelocity}, Ts = {t}')
    plt.draw()

# Utility for updating the value of Disturbance


def updateDisturbance(text):
    global disturbanceForce, T, V

    def disturbanceForce(t):
        return eval(text)
    T, V = rungeKutta(disturbanceForce, theta, (0, vo), n, flag = flag)
    l.set_ydata(V)
    ax.set_ylim(min(V), max(V))
    l.set_xdata(T)
    ax.set_xlim(min(T), max(T))
    axbox4 = plt.axes([0.1, 0.2, 0.8, 0.075])
    axbox3 = plt.axes([0.1, 0.125, 0.8, 0.075])
    TextBox(axbox4, "", initial=f'Max velocity: {max(V)}')
    saturationVelocity, t = getSaturationVelocity(T, V)
    if (saturationVelocity == None):
        TextBox(axbox3, "", initial=f'No Saturation in the given range')
    else:
        TextBox(
            axbox3, "", initial=f'Saturation velocity: {saturationVelocity}, Ts = {t}')
    plt.draw()
# Utility Function for updating the value of Kc

def updateKc(val):
    global Kc, V, T
    Kc = val
    T, V = rungeKutta(disturbanceForce, theta, (0, vo), n, flag = flag, Kc = Kc)
    l.set_ydata(V)
    ax.set_ylim(min(V)- 10, max(V) + 10)
    l.set_xdata(T)
    ax.set_xlim(min(T) - 10, max(T) + 10)
    axbox4 = plt.axes([0.1, 0.2, 0.8, 0.075])
    axbox3 = plt.axes([0.1, 0.125, 0.8, 0.075])
    TextBox(axbox4, "", initial=f'Max velocity: {max(V)}')
    saturationVelocity, t = getSaturationVelocity(T, V)
    if (saturationVelocity == None):
        TextBox(axbox3, "", initial=f'No Saturation in the given range')
    else:
        TextBox(
            axbox3, "", initial=f'Saturation velocity: {saturationVelocity}, Ts = {t}')
    plt.draw()

# Utility for updating the initial velocity


def updateVelocity(text):
    text = float(text)
    global vo
    vo = text * 5/18
    T, V = rungeKutta(disturbanceForce, theta, (0, vo), n)
    l.set_ydata(V)
    ax.set_ylim(min(V), max(V))
    l.set_xdata(T)
    ax.set_xlim(min(T), max(T))
    axbox4 = plt.axes([0.1, 0.2, 0.8, 0.075])
    TextBox(axbox4, "", initial=f'Max velocity: {max(V)}')
    axbox3 = plt.axes([0.1, 0.125, 0.8, 0.075])
    saturationVelocity, t = getSaturationVelocity(T, V)
    if (saturationVelocity == None):
        TextBox(axbox3, "", initial=f'No Saturation in the given range')
    else:
        TextBox(
            axbox3, "", initial=f'Saturation velocity: {saturationVelocity}, Ts = {t}')
    plt.draw()

# Utility for getting the saturation velocity


def getSaturationVelocity(T, V, tol=1e-5):
    for i in range(len(T) - 2, -1, -1):
        if abs((V[i+1] - V[i])/V[i+1]) < tol:
            for j in range(i, -1, -1):
                if abs((V[j+1] - V[j])/V[j+1]) > tol:
                    return V[j+1], T[j+1]
    return None, None


#########################################################################

########################## With Control ##########################
# Initial values
flag = True # with control
theta = [[50], [0]]
n = 1000
T, V = rungeKutta(disturbanceForce, theta, (0, vo), n, flag = flag)

# Plotting
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.5)  # Leaving space at the bottom
l, = plt.plot(T, V)
plt.title("Velocity (vs) time")
plt.xlabel("time(s)")
plt.ylabel("velocity(Km/hr)")
axbox4 = plt.axes([0.1, 0.2, 0.8, 0.075])
TextBox(axbox4, "", initial=f'Max velocity: {max(V)}')
axbox3 = plt.axes([0.1, 0.125, 0.8, 0.075])
saturationVelocity, t = getSaturationVelocity(T, V)
if (saturationVelocity == None):
    TextBox(axbox3, "", initial=f'No Saturation in the given range')
else:
    TextBox(
        axbox3, "", initial=f'Saturation velocity: {saturationVelocity}, Ts = {t}')
axbox5 = plt.axes([0.1, 0.275, 0.3, 0.075])
axbox2 = plt.axes([0.1, 0.05, 0.5, 0.075])
axbox6 = plt.axes([0.6, 0.275, 0.3, 0.075])
voText = TextBox(axbox6, 'Vo', initial=str(vo))
voText.on_submit(updateVelocity)
axslide1 = plt.axes([0.1, 0.35, 0.8, 0.03])
Kc_slider = Slider(
    ax=axslide1,
    label='Kc',
    valmin=1e-3,
    valmax=10,
    valinit=Kc,
)
Kc_slider.on_changed(updateKc)
distForceText = TextBox(axbox5, 'Disturbance', initial='0')
distForceText.on_submit(updateDisturbance)
thetaN = TextBox(axbox2, 'End Time', initial=str(n))
thetaN.on_submit(updateN)
plt.show()

####################################################################
########## Finding the best value of Kc #############
# Best Kc - Kc that minimizes the error at saturation velocity
mError = 1
bestKc = -1
for kc in [0 + 1e-2 * i for i in range(1, 1000)]:
    T, V = rungeKutta(disturbanceForce, theta, (0, vo), n, flag = flag, Kc = kc)
    Vs, Ts = getSaturationVelocity(T, V, tol=1e-5)
    if (Vs != None):
        error = abs((vd(0) - Vs * 5/18)/vd(0))
    if mError > error:
        mError = error
        bestKc = kc
print("Kc Value that minimizes the error:", bestKc)

