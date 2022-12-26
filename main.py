import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation
import webbrowser

# Animation frames per second
FPS = 30

# Number of bobs
while True:
    try:
        n = int(input('How many bobs does the pendulum have?: '))
        if n > 0:
            break
        else:
            print('Input must be strictly positive.')
    except:
        print('Input must be an integer.')

# Get initial conditions
def request_float(message: str, positive: bool) -> float:
    while True:
        try:
            f = float(input(message))
            if not(positive):
                return f
            if f > 0:
                return f
            print('Input must be strictly positive.')
        except:
            print('Input must be a float.')

print("Input each bob's parameters.")
theta_0 = []
theta_d_0 = []
length = []
mass = []
for i in range(n):
    theta_0.append(request_float(f'Initial angle wrt vertical of pendulum {i+1} (rad): ', False))
    theta_d_0.append(request_float(f'Initial angular velocity of pendulum {i+1} (rad/s): ', False))
    length.append(request_float(f'Length of pendulum {i+1} (m): ', True))
    mass.append(request_float(f'Mass of pendulum {i+1} (kg): ', True))

gravity = 9.8
t0 = request_float('Initial time (s): ', False)
tf = request_float('Final time (s): ', False)
T = np.abs(tf - t0)
time = np.linspace(t0, tf, int(FPS * T) + 1)

# Second derivative of theta
def theta_dd(
    g: float, 
    m: list[float], 
    l: list[float], 
    theta: list[float], 
    theta_d: list[float]
    ) -> np.array:
    rows_A = []
    for k in range(n):
        rows_A.append([l[j] * np.cos(theta[k] - theta[j]) * np.sum(m[np.max([j, k]):]) for j in range(n)])
    A = np.array(rows_A)
    rows_B = []
    for k in range(n):
        rows_B.append([l[j] * np.sin(theta[k] - theta[j]) * np.sum(m[np.max([j, k]):]) for j in range(n)])
    B = np.array(rows_B)
    c = np.array([[np.sin(theta[i]) * np.sum(m[i:])] for i in range(n)])
    theta_d_2 = np.array([[theta_d[i] ** 2] for i in range(n)])
    ans = - np.linalg.inv(A) @ (g * c + B @ theta_d_2)
    return ans.ravel()

# Define new vector quantity S = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
S0 = (*theta_0, *theta_d_0)

def dSdt(S: list[float], t: float, g: float, m: list[float], l: list[float]) -> list[float]:
    theta = S[0:n]
    theta_d = S[n:2*n] 
    return [
        *theta_d,
        *theta_dd(g, m, l, theta, theta_d)
    ]

# Solution
sol = odeint(dSdt, y0=S0, t=time, args=(gravity, mass, length))

length = np.reshape(length, (n, 1))
x = np.zeros((n + 1, len(time))) # Start with row of zeros representing anchor point
y = np.zeros((n + 1, len(time)))
for i in range(n):
    i += 1
    x[i] = np.sum(length[0:i] * np.sin(sol.T[0:i]), axis=0)
    y[i] = - np.sum(length[0:i] * np.cos(sol.T[0:i]), axis=0)

# Animation
def animate(t):
    ln.set_data(x.T[t], y.T[t])

fig = plt.figure()
ln, = plt.plot([], [], 'o-')
lim = np.sum(length) * 1.1
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.title(f'{n}-pendulum')
ani = animation.FuncAnimation(fig, animate, frames=len(time), interval=50)
ani.save('pendulum.gif', writer='pillow', fps=FPS)
webbrowser.open('pendulum.gif')