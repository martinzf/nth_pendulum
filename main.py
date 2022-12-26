import numpy as np
import sympy as smp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import animation

while True:
    try:
        n = int(input('How many bobs does the pendulum have?: '))
        if n > 0:
            break
        else:
            print('Number of bobs must be strictly positive')
    except:
        print('Input must be an integer')

# Lengths
symbols = (fr'l_{i}' for i in range(n))
l = smp.symbols(' '.join(symbols))

# Masses
symbols = (fr'm_{i}' for i in range(n))
m = smp.symbols(' '.join(symbols))

# Time & gravity
t, g = smp.symbols('t g')

# Angles
symbols = (fr'\theta_{i}' for i in range(n))
theta = smp.symbols(' '.join(symbols), cls=smp.Function)
theta = [theta[i](t) for i in range(n)]
theta_d = [smp.diff(theta[i], t) for i in range(n)] # 1st derivative

# Matrices
def summation(x: tuple[smp.Symbol], a: int, b: int) -> smp.Symbol:
    # Sums elements in tuple of SymPy symbols
    sum = 0
    for i in np.arange(a, b):
        sum += x[i]
    return sum

rows_A = []
for k in range(n):
    rows_A.append([l[j] * smp.cos(theta[k] - theta[j]) * summation(m, np.max([j, k]), n) for j in range(n)])
A = smp.Matrix(rows_A)

rows_B = []
for k in range(n):
    rows_B.append([l[j] * smp.sin(theta[k] - theta[j]) * summation(m, np.max([j, k]), n) for j in range(n)])
B = smp.Matrix(rows_B)

c = smp.Matrix([[smp.sin(theta[i]) * summation(m, i, n)] for i in range(n)])

theta_d_2 = smp.Matrix([[theta_d[i] ** 2] for i in range(n)])

# Explicit system of ODEs
theta_dd = - (A.inv() * (g * c + B * theta_d_2))
for i, row in enumerate(theta_dd):
    theta_dd[i] = row.simplify()

# Obtain numerical functions
theta_dd_f = smp.lambdify((t, g, *m, *l, *theta, *theta_d), theta_dd.T.tolist()[0])

# Get initial conditions
def request_float():
    pass
theta_0 = [1, 1]
theta_d_0 = [1, 1]
length = np.array([1, 1])
mass = [1, 1]
gravity = 9.8
t0 = 0
tf = 10
time = np.linspace(0, 10, 101)

# Define new vector quantity S = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
S0 = (*theta_0, *theta_d_0)

def dSdt(S: list[float], t: float, g: float, m: list[float], l: list[float]) -> list[float]:
    theta = S[0:n]
    theta_d = S[n:2*n] 
    return [
        *theta_d,
        *theta_dd_f(t, g, *m, *l, *theta, *theta_d)
    ]

# Solution
sol = odeint(dSdt, y0=S0, t=time, args=(gravity, mass, length))

length = np.array([[1], [1]])
x = np.zeros((n, 101))
y = np.zeros((n, 101))
for i in range(n):
    x[i] = np.sum(length[0:i+1] * np.sin(sol.T[0:i+1]))
    y[i] = - np.sum(length[0:i+1] * np.sin(sol.T[0:i+1]))

# Animation
def animate(t):
    ln.set_data(x.T[t], y.T[t])

fig = plt.figure()
ln, = plt.plot([], [], 'o--')
plt.xlim(-4, 4)
plt.ylim(-4, 4)
ani = animation.FuncAnimation(fig, animate, frames=100, interval=50)
ani.save('pendulum.gif', writer='pillow', fps=25)