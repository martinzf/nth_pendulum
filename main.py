import numpy as np
from scipy.integrate import odeint
import scipy.sparse.linalg as la
import numba as nb
import matplotlib.pyplot as plt
from matplotlib import animation
import webbrowser

# Animation frames per second
FPS = 30

# Get user input
def request(type: callable, prompt: str, positive: bool) -> float:
    while True:
        try:
            answer = type(input(prompt))
            if not(positive):
                return answer
            if answer > 0:
                return answer
            print('Input must be strictly positive.')
        except ValueError:
            print(f'Input must be {type}.') 

def get_data() -> tuple[float, list[float], np.array]:
    # Number of bobs
    n = request(int, 'Input the number of links the n-pendulum has: ', True)
    # Initial conditions
    print("Input each pendulum's parameters.")
    print('Initial wrt vertical of each pendulum (rad): ')
    theta_0 = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    print('Initial angular velocity of each pendulum (rad/s): ')
    theta_d_0 = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    print('Length of each pendulum (m): ')
    length = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    print('Mass of each pendulum (kg): ')
    mass = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    # Gravity
    gravity = request(float, 'Gravity (m/s^2): ', True)
    # Duration
    T = request(float, 'Duration (s): ', False)
    time = np.linspace(0, T, int(FPS * T) + 1)
    return n, theta_0, theta_d_0, length, mass, gravity, time

# Second derivative of theta system of equations
@nb.njit('i8, f8, f8[:], f8[:], f8[:], f8[:]', parallel=True, nogil=True, fastmath=True)
def systeq(n, g, m, l, theta, theta_d) -> tuple[np.array]:
    A = np.empty((n, n), dtype=np.float64)
    B = np.empty((n, n), dtype=np.float64)
    c = np.empty(n, dtype=np.float64)
    for k in nb.prange(n):
        c[k] = np.sin(theta[k]) * np.sum(m[k:])
        for j in nb.prange(n):
            i = int(max(j, k))
            A[k, j] = l[j] * np.cos(theta[k] - theta[j]) * np.sum(m[i:])
            B[k, j] = l[j] * np.sin(theta[k] - theta[j]) * np.sum(m[i:])
    return A, - (g * c + B @ theta_d ** 2)

# Define new vector quantity S = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
def dSdt(S: np.array, t: float, n:int, g: float, m: np.array, l: np.array) -> np.array:
    theta = S[0:n]
    theta_d = S[n:2*n] 
    x0 = np.gradient(theta_d) if n > 1 else theta_d
    theta_dd = la.gmres(*systeq(n, g, m, l, theta, theta_d), x0=x0)[0]
    return np.hstack((theta_d, theta_dd))

# Return to cartesian coordinates
def theta2xy(length: np.array, theta: np.array) -> tuple[np.array]:
    length = np.reshape(length, (n, 1))
    x = np.zeros((n + 1, len(time))) # Start with row of zeros representing anchor point
    y = np.zeros((n + 1, len(time)))
    s = length * np.sin(theta)
    c = length * np.cos(theta)
    for i in range(n):
        i += 1
        x[i] = np.sum(s[:i], axis=0)
        y[i] = - np.sum(c[:i], axis=0)
    return x, y

# Animation
def animate(t: int):
    ln.set_data(x.T[t], y.T[t])
    return ln,

if __name__ == '__main__':
    # Data
    n, theta_0, theta_d_0, length, mass, gravity, time = get_data()
    # Vector quantity S = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
    S0 = [*theta_0, *theta_d_0]
    # Solution
    sol = odeint(dSdt, y0=S0, t=time, args=(n, gravity, mass, length))
    # Return to cartesian coordinates
    x, y = theta2xy(length, sol.T[:n])
    # Animation
    fig = plt.figure()
    ln, = plt.plot([], [], 'o-')
    lim = np.sum(length) * 1.1
    plt.xlim(- lim, lim)
    plt.ylim(- lim, lim)
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title(f'{n}-pendulum')
    ani = animation.FuncAnimation(fig, animate, frames=len(time), interval=50, blit=True)
    ani.save('pendulum.gif', writer='pillow', fps=FPS)
    webbrowser.open('pendulum.gif')