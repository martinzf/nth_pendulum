import numpy as np
from scipy.integrate import odeint
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use('fast')

# Animation frames per second
DT = .03

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
    d_theta0 = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    print('Length of each pendulum (m): ')
    l = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    print('Mass of each pendulum (kg): ')
    m = np.array([request(float, f'{i + 1}: ', False) for i in range(n)])
    # Gravity
    g = request(float, 'Gravity (m/s^2): ', True)
    # Duration
    T = request(float, 'Duration (s): ', False)
    t = np.arange(0, T, DT)
    return n, theta_0, d_theta0, l, m, g, t

# ODE
def simple_pendulum(y, _, l, g):
    theta = y[:len(y)//2]
    d_theta = y[len(y)//2:]
    dd_theta = - g / l * np.sin(theta)
    return np.concatenate((d_theta, dd_theta))

def n_pendulum(y, _, l, g, main_diag, off_diag, mu):
    theta = y[:len(y)//2]
    d_theta = y[len(y)//2:]
    phase = theta[:-1] - theta[1:]
    Zinv = sp.diags((main_diag, np.exp(1j * phase) * off_diag, np.exp(-1j * phase) * off_diag), 
                           offsets=(0, 1, -1))
    RZinv_inv = sp.linalg.inv(np.real(Zinv))
    A = np.imag(Zinv) * RZinv_inv * (l * d_theta ** 2)
    B = - g * (np.real(Zinv) + np.imag(Zinv) * RZinv_inv * np.imag(Zinv)) * (mu * np.sin(theta))
    dd_theta = (A + B) / l
    return np.concatenate((d_theta, dd_theta))

# Return to cartesian coordinates
def theta2xy(l: np.array, theta: np.array) -> tuple[np.array]:
    l = np.reshape(l, (n, 1))
    x = np.zeros((n + 1, len(t))) # Start with row of zeros representing anchor point
    y = np.zeros((n + 1, len(t)))
    s = l * np.sin(theta)
    c = l * np.cos(theta)
    for i in range(n):
        i += 1
        x[i] = np.sum(s[:i], axis=0)
        y[i] = - np.sum(c[:i], axis=0)
    return x, y

# Animation
def animate(i: int):
    ln.set_data(x.T[i], y.T[i])
    return ln,

if __name__ == '__main__':
    # Data
    n, theta0, d_theta0, l, m, g, t = get_data()
    # Vector quantity y = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
    y0 = np.concatenate((theta0, d_theta0))
    # Solution
    main_diag = 1/m + np.concatenate(([0], 1/m[:-1]))
    off_diag = - 1/m[:-1]
    mu = np.flip(np.cumsum(np.flip(m)))
    if n == 1:
        sol = odeint(simple_pendulum, y0, t, args=(l, g))
    else:
        sol = odeint(n_pendulum, y0, t, args=(l, g, main_diag, off_diag, mu))
    # Return to cartesian coordinates
    x, y = theta2xy(l, sol.T[:n])
    # Animation
    fig = plt.figure()
    ln, = plt.plot([], [], 'o-')
    lim = np.sum(l) * 1.1
    plt.xlim(- lim, lim)
    plt.ylim(- lim, lim)
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title(f'{n}-Pendulum')
    ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=DT*1e3, blit=True)
    plt.show()