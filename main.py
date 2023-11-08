import numpy as np
from scipy.integrate import odeint
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import animation
plt.style.use('fast')
import gui

# Animation frames per second
DT = .03

# ODE
def simple_pendulum(y, _, l, g):
    theta = y[:len(y)//2]
    d_theta = y[len(y)//2:]
    dd_theta = - g / l * np.sin(theta)
    return np.concatenate((d_theta, dd_theta))

def n_pendulum(y, _, n, l, g, Zinv, off_diag, mu):
    theta = y[:n]
    d_theta = y[n:]
    Zinv.setdiag(np.exp(- 1j * np.diff(theta)) * off_diag, 1)
    Zinv.setdiag(np.exp(1j * np.diff(theta)) * off_diag, - 1)
    RZinv_inv = sp.linalg.factorized(np.real(Zinv) + sp.csc_matrix((n, n)))
    dd_theta = 1 / l * (np.imag(Zinv) * RZinv_inv(l * d_theta ** 2) -
               g * (np.real(Zinv) * (mu * np.sin(theta)) +
               np.imag(Zinv) * RZinv_inv(np.imag(Zinv) * (mu * np.sin(theta)))))
    return np.concatenate((d_theta, dd_theta))

if __name__ == '__main__':
    # Data
    n, theta0, d_theta0, l, m, g, dur = gui.setup()
    t = np.arange(0, dur, DT)
    # Vector quantity y = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
    y0 = np.concatenate((theta0, d_theta0))
    # Z^-1 matrix
    off_diag = - 1/m[:-1]
    Zinv = sp.diags((np.exp(- 1j * np.diff(theta0)) * off_diag, 
                     1/m + np.concatenate(([0], 1/m[:-1])),
                     np.exp(1j * np.diff(theta0)) * off_diag), 
                     offsets=(1, 0, -1), format='csc')
    # Solution
    mu = np.flip(np.cumsum(np.flip(m)))
    if n == 1:
        y = odeint(simple_pendulum, y0, t, args=(l, g))
    else:
        y = odeint(n_pendulum, y0, t, args=(n, l, g, Zinv, off_diag, mu))
    # Conversion to cartesian coordinates
    x = np.cumsum(l[:, np.newaxis] * np.sin(y.T[:n]), axis=0)
    y = - np.cumsum(l[:, np.newaxis] * np.cos(y.T[:n]), axis=0)
    x = np.vstack((np.zeros(len(t)), x)) # Start with row of zeros representing anchor point
    y = np.vstack((np.zeros(len(t)), y))
    # Animation
    fig = plt.figure(figsize=(8, 8))
    ln, = plt.plot([], [], 'o-')
    lim = np.sum(l) * 1.1
    plt.xlim(- lim, lim)
    plt.ylim(- lim, lim)
    plt.xlabel('X (m)')
    plt.ylabel('Y (m)')
    plt.title(f'{n}-Pendulum')
    def animate(i):
        ln.set_data(x.T[i], y.T[i])
        return ln,
    ani = animation.FuncAnimation(fig, animate, frames=len(t), interval=DT*1e3, blit=True)
    plt.show()