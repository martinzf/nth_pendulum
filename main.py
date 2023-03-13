import numpy as np
from scipy.integrate import odeint
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
    theta_0 = []
    theta_d_0 = []
    length = []
    mass = []
    for i in range(n):
        theta_0.append(request(float, f'Initial angle wrt vertical of pendulum {i + 1} (rad): ', False))
        theta_d_0.append(request(float, f'Initial angular velocity of pendulum {i + 1} (rad/s): ', False))
        length.append(request(float, f'Length of pendulum {i + 1} (m): ', True))
        mass.append(request(float, f'Mass of pendulum {i + 1} (kg): ', True))
    # Gravity
    gravity = request(float, 'Gravity (m/s^2): ', True)
    # Duration
    T = request(float, 'Duration (s): ', False)
    time = np.linspace(0, T, int(FPS * T) + 1)
    return n, theta_0, theta_d_0, length, mass, gravity, time

# Second derivative of theta
def theta_dd(
    n: int,
    g: float, 
    m: list[float], 
    l: list[float], 
    theta: list[float], 
    theta_d: list[float]
    ) -> np.array:
    A = np.empty((n, n), np.float64)
    B = np.empty((n, n), np.float64)
    for k in range(n):
        for j in range(n):
            A[k, j] = l[j] * np.cos(theta[k] - theta[j]) * np.sum(m[np.max([j, k]):])
            B[k, j] = l[j] * np.sin(theta[k] - theta[j]) * np.sum(m[np.max([j, k]):])
    c = np.array([[np.sin(theta[i]) * np.sum(m[i:])] for i in range(n)])
    theta_d_2 = np.array([[theta_d[i] ** 2] for i in range(n)])
    ans = - np.linalg.inv(A) @ (g * c + B @ theta_d_2)
    return ans.ravel()

# Define new vector quantity S = (theta1, theta2, ..., thetan, d/dt theta1, d/dt theta2, ..., d/dt thetan)
def dSdt(S: list[float], t: float, n:int, g: float, m: list[float], l: list[float]) -> list[float]:
    theta = S[0:n]
    theta_d = S[n:2*n] 
    return [
        *theta_d,
        *theta_dd(n, g, m, l, theta, theta_d)
    ]

# Return to cartesian coordinates
def theta2xy(length, theta):
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
    ani = animation.FuncAnimation(fig, animate, frames=len(time), interval=50)
    ani.save('pendulum.gif', writer='pillow', fps=FPS)
    webbrowser.open('pendulum.gif')