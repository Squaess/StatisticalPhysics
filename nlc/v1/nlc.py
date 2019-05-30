import numpy as np
import math
import matplotlib.pyplot as plt
import random


L = 50

# strength of NLC-NLC orientational interaction
xi = 1

# Nie mam pojęcia, ile to powinno być
k = 1

# To jest krypna wartość wzoru
delta_angle = 180


def to_rad(angle):
    return (angle * math.pi) / 180


def p2(angle):
    return (1.5 * (math.cos(to_rad(angle)) ** 2 - 0.5))


def calc_energy(state, base_value, point):
    """
    state - macierz
    base_value = wartość w punkcie point
    point - (x, y)
    """
    x, y = point
    return - xi * sum([
        p2(base_value - state[(x + 1) % L, y]),
        p2(base_value - state[(x - 1) % L, y]),
        p2(base_value - state[x, (y + 1) % L]),
        p2(base_value - state[x, (y - 1) % L])
    ])


def one_mc_step(matrix, T):
    # for i in range(L):
        # for j in range(L):
    for i, j in get_next_point():
        # print(f'{i}, {j}')
        phi_old = matrix[i, j]
        phi_new = phi_old + (np.random.uniform() - 0.5) * delta_angle
        # normalizacja tych kątów
        if phi_new > 90:
            phi_new -= 90
        if phi_new < -90:
            phi_new += 90
        E_new = calc_energy(matrix, phi_new, (i, j))
        E_old = calc_energy(matrix, phi_old, (i, j))

        if E_new < E_old:
            matrix[i, j] = phi_new
        else:
            if np.random.rand() < np.exp((E_old - E_new) / (k * T)):
                matrix[i, j] = phi_new


def get_next_point():
    for i in range(L*L):
        yield random.randrange(0, L), random.randrange(0, L)


def nlc_experiment(mc_steps, T):

    Q = np.zeros((2, 2))
    steps_to_ignore = 3000
    assert steps_to_ignore < mc_steps

    matrix = np.random.randint(-90, 90, (L, L))

    # print(matrix)

    counter = 0
    for step in range(mc_steps):
        one_mc_step(matrix, T)
        if step > steps_to_ignore and step % 100 == 0:
            print(f'Step: {step}')
            update_q(Q, matrix)
            counter += 1

    # print(matrix)
    # print(Q)
    eigenvalues = np.linalg.eigvals(Q)
    # print(eigenvalues)
    S = max(eigenvalues)
    # return S
    # return (S / (mc_steps - steps_to_ignore) / 100)
    return S / counter



def update_q(Q, matrix):
    for i in range(L):
        for j in range(L):
            phi = matrix[i, j]
            Q[0, 0] = Q[0, 0] + (2.0 * math.cos(to_rad(phi)) ** 2 - 1.0) / L ** 2
            Q[0, 1] = Q[0, 1] + (2.0 * math.cos(to_rad(phi)) * math.sin(to_rad(phi)) ) / L ** 2
            Q[1, 1] = - Q[0, 0]
            Q[1, 0] = Q[0, 1]


temperatures = np.arange(0.1, 3.7, 0.5)

S_result = []
for temperature in temperatures:
    s = nlc_experiment(100000, temperature)
    S_result.append(s)

print(S_result)

plt.plot(temperatures, S_result, '-o')
plt.xlabel("Temperature ")
plt.ylabel("Order parameter S")
plt.savefig('order_parameter_s.png')
