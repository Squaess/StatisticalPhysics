from numba import njit
import numpy as np
import csv

@njit(parallel=True)
def main(N, L):
    for i in range(1, 51):
        r2 = 0
        for _ in range(10):
            r2 += diffusion(i, N, L)
        r2 = r2/10
        yield (i, r2/(4*i))



@njit(parallel=True)
def diffusion(number_of_steps, N, L):
    A = np.ones((L,L), dtype=np.bool_)
    x = np.zeros(N, dtype=np.int64)
    y = np.zeros(N, dtype=np.int64)

    dx = np.zeros(N, dtype=np.int64)
    dy = np.zeros(N, dtype=np.int64)

    i = N-1
    while i >= 0:
        xt = np.random.randint(L)
        yt = np.random.randint(L)
        if A[xt, yt]:
            A[xt, yt] = False
            x[i], y[i] = xt, yt
            i = i-1

    for _ in range(number_of_steps):
        for i in range(N):
            direction = np.random.randint(4)
            xt, yt = x[i], y[i]
            tdx, tdy = 0, 0
            if direction == 0:
                xt = (x[i] + 1) % L 
                tdx = 1
            elif direction == 1:
                xt = (x[i] - 1) % L 
                tdx = -1
            elif direction == 2:
                yt = (y[i] + 1) % L 
                tdy = 1
            else:
                yt = (y[i] - 1) % L 
                tdy = -1
            if A[xt, yt]:
                A[xt,yt] = False
                A[x[i], y[i]] = True
                dx[i] += tdx
                dy[i] += tdy
                x[i], y[i] = xt, yt

    R2 = np.mean(np.power(dx,2) + np.power(dy,2))
    return R2

if __name__ == "__main__":
    N = 40
    L = 20
    c = main(N, L)

    with open("result2.csv", "w", newline='') as f:
        writer = csv.writer(f, delimiter=",")
        writer.writerow(["C", N/L**2])
        for i,d in c:
            writer.writerow([i, d])