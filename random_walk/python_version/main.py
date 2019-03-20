from numba import njit
import numpy as np
import time

@njit(parallel=True)
def run_numba(steps_array):

    ret_array = np.zeros((steps_array.shape[0],2), dtype=np.double)
    # ilosc pijakow
    K = np.int(30000)
    x = 0
    xNa = 0.
    xN2a = 0.

    for i in range(steps_array.shape[0]):
        N = np.int(steps_array[i])
        xNa = 0
        xN2a = 0
        for j in range(K):
            x = 0
            for step in range(N):
                if np.random.rand() < 0.5:
                    x -= 1
                else:
                    x += 1
            xNa += x
            xN2a += np.square(x)
        s = np.square(np.divide(xNa, K))
        o = np.divide(xN2a, K)
        u = np.sqrt(o - s)
        ret_array[i,0] = N
        ret_array[i,1] = u

    return ret_array

steps_array = np.logspace(1.0, 4.0, 10)
# steps_array = np.arange(1, 1001)
print(steps_array)
start = time.time()
ret = run_numba(steps_array)
end = time.time()
print("Elapsed (with compilation) = %s" % (end - start))
start = time.time()
ret = run_numba(steps_array)
end = time.time()
print("Elapsed (after compilation) = %s" % (end - start))
np.savetxt("data.csv", ret, delimiter=",")