import numpy as np
from tqdm import tqdm

class Model:
    def __init__(self, T, L):
        self._T = T
        self._lattice = np.ones((L,L), dtype=np.int)
        self._L = L

    def _calc_de(self, row, col):
        left_x = (col-1)%self._L
        right_x = (col+1)%self._L
        up_y = (row-1)%self._L
        down_y = (row+1)%self._L
        left = self._lattice[row,left_x]
        right = self._lattice[row, right_x]
        up = self._lattice[up_y, col]
        down = self._lattice[down_y, col]
        s = left + right + up + down
        de = 2 * self._lattice[row, col] * s
        return de

    def change(self, row, col):
        self._lattice[row, col] *= -1

    def trial(self, row, col):
        de = self._calc_de(row, col)
        if de < 0:
            self.change(row, col)
        else:
            if np.random.rand() <= np.exp(-de/self._T):
                self.change(row, col)

    def mcs(self):
        for row in range(self._lattice.shape[0]):
            for col in range(self._lattice.shape[1]):
                self.trial(row, col)

    def cal_m(self):
        return np.mean(self._lattice)

def main():
    m = Model(1.7, 10)
    for i in tqdm(range(30000)):
        m.mcs()
    with open('flips.csv', "w") as f:
        f.write("k,m\n")
        for i in tqdm(range(2200000)):
            m.mcs()
            if i%100 == 0:
                ma = m.cal_m()
                f.write(f"{i},{ma}\n")

    # m._trial(0,0)
    # print(m._lattice)

if __name__ == '__main__':
    main()
