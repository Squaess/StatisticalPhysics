import numpy as np
import csv

np.random.seed(442)

def initiliaze_board(board, n):
    i = n
    atoms = []
    while i >= 1:
        x, y = np.random.randint(0, len(board), 2)
        if not board[x, y]:
            board[x, y] = 1
            atoms.append((x, y))
            i -= 1

    return atoms
    
def move_right(atoms, i, board):
    x, y = 0, 0
    if not board[(atoms[i][0]+1) % len(board), atoms[i][1]]:
        board[(atoms[i][0]+1) % len(board), atoms[i][1]] = 1
        board[atoms[i][0], atoms[i][1]] = 0
        atoms[i] = ((atoms[i][0]+1) % len(board), atoms[i][1])
        x = 1
    return x, y

def move_left(atoms, i, board):
    x, y = 0, 0
    if not board[(atoms[i][0]-1) % len(board), atoms[i][1]]:
        board[(atoms[i][0]-1) % len(board), atoms[i][1]] = 1
        board[atoms[i][0], atoms[i][1]] = 0
        atoms[i] = ((atoms[i][0]-1) % len(board), atoms[i][1])
        x = -1
    return x, y

def move_up(atoms, i, board):
    x, y = 0, 0
    if not board[atoms[i][0], (atoms[i][1]+1) % len(board)]:
        board[atoms[i][0], (atoms[i][1]+1) % len(board)] = 1
        board[atoms[i][0], atoms[i][1]] = 0
        atoms[i] = (atoms[i][0], (atoms[i][1]+1) % len(board))
        y = 1
    return x, y

def move_down(atoms, i, board):
    x, y = 0, 0
    if not board[atoms[i][0], (atoms[i][1]-1) % len(board)]:
        board[atoms[i][0], (atoms[i][1]-1) % len(board)] = 1
        board[atoms[i][0], atoms[i][1]] = 0
        atoms[i] = (atoms[i][0], (atoms[i][1]-1) % len(board))
        y = -1
    return x, y

def move(atoms, i, board):
    r = np.random.randint(4)
    if r == 0:
        x, y = move_right(atoms, i, board)
    elif r == 1:
        x, y = move_left(atoms, i, board)
    elif r == 2:
        x, y = move_up(atoms, i, board)
    else:
        x, y = move_down(atoms, i, board)
    return x,y


def diffusion_mcs(number_of_steps, density, L=20):
    n = int(density * L * L)
    board = np.zeros((L,L))
    atoms = initiliaze_board(board, n)
    # print(atoms)
    displacement = [[0, 0] for _ in range(n)]
    for _ in range(number_of_steps):
        for i in range(n):
            x, y = move(atoms, i, board)
            displacement[i][0] += x
            displacement[i][1] += y
    # print(atoms)
    # print(displacement)
    return sum([dx**2 + dy**2 for dx,dy  in displacement])/len(displacement)

if __name__ == "__main__":
    dens = np.arange(0.1, 1, 0.1)
    s = 10
    with open("res_pv2.csv", "w", newline='') as f:
        writer = csv.writer(f, delimiter=',')
        for d in dens:
            print("d = ", d)
            for i in range(1, 51):
                r2 = 0 
                for _ in range(s):
                    r2 += diffusion_mcs(i, d)
                r2 /= s
                coef = r2/(4*i)
                writer.writerow([d, i, coef])
    print("Dupa")