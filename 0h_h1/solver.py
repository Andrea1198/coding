import numpy as np
# from main import solve

def fill():
    return 0 not in grid

def lines():
    for i in range(N):
        count = 0
        for j in range(N):
            if grid[i,j] == 1: count += 1
            elif grid[i,j] == 2: count -= 1
        if count != 0: return False
        count = 0
        for j in range(N):
            if grid[j,i] == 1: count += 1
            elif grid[j,i] == 2: count -= 1
        if count != 0: return False
    return True

def same(lin1, lin2):
    for i in range(N):
        if lin1[i] != lin2[i]: return False
    return True
    
def same_line():
    for i in range(N-1):
        for j in range(i+1,N):
            if same(grid[j,:], grid[i,:]): return False
            if same(grid[:,j], grid[:,i]): return False
    return True

def threes():
    for i in range(N):
        count = 0
        for j in range(N):
            if grid[i,j] == 1: count += 1
            elif grid[i,j] == 2: count = 0
            if count == 3: return False
        count = 0
        for j in range(N):
            if grid[j,i] == 1: count += 1
            elif grid[j,i] == 2: count = 0
            if count == 3: return False
        count = 0
        for j in range(N):
            if grid[i,j] == 2: count += 1
            elif grid[i,j] == 1: count = 0
            if count == 3: return False
        count = 0
        for j in range(N):
            if grid[j,i] == 2: count += 1
            elif grid[j,i] == 1: count = 0
            if count == 3: return False
    return True
    
def check_win():
    if fill() and lines() and same_line() and threes():
        return True
    # elif fill():
    #     if not lines(): char = "lines"
    #     elif not same_line(): char = "same line"
    #     elif not threes(): char = "threes"
    #     print("There is some error", char)
    return False

def solve(i = 0, depth = 0):
    if i == len(indexes[0]):
        if check_win():
            return True
        else:
            return False
    else:
        grid[indexes[0][i], indexes[1][i]] = 1
        if solve(i+1, depth+1):
            return True
        grid[indexes[0][i], indexes[1][i]] = 2
        if solve(i+1, depth):
            return True
        return False

def call_solve(grid_,N_):
    global grid, indexes,N 
    N=N_
    grid = grid_
    indexes = np.where(grid==0)
    solve()

def main():
    global grid, indexes, N
    N = 4
    grid = np.array([
        [1,2,0,1],
        [2,2,1,0],
        [1,1,2,0],
        [2,1,1,0]    
    ])

    indexes = np.where(grid == 0)
    for i,j in zip(indexes[0], indexes[1]):
        grid[i,j] = 1 

    print(solve())
    print(grid)

if __name__ == "__main__":
    main()