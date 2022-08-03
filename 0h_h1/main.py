import pygame 
import numpy as np 
from numba import njit
from solver import call_solve 

N = 4
grid = np.zeros((N, N))
YELLOW = (255, 255, 0)
BLUE = (0, 0, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
width = 600 
height = 600  
dx = width // N
dy = height // N
screen  = pygame.display.set_mode((width, height))

def create_game():
    tmp = grid
    i = np.random.randint(0,N)
    j = np.random.randint(0,N)
    while call_solve(tmp,N):
        tmp = grid
        i = np.random.randint(0,N)
        j = np.random.randint(0,N)
        while grid[i,j] != 0:
            i = np.random.randint(0,N)
            j = np.random.randint(0,N)
        
        grid[i,j] = np.random.randint(1,3)
    grid[i,j] = 0
    pass
def eliminate_max():
    return 

def update_grid(i,j):
    grid[i,j] += 1
    grid[i,j] %= 3


def count_wins():
    pass

def show_grid():
    for i in range(N):
        for j in range(N):
            if grid[i,j] == 1:
                pygame.draw.rect(screen, YELLOW, (i*dx, j*dy, dx, dy))
            elif grid[i,j] == 2:
                pygame.draw.rect(screen, BLUE, (i*dx, j*dy, dx, dy))
        pygame.draw.line(screen, WHITE, (0, i*dy), (width, i*dy))
        pygame.draw.line(screen, WHITE, (i*dx, 0), (i*dx, height))

# def create_grid():
#     return np.zeros((N, N))

def update_window():
    show_grid()
    pygame.display.update()

def start():
    pygame.init()
    pygame.display.set_caption("0h h1")
    running = True
    create_game()
    eliminate_max() 
    while running:
        screen.fill(BLACK)
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False 
            if event.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                i = pos[0]//dx
                j = pos[1]//dy
                update_grid(i,j)
        # if check_win():
        #     print("You won")
        #     running = False
        update_window()

def main():
    if N%2 != 0: 
        print("N must be even")
        return
    # grid = create_grid()
    start()
if __name__ == "__main__":
    main()
