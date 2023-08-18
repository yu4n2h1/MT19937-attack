from sage.all import Matrix,Zmod,vector
from tqdm import tqdm

def read_Matrix(path):
    """读取目录下的两个矩阵"""
    with open(path) as f:
        data = [[int(y) for y in x.strip().split(',')] for x in f.readlines()]
        print(len(data))
        M = Matrix(Zmod(2),624*32,624*32)
        print("矩阵读取开始")
        for i in tqdm(range(624*32)):
            for y in data[i]:
                M[i,y] = 1
    return M


def state_to_vector(state):
    origin_state = vector(Zmod(2),624*32)

    for i in (range(624)):
        for j in range(32):
            if (state[i] & (1<<j)) > 0:
                origin_state[i*32+j] = 1
    
    return origin_state

if __name__ == "__main__":
    E = read_Matrix("simulate/extract_Matix.txt")
    print(E[0])