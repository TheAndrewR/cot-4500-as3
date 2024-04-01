def euler(func, ip, range1, range2, iterations):
    h = (range2 - range1) / iterations
       
    t = range1
    y = ip
       
    for _ in range(iterations):
        y = y + h * func(t, y)
        t = t + h
    
    return y

def func(t, y):
    return t - y**2

ip = 1

range1 = 0
range2 = 2
iterations = 10

result = euler(func, ip, range1, range2, iterations)
print("Q1 result:", result)

def rungekutta(func, ip, range1, range2, iterations):   
    h = (range2 - range1) / iterations
       
    t = range1
    y = ip
      
    for _ in range(iterations):
        k1 = h * func(t, y)
        k2 = h * func(t + h/2, y + k1/2)
        k3 = h * func(t + h/2, y + k2/2)
        k4 = h * func(t + h, y + k3)
        
        y += (k1 + 2*k2 + 2*k3 + k4) / 6
        t += h
    
    return y

def func(t, y):
    return t - y**2

ip = 1

range1 = 0
range2 = 2
iterations = 10

result = rungekutta(func, ip, range1, range2, iterations)
print("Q2 result:", result)

def gaussianelimination(matrix):
    n = len(matrix)

    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(matrix[k][i]) > abs(matrix[max_row][i]):
                max_row = k
        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]

        for k in range(i + 1, n):
            factor = matrix[k][i] / matrix[i][i]
            for j in range(i, n + 1):
                matrix[k][j] -= factor * matrix[i][j]

    x = [0] * n
    for i in range(n - 1, -1, -1):
        x[i] = matrix[i][n] / matrix[i][i]
        for k in range(i - 1, -1, -1):
            matrix[k][n] -= matrix[k][i] * x[i]

    return x

augmentedMatrix = [
    [2, -1, 1, 6],
    [1, 3, 1, 0],
    [-1, 5, 4, -3]
]

answer = gaussianelimination(augmentedMatrix)

print("Q3 result:", answer)

def LUFactor(matrix):
    n = len(matrix)
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]

    for j in range(n):
        L[j][j] = 1.0
        for i in range(j + 1):
            s1 = sum(U[k][j] * L[i][k] for k in range(i))
            U[i][j] = matrix[i][j] - s1
        for i in range(j, n):
            s2 = sum(U[k][j] * L[i][k] for k in range(j))
            L[i][j] = (matrix[i][j] - s2) / U[j][j]

    return L, U

def determinant(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    sign = 1
    det = 0
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in matrix[1:]]
        det += sign * matrix[0][j] * determinant(minor)
        sign *= -1
    return det

matrix = [
    [1, 1, 0, 3],
    [2, 1, -1, 1],
    [3, -1, -1, 2],
    [-1, 2, 3, -1]
]

L, U = LUFactor(matrix)

det = determinant(matrix)
print("Determinant:", det)

print("L matrix:")
for row in L:
    print(row)

print("U matrix:")
for row in U:
    print(row)
    

def DD(matrix):
    n = len(matrix)
    for i in range(n):
        diagonal = abs(matrix[i][i])
        rs = sum(abs(matrix[i][j]) for j in range(n) if j != i)
        if diagonal <= rs:
            return False
    return True

matrix = [
    [9, 0, 5, 2, 1],
    [3, 9, 1, 2, 1],
    [0, 1, 7, 2, 3],
    [4, 2, 3, 12, 2],
    [3, 2, 4, 0, 8]
]

result = DD(matrix)
if result:
    print("True")
else:
    print("False")
    
    
    
def PD(matrix):
    n = len(matrix)
    for i in range(1, n+1):
        minor = [[matrix[row][col] for col in range(i)] for row in range(i)]
        determinant = calculateD(minor)
        if determinant <= 0:
            return False
    return True

def calculateD(matrix):
    n = len(matrix)
    if n == 1:
        return matrix[0][0]
    sign = 1
    det = 0
    for j in range(n):
        minor = [row[:j] + row[j+1:] for row in matrix[1:]]
        det += sign * matrix[0][j] * calculateD(minor)
        sign *= -1
    return det

matrix = [
    [2, 2, 1],
    [2, 3, 0],
    [1, 0, 2]
]

result = PD(matrix)
if result:
    print("True")
else:
    print("False")
