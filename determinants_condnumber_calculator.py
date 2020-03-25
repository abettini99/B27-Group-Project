import numpy as np
from io import StringIO

for i in ['asymmetric_aperiodicroll.txt', 'asymmetric_dutchroll.txt', 'asymmetric_dutchrollYD.txt', 'asymmetric_spiral.txt', 'symmetric_phugoid.txt', 'symmetric_shortperiod.txt']:
    Matrix = []
    matrix_content = open("matrices/" + i)
    lines = matrix_content.readlines()
    for line in lines:
        line = np.genfromtxt(StringIO(line), delimiter=",")
        Matrix.append(line)
    Matrix = np.array(Matrix)
    
    print(i, 'det A matrix = ', np.linalg.det(Matrix))
    print(i, 'Cond no. A matrix = ', np.linalg.cond(Matrix))
    
    