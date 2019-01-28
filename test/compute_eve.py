# Compute 5 eigenvectors corresponding to 5 greatest eigenvalues of given matrix
# Format is .mtx (Market Matrix)
# Eigenvectors are writen as 5 line in file <matrix-path>.eve

from scipy.io import mmread,mmwrite
import numpy as np

filename = "matrices/bcspwr01"

A = mmread(filename + ".mtx")

# eig[0] are eigenvalues
# aig[1] are eigenvectors (a.k.a eve)
eig = np.linalg.eig(A.todense())

# Sort both by eigenvalue's absolute value
idx = abs(eig[0]).argsort()[::-1]   
eigenValues = eig[0][idx]
eigenVectors = eig[1][:,idx]

# Write in a file 5 first eves
with open(filename + '.eve', 'w') as f:
    for line in eigenVectors[:5]:
        np.savetxt(f, line, fmt='%f')
