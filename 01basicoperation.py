from qutip import *
import numpy as np
import matplotlib.pyplot as plt


# 0态bra
print(Qobj())

# bra
x=np.array([[1,2,3,4,5]])
print(x)
print(Qobj(x))
# ket
y=np.array([[1],[2],[3],[4],[5]])
print(y)
print(Qobj(y))

