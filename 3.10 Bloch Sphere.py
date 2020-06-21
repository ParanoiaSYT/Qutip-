from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import mayavi


b=Bloch3d()
print(b)
vec=[0,1,0]
b.add_vectors(vec)
vec2=[1,1,0]
b.add_vectors(vec2)
up=basis(2,0)
b.add_states(up)
pnt=[1/np.sqrt(3),1/np.sqrt(3),1/np.sqrt(3)]
b.add_points(pnt)

b.show()
plt.show()