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

# 随机的4乘4的算子
r=np.random.rand(4,4)
print(Qobj(r))

# 0~4五个态，处于3态的state
print(basis(5,3))

# 相干态，其中0.5-0.5j是本征值
print(coherent(5,0.5-0.5j))

# 四阶湮灭算子
print(destroy(4))

print(sigmaz())

# Higher spin operators
print(jmat(1/2.0,'z'))


print('===========================')
# 可以显示属性(data,dims,shape,isherm,type)
q=destroy(4)
print(q.dims)
print(q.shape)
print(q.type)


print('===========================')
# 运算规则和矩阵一样
q=destroy(4)
print(q+5)
# print(q**2)
print(q**3)

x=sigmax()
print(x*x)
print(x/np.sqrt(2))


print('===========================')
print(basis(5,3))
# 共轭
print(basis(5,3).dag())

# 密度矩阵
print(coherent_dm(5,1))

# 对角项
print(coherent_dm(5,1).diag())

# 完整复数形式
print(coherent_dm(5,1).full())

# L2 norm就是欧几里德距离
print(coherent_dm(5,1).norm())

