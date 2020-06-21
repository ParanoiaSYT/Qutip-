from qutip import *
import numpy as np
import matplotlib.pyplot as plt


vac=basis(5,0)
print(vac)

a=destroy(5)
print(a)
print(a*vac,a.dag()*vac)

c=create(5)
print(c)


print('==========================')
# unit()进行归一化
print((c**2*vac))
print((c**2*vac).unit())
print(c*a*(c**2*vac).unit())
print((c**2*vac).unit())


# 粒子数算子
n=num(5)
print(n)
ket=basis(5,2)
print(n*ket)

ket=(basis(5,0)+basis(5,1)).unit()
print(ket)
print(n*ket)

# 位移算子、挤压算子？？
d=displace(5,1j)
s=squeeze(5,0.25+0.25j)
print(d*vac)
print(d*s*vac)


print('==========================')
# 两能级系统
spin=basis(2,0)
vac=basis(2,0)
print(vac)
c=create(2)
print(c*vac)
# sigmap就是二乘二0100矩阵
print(sigmap()*spin)
print(sigmaz()*spin)

spin2=basis(2,1)
print(spin2)


print('==========================')
# 平均值
vac=basis(5,0)
two=basis(5,2)
c=create(5)
N=num(5)
print(expect(N,vac))
print(expect(N,two))
# 自旋平均值
up=basis(2,0)
print(up)
down=basis(2,1)
print(expect(sigmaz(),up))
print(expect(sigmaz(),down))
# 自旋张量期望值
spin1=basis(2,0)
spin2=basis(2,1)
two_spin=tensor(spin1,spin2)
print(two_spin)
print(qeye(2))   # 二维1向量
sz1=tensor(sigmaz(),qeye(2))
sz2=tensor(qeye(2),sigmaz())
print(expect(sz1,two_spin))
print(expect(sz2,two_spin))

print('==========================')
# 3.3.6 超级运算符和矢量化运算符
psi=basis(2,0)
# 密度矩阵形式
rho=ket2dm(psi)
print(rho)
# 密度矩阵变为密度矩阵矢量式
vec_rho=operator_to_vector(rho)
print(vec_rho)
# 密度矩阵矢量式变为密度矩阵
rho2=vector_to_operator(vec_rho)
print(rho2)
# norm()应该是标准化
print((rho-rho2).norm())

A=Qobj(np.arange(4))
print(A)
A=Qobj(np.arange(4).reshape((2,2)))
# reshape里面赋值是一个元组
print(A)
print(operator_to_vector(A))
print('==========================')


# 希尔伯特空间的线性映射可以用矩阵表示，称为超算子
X=sigmax()
S=spre(X)*spost(X.dag())
print(spre(X))
# spre()效果将
print(spost(X.dag()))
print(S)
S2=to_super(X)
print(S2)
print('==========================')
# to_super()效果将X放在反对角线上


H=10*sigmaz()
c1=destroy(2)
L=liouvillian(H,[c1])
print(L)
