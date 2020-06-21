from qutip import *
import numpy as np
import matplotlib.pyplot as plt


# 列表形式或参量形式均可
S1=tensor(basis(2,0),basis(2,0))
S1=tensor([basis(2,0),basis(2,0)])
print(S1)

print(identity(3))
print(qeye(3))
# 上面两个好像差不多
print(tensor(sigmaz(),identity(2)))
print(tensor(identity(2),sigmaz()))


# 用张量来构建复合哈密顿量
# 两个耦合比特
H=tensor(sigmaz(),identity(2))+tensor(identity(2),sigmaz())+0.05*tensor(sigmax(),sigmax())
print(H)
# 三个耦合比特
H=tensor(sigmaz(),identity(2),identity(2))+\
    tensor(identity(2),sigmaz(),identity(2))+\
    tensor(identity(2),identity(2),sigmaz())+\
    0.5*tensor(sigmax(),sigmax(),identity(2))+\
    0.25*tensor(identity(2),sigmax(),sigmax())
print(H)

# 两比特和一个腔耦合(JC哈密顿)
N=10
omega_a=1.0
omega_c=1.25
g=0.05
a=tensor(identity(2),destroy(N))
sm=tensor(destroy(2),identity(N))
sz=tensor(sigmaz(),identity(N))
H=0.5*omega_a*sz+omega_c*a.dag()*a+g*(a.dag()*sm+a*sm.dag())
# Here N is the number of Fock states included in the cavity mode.
print(H)

# partial trace局部跟踪是通过平均（跟踪）消除某些自由度来减小希尔伯特空间的维数的操作。
# 因此，从这个意义上讲，它是张量积的反函数。
psi=tensor(basis(2,0),basis(2,1))
# 选择追踪0态或者1态
print(psi.ptrace(0))
print(psi.ptrace(1))

psi=tensor((basis(2,0)+basis(2,1)).unit(),basis(2,0))
print(psi)
print(psi.ptrace(0))

rho=tensor(ket2dm((basis(2,0)+basis(2,1)).unit()),fock_dm(2,0))
print(ket2dm((basis(2,0)+basis(2,1)).unit()))
print(fock_dm(2,0))
print(rho)
print(rho.ptrace(0))

# 3.4.4 Superoperators and Tensor Manipulations
A=qeye(2)
print(A)
B=qeye(3)
print(to_super(tensor(A,B)))
print(tensor(to_super(A),to_super(B)))
# to_super()效果将X放在反对角线上
# 对于张量tensor，里面的参量类型有oper和super，效果不同
# 如果是super类型就可以通过supertensor来达到一样效果,to_super(tensor(A, B)) == super_tensor(to_super(A), to_super(B))


# 可以通过composite来转换