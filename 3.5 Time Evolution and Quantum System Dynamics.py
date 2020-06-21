from qutip import *
import qutip.solver
import numpy as np
import matplotlib.pyplot as plt


result=qutip.solver.Result()
print(result)


# 对于一个自旋1/2系统，幺正演变
H=2*np.pi*0.1*sigmax()
psi0=basis(2,0)
times=np.linspace(0.0,10.0,100)
# 从0到10,取100个点
print(times)
# unitary evolution幺正变换
result=sesolve(H,psi0,times,[sigmaz(),sigmay()])
# 给哈密顿量、初态和时间求解
print(result)
# result.except返回一个列表，里面的值是result第四个参数（列表中各项的平均值）
fig,ax=plt.subplots()
ax.plot(result.times,result.expect[0])
ax.plot(result.times,result.expect[1])
ax.set_xlabel('Time')
ax.set_ylabel('Exceptation values')
ax.legend(("Sigma-Z","Sigma-Y"))
# plt.show()


# 仍以自旋1/2为例
# 对于耗散过程系统，加入耗散项(第四项),第五项是要求的平均
H=2*np.pi*0.1*sigmax()
psi0=basis(2,0)
times=np.linspace(0.0,10.0,100)
result=mesolve(H,psi0,times,[np.sqrt(0.05)*sigmax()],[sigmaz(),sigmay()])
fig,ax=plt.subplots()
ax.plot(result.times,result.expect[0])
ax.plot(result.times,result.expect[1])
ax.set_xlabel('Time')
ax.set_ylabel('Exceptation values')
ax.legend(("Sigma-Z","Sigma-Y"))
# plt.show()


# 腔和原子相互作用体系
# 一个两能级原子与一个泄露单模腔通过偶极相互作用耦合
times=np.linspace(0.0,10.0,200)
psi0=tensor(fock(2,0),fock(10,5))
print(fock(10,5))
a=tensor(qeye(2),destroy(10))
sm=tensor(destroy(2),qeye(10))
H=2*np.pi*a.dag()*a+2*np.pi*sm.dag()*sm+\
  2*np.pi*0.25*(sm*a.dag()+sm.dag()*a)
result=mesolve(H,psi0,times,[np.sqrt(0.1)*a],[a.dag()*a,sm.dag()*sm])
# 第四项应该就是腔的耗散
plt.figure()
plt.plot(times,result.expect[0])
plt.plot(times,result.expect[1])
plt.xlabel('Time')
plt.ylabel('Expectation values')
plt.legend(("cavity","atom excitation probability"))


# 3.5.3 蒙特卡罗解法(出现报错)
# 两能级系统与一个泄露腔耦合
times=np.linspace(0.0,10.0,200)
psi0=tensor(fock(2,0),fock(10,5))
a=tensor(qeye(2),destroy(10))
sm=tensor(destroy(2),qeye(10))
H=2*np.pi*a.dag()*a+2*np.pi*sm.dag()*sm+\
  2*np.pi*0.25*(sm*a.dag()+sm.dag()*a)
# data=mcsolve(H,psi0,times,[np.sqrt(0.1)*a],[a.dag()*a,sm.dag()*sm])
# mcsolve出错
data=mesolve(H,psi0,times,[np.sqrt(0.1)*a],[a.dag()*a,sm.dag()*sm])
plt.figure()
plt.plot(times,data.expect[0],times,data.expect[1])
plt.title('Monte Carlo time evolution')
plt.xlabel('Time')
plt.ylabel('Expectation values')
plt.legend(("cavity photon number","atom excitation probabability"))
# plt.show()


# 3.5.6 时间演化哈密顿量
def H1_coeff(t, args):
    A = args['A']
    sig = args['sigma']
    return A * np.exp(-(t / sig) ** 2)

ustate=basis(3,0)
excited=basis(3,1)
ground=basis(3,2)
N=2
sigma_ge=tensor(qeye(N),ground*excited.dag())
sigma_ue=tensor(qeye(N),ustate*excited.dag())
a=tensor(destroy(N),qeye(3))
ada=tensor(num(N),qeye(3))
c_ops=[]    # Build collapse operators
kappa=1.5   # Cavity decay rate
c_ops.append(np.sqrt(kappa)*a)
gamma=6
c_ops.append(np.sqrt(5*gamma/9)*sigma_ue)
c_ops.append(np.sqrt(4*gamma/9) * sigma_ge) # 4/9 e->g
t = np.linspace(-15, 15, 100) # Define time vector
psi0 = tensor(basis(N, 0), ustate) # Define initial state
state_GG = tensor(basis(N, 1), ground) # Define states onto which to ˓→project
sigma_GG = state_GG * state_GG.dag()
state_UU = tensor(basis(N, 0), ustate)
sigma_UU = state_UU * state_UU.dag()
g=5
H0 = -g * (sigma_ge.dag() * a + a.dag() * sigma_ge) # time-independent ˓→term
H1 = (sigma_ue.dag() + sigma_ue) # time-dependent term
H = [H0,[H1,H1_coeff]]
output = mesolve(H, psi0, t, c_ops, [ada, sigma_UU, sigma_GG])
