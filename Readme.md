# CFD FDM Practice

> 学习 <https://github.com/surajp92/CFD_Julia>

## 1D Heat Equation

使用有限差分法求方程 $$ \dfrac{\partial u}{\partial t} = \alpha \dfrac{\partial^2 u}{\partial x^2}. \tag{1}$$的近似解.

### FTCS

1. 时间项离散
   将$u^{n+1}_i$在时间步n处向前展开$$ u^{n+1}_i=u^n_i+\dfrac{\partial u^n_i}{\partial t}\Delta t+\dfrac{\partial^2 u^n_i}{\partial t^2}\dfrac{\Delta t^2}{2!}+o(\Delta t^3). \tag{2} $$
   得到前向欧拉(FE)格式$$ \dfrac{\partial u^n_i}{\partial t}=\dfrac{u^{n+1}_i-u^n_i}{\Delta t}+o(\Delta t) \tag{3}$$
2. 扩散项离散
    将$u^n_{i+1}$在空间坐标i+1处向后展开$$ u^n_{i+1}=u^n_i+\dfrac{\partial u^n_i}{\partial x}\Delta x+\dfrac{\partial^2 u^n_i}{\partial x^2}\dfrac{\Delta x^2}{2!}+\dfrac{\partial^3 u^n_i}{\partial x^3}\dfrac{\Delta x^3}{3!}+o(\Delta x^4). \tag{4} $$
    将$u^n_{i-1}$在空间坐标i-1处向前展开$$ u^n_{i-1}=u^n_i-\dfrac{\partial u^n_i}{\partial x}\Delta x+\dfrac{\partial^2 u^n_i}{\partial x^2}\dfrac{\Delta x^2}{2!}-\dfrac{\partial^3 u^n_i}{\partial x^3}\dfrac{\Delta x^3}{3!}+o(\Delta x^4). \tag{5} $$
    (4)和(5)两式相加，得到中心差分格式$$\dfrac{\partial^2 u^n_i}{\partial x^2}=\dfrac{u^n_{i+1}-2u^n_i+u^n_{i-1}}{\Delta x^2}+o(\Delta x^2). \tag{6}$$
3. 将(3)和(6)带入方程(1)，得到离散方程$$ \dfrac{u^{n+1}_i-u^n_i}{\Delta t}=\alpha \dfrac{u^n_{i+1}-2u^n_i+u^n_{i-1}}{\Delta x^2}. \tag{7}$$易知该离散格式在时间上有一阶精度，空间上有二阶精度。

### 3rd Runge-Kutta

略

### Crank-Nicolson

将$u^{n+\frac{1}{2}}_i$在时间步n处向前展开$$ u^{n+\frac{1}{2}}_i=u^n_i+\dfrac{\partial u^n_i}{\partial t}\dfrac{\Delta t}{2}+\dfrac{\partial^2 u^n_i}{\partial t^2}\dfrac{\Delta t^2}{2!\cdot 2^2}+o(\Delta t^3). \tag{8} $$ 得到$$ \dfrac{\partial u^n_i}{\partial t}=2\dfrac{u^{n+\frac{1}{2}}_i-u^n_i}{\Delta t}-\dfrac{\partial^2 u^n_i}{\partial t^2}\dfrac{\Delta t}{2!\cdot 2}+o(\Delta t^2) \tag{9}$$
得到离散方程$$ 2\dfrac{u^{n+\frac{1}{2}}_i-u^n_i}{\Delta t}-\dfrac{\partial^2 u^n_i}{\partial t^2}\dfrac{\Delta t}{2!\cdot 2}+o(\Delta t^2)=\alpha \dfrac{u^n_{i+1}-2u^n_i+u^n_{i-1}}{\Delta x^2}+o(\Delta x^2). \tag{10}$$
将$u^{n+\frac{1}{2}}_i$在时间步n+1处向后展开$$ u^{n+\frac{1}{2}}_i=u^{n+1}_i-\dfrac{\partial u^{n+1}_i}{\partial t}\dfrac{\Delta t}{2}+\dfrac{\partial^2 u^{n+1}_i}{\partial t^2}\dfrac{\Delta t^2}{2!\cdot 2^2}+o(\Delta t^3). \tag{11} $$ 得到$$ \dfrac{\partial u^n_i}{\partial t}=2\dfrac{u^{n+\frac{1}{2}}_i-u^n_i}{\Delta t}+\dfrac{\partial^2 u^n_i}{\partial t^2}\dfrac{\Delta t}{2!\cdot 2}+o(\Delta t^2) \tag{12}$$
得到离散方程$$ 2\dfrac{u^{n+1}_i-u^{n+\frac{1}{2}}_i}{\Delta t}+\dfrac{\partial^2 u^n_i}{\partial t^2}\dfrac{\Delta t}{2!\cdot 2}+o(\Delta t^2)=\alpha \dfrac{u^{n+1}_{i+1}-2u^{n+1}_i+u^{n+1}_{i-1}}{\Delta x^2}+o(\Delta x^2). \tag{13}$$
将(10)和(13)两式相加，得到Crank-Nicolson格式$$ \dfrac{u^{n+1}_i-u^n_i}{\Delta t}=\alpha(\dfrac{u^n_{i+1}-2u^n_i+u^n_{i-1}}{\Delta x^2}+\dfrac{u^{n+1}_{i+1}-2u^{n+1}_i+u^{n+1}_{i-1}}{\Delta x^2}) \tag{14}$$易知CN格式在时间和空间上均为二阶精度，且为隐式格式。

# ↓ TODO ↓

写公式太累了，不写了。

## Inviscid Burgers Equation: Non-Conservative Form

## Inviscid Burgers Equation: Conservative Form

## 2D Poisson Equation

## Incompressible 2D NS Equation
