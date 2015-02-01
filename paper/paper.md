# 相对论重离子碰撞过程中背景磁场的计算与讨论

在相对论重离子碰撞过程中，带正电的原子核以极高的速度相互碰撞。对于非对心碰撞，会产生一个很强的磁场。这种背景磁场的存在是形成手征磁效应（CME: Chiral Magnetic Effect）的关键因素。手征磁效应指的是：一个带有非零手征性的系统，在磁场的作用下，产生一个沿磁场方向的电流的现象。在实验上观测到手征磁效应是具有非平凡拓扑结构的胶子构形的直接证据，并且也可作为QCD（Quantum Chromo Dynamics）中P与CP对称性破坏的一个信号。

本文主要讨论非对心重离子碰撞中背景磁场公式的推导，并对背景磁场的一些性质进行了研究。

## 背景磁场计算公式的推导

我们首先考虑匀速运动点电荷产生的磁场。匀速运动点电荷的场可以利用李纳-维谢尔势或者直接通过洛伦兹变换求得。下面分别用这两种方法来进行求解。

### 李纳-维谢尔势求点电荷磁场

设点电荷的位置矢量为$\vec{w}$，场点的位置矢量为$\vec{x}$，所带的电荷量为$q$，任意运动点电荷的李纳-维谢尔势为
$$V(\vec{x}, t) = \frac{1}{4 \pi \varepsilon_0} \frac{qc}{r c - \vec{r} \cdot \vec{v}}, \quad \vec{A}(\vec{x}, t) = \frac{\vec{v}}{c^2} V(\vec{x}, t)$$
其中$\vec{r}$是从推迟位置到场点$\vec{x}$的矢量，即$\vec{r} = \vec{x} - \vec{w}(t_r)$。$\vec{v}$是在推迟时间时速度的取值，即$\vec{v} = \dot{\vec{w}}(t_r)$。

我们利用下列势与场的关系来计算磁场：
$$\vec{E} = - \nabla V - \frac{\partial \vec{A}}{\partial t},   \vec{B} = \nabla \times \vec{A}$$
由于计算过程较复杂，故不在此推导，详细推导过程可参考文献[]。求得的结果为
$$\vec{E}(\vec{x}, t) = \frac{q}{4 \pi \varepsilon_0} \frac{r}{(\vec{r} \cdot \vec{u})^3} [(c^2 - v^2) \vec{u} + \vec{r} \times (\vec{u} \times \vec{a})]$$
$$\vec{B}(\vec{x}, t) = \frac{1}{c} \hat{\vec{r}} \times \vec{E}(\vec{x}, t)$$
其中$\vec{u} = c \hat{\vec{r}} - \vec{v}$，加速度$\vec{a}$也是在推迟时间取值的。

对于匀速运动点电荷，因此$\vec{a} = 0$，$\vec{w} = \vec{v} t$。将其带入之前求得的结果中，有
$$\vec{E}(\vec{x}, t) = \frac{q}{4 \pi \varepsilon_0} \frac{1 - v^2 / c^2}{(1 - v^2 \sin^2 \theta / c^2)^{3/2}} \frac{\hat{\vec{R}}}{R^2}$$
$$\vec{B} = \frac{1}{c} (\hat{\vec{r}} \times \vec{E}) = \frac{1}{c^2} (\vec{v} \times \vec{E})$$
其中$\vec{R} = \vec{x} - \vec{v} t$，是点电荷现在的位置到场点$\vec{x}$处的矢量。$\theta$是$\vec{R}$和$\vec{v}$之间的夹角。

### 洛伦兹变换求点电荷磁场



## 数值计算结果

