> [!NOTE] 
> For proper equation rendering, please view this documentation in day mode instead of night mode.

This repository provides the MATLAB implementation of the factorization method presented in the following paper: 

[Pu-Zhao Kow](https://puzhaokow1993.github.io/homepage/) and [Jenn-Nan Wang](http://www.math.ntu.edu.tw/~jnwang/), *Reconstruction of an impenetrable obstacle in anisotropic inhomogeneous background* IMA J. Appl. Math. **86** (2021), no. 2, 320--348, Inst. Math. Appl., [MR4246858](https://mathscinet.ams.org/mathscinet-getitem?mr=4246858), [Zbl:1471.76063](https://zbmath.org/1471.76063), [doi:10.1093/imamat/hxab002](https://doi.org/10.1093/imamat/hxab002). A corrigendum is available in [doi:10.1093/imamat/hxab046](https://doi.org/10.1093/imamat/hxab046) 

We will not going to explain all notations here, please refer to our paper for more details. 

Let ![A(x)=(a_{ij}(x))](https://latex.codecogs.com/png.image?\dpi{110}A(x)=(a_{ij}(x))) be a real-symmetric matrix with ![C^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}C^{\infty}) entries, which satisfies the uniform elliptic condition: there exists a constant ![0<c<1](https://latex.codecogs.com/png.image?\dpi{110}0<c<1) such that 
<div align="center">
  
![c|\xi|^{2}\le\sum_{ij}a_{ij}(x)\xi_{i}\xi_{j}\le%20c^{-1}|\xi|^{2}](https://latex.codecogs.com/png.image?\dpi{110}c|\xi|^{2}\le\sum_{ij}a_{ij}(x)\xi_{i}\xi_{j}\le%20c^{-1}|\xi|^{2})
</div>

for all ![x\in\mathbb{R}^{d}](https://latex.codecogs.com/png.image?\dpi{110}x\in\mathbb{R}^{d}),![\xi\in\mathbb{R}^{d}](https://latex.codecogs.com/png.image?\dpi{110}\xi\in\mathbb{R}^{d}), ![d=2](https://latex.codecogs.com/png.image?\dpi{110}d=2) or ![3](https://latex.codecogs.com/png.image?\dpi{110}3). Moreover, we assume that ![{\rm%20supp}(A-I)](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20supp}(A-I)) is compact in ![\mathbb{R}^{d}](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^{d}). Suppose that the acoustic refraction index ![n\in%20L^{\infty}(\mathbb{R}^{d})](https://latex.codecogs.com/png.image?\dpi{110}n\in%20L^{\infty}(\mathbb{R}^{d})) with ![n\ge%20c](https://latex.codecogs.com/png.image?\dpi{110}n\ge%20c), such that ![{\rm%20supp}(n-1)](https://latex.codecogs.com/png.image?\dpi{110}{\rm%20supp}(n-1)) is compact in ![\mathbb{R}^{d}](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^{d}). 

We model the (unknown) impenetrable obstacle by the open bounded domain ![D](https://latex.codecogs.com/png.image?\dpi{110}D) with Lipschitz boudary ![\partial%20D](https://latex.codecogs.com/png.image?\dpi{110}\partial%20D) such that ![\overline{D}\subset%20B_{R}](https://latex.codecogs.com/png.image?\dpi{110}\overline{D}\subset%20B_{R}) and ![\mathbb{R}^{d}\setminus\overline{D}](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^{d}\setminus\overline{D}) is connected. Let ![u_{D}^{\rm%20to}(x,\hat{z})=u_{A,n}^{\rm%20inc}(x,\hat{z})+u_{D}^{\rm%20sc}(x,\hat{z})](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\rm%20to}(x,\hat{z})=u_{A,n}^{\rm%20inc}(x,\hat{z})+u_{D}^{\rm%20sc}(x,\hat{z})) with an appropriate incident field ![u_{A,n}^{\rm%20inc}(x,\hat{z})](https://latex.codecogs.com/png.image?\dpi{110}u_{A,n}^{\rm%20inc}(x,\hat{z})) and the corresponding scattered field ![u_{D}^{\rm%20sc}(x,\hat{z})](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\rm%20sc}(x,\hat{z})) satisfy the following acoustic equation: 
<div align="center">
  
![\nabla\cdot(A(x)\nabla%20u_{D}^{\rm%20to})+k^{2}n(x)u_{D}^{\rm%20to}=0](https://latex.codecogs.com/png.image?\dpi{110}\nabla\cdot(A(x)\nabla%20u_{D}^{\rm%20to})+k^{2}n(x)u_{D}^{\rm%20to}=0) in ![\mathbb{R}^{d}\setminus\overline{D}](https://latex.codecogs.com/png.image?\dpi{110}\mathbb{R}^{d}\setminus\overline{D})
</div>
<div align="center">
  
![\mathcal{B}u_{D}^{\rm%20to}=0](https://latex.codecogs.com/png.image?\dpi{110}\mathcal{B}u_{D}^{\rm%20to}=0) on ![\partial%20D](https://latex.codecogs.com/png.image?\dpi{110}\partial%20D)
</div>
<div align="center">
  
![u_{D}^{\rm%20sc}](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\rm%20sc}) satisfies the Sommerfeld radiation condition at ![|x|\rightarrow\infty](https://latex.codecogs.com/png.image?\dpi{110}|x|\rightarrow\infty), 
</div>

where ![\mathcal{B}u_{D}^{\rm%20to}](https://latex.codecogs.com/png.image?\dpi{110}\mathcal{B}u_{D}^{\rm%20to}) is either Dirichlet or impedance boundary conditions. Let ![u_{D}^{\infty}(\hat{x},\hat{z})](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\infty}(\hat{x},\hat{z})) be the far-field pattern of ![u_{D}^{\rm%20sc}](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\rm%20sc}). The inverse problem is to determine ![D](https://latex.codecogs.com/png.image?\dpi{110}D) from the knowledge of ![u_{D}^{\infty}(\hat{x},\hat{z})](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\infty}(\hat{x},\hat{z})) for all ![\hat{x},\hat{z}\in\mathcal{S}^{d-1}](https://latex.codecogs.com/png.image?\dpi{110}\hat{x},\hat{z}\in\mathcal{S}^{d-1}). 

# Simulation of measurements # 

In our simulation, we take ![A\equiv%20I](https://latex.codecogs.com/png.image?\dpi{110}A\equiv%20I) and 
<div align="center">
  
![n(x,y)=1+e^{-\frac{1}{1+x^{2}+y^{2}}}](https://latex.codecogs.com/png.image?\dpi{110}n(x,y)=1+e^{-\frac{1}{1+x^{2}+y^{2}}}) for ![x^{2}+y^{2}<1](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}<1), otherwise ![n(x,y)=1](https://latex.codecogs.com/png.image?\dpi{110}n(x,y)=1)
</div>

which has jump discontinuities at ![x^{2}+y^{2}=1](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}=1). 

We first need to find suitable incident fields to stimulate the impenetrable obstacle (without considering the impenetrable obstacle ![D](https://latex.codecogs.com/png.image?\dpi{110}D)). We restrict the computational domain in the ball ![x^{2}+y^{2}<4](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}<4) and approximate the Sommerfeld radiation condition by the impedance condition: 
<div align="center">
  
![\Delta\tilde{u}_{\rm%20ref}^{\rm%20sc}+k^{2}n\tilde{u}_{\rm%20ref}^{\rm%20sc}=k^{2}(1-n)u_{\rm%20ref}^{\rm%20inc}](https://latex.codecogs.com/png.image?\dpi{110}\Delta\tilde{u}_{\rm%20ref}^{\rm%20sc}+k^{2}n\tilde{u}_{\rm%20ref}^{\rm%20sc}=k^{2}(1-n)u_{\rm%20ref}^{\rm%20inc}) for ![x^{2}+y^{2}<4](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}<4)
</div>
<div align="center">
  
![\partial_{r}\tilde{u}_{\rm%20ref}^{\rm%20sc}-\mathbf{i}k\tilde{u}_{\rm%20ref}^{\rm%20sc}=0](https://latex.codecogs.com/png.image?\dpi{110}\partial_{r}\tilde{u}_{\rm%20ref}^{\rm%20sc}-\mathbf{i}k\tilde{u}_{\rm%20ref}^{\rm%20sc}=0) for ![x^{2}+y^{2}=4](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}=4)
</div>

We solve this boundary value problem FEM with mesh size ![\le0.1](https://latex.codecogs.com/png.image?\dpi{110}\le0.1). Then we can approximate the total field by ![\tilde{u}_{\rm%20ref}^{\rm%20to}=\tilde{u}_{\rm%20ref}^{\rm%20sc}+u_{\rm%20ref}^{\rm%20inc}](https://latex.codecogs.com/png.image?\dpi{110}\tilde{u}_{\rm%20ref}^{\rm%20to}=\tilde{u}_{\rm%20ref}^{\rm%20sc}+u_{\rm%20ref}^{\rm%20inc}). 

We now stimulate the impenetrable obstacle ![D](https://latex.codecogs.com/png.image?\dpi{110}D) by using the incident field 
<div align="center">
  
![\tilde{u}_{I,n}^{\rm%20inc}(x,\hat{z})=\overline{\tilde{u}_{\rm%20ref}^{\rm%20to}(x,-\hat{z})}](https://latex.codecogs.com/png.image?\dpi{110}\tilde{u}_{I,n}^{\rm%20inc}(x,\hat{z})=\overline{\tilde{u}_{\rm%20ref}^{\rm%20to}(x,-\hat{z})}) 
</div>

and the approximated scattered field is obtained by solving 
<div align="center">
  
![\Delta\tilde{u}_{D}^{\rm%20sc}+k^{2}n\tilde{u}_{D}^{\rm%20sc}=k^{2}(1-n)\tilde{u}_{I,n}^{\rm%20inc}](https://latex.codecogs.com/png.image?\dpi{110}\Delta\tilde{u}_{D}^{\rm%20sc}+k^{2}n\tilde{u}_{D}^{\rm%20sc}=k^{2}(1-n)\tilde{u}_{I,n}^{\rm%20inc}) for ![x^{2}+y^{2}<4](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}<4) with ![(x,y)\notin%20D](https://latex.codecogs.com/png.image?\dpi{110}(x,y)\notin%20D)
</div>
<div align="center">
  
![\mathcal{B}(\tilde{u}_{D}^{\rm%20sc}+\tilde{u}_{I,n}^{\rm%20sc})=0](https://latex.codecogs.com/png.image?\dpi{110}\mathcal{B}(u_{D}^{\rm%20sc}+\tilde{u}_{I,n}^{\rm%20sc})=0) on ![\partial%20D](https://latex.codecogs.com/png.image?\dpi{110}\partial%20D)
</div>
<div align="center">
  
![\partial_{r}\tilde{u}_{D}^{\rm%20sc}-\mathbf{i}k\tilde{u}_{D}^{\rm%20sc}=0](https://latex.codecogs.com/png.image?\dpi{110}\partial_{r}\tilde{u}_{D}^{\rm%20sc}-\mathbf{i}k\tilde{u}_{D}^{\rm%20sc}=0) for ![x^{2}+y^{2}=4](https://latex.codecogs.com/png.image?\dpi{110}x^{2}+y^{2}=4)
</div>

and the far-field pattern ![u_{D}^{\infty}](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\infty}) of ![u_{D}^{\rm%20sc}](https://latex.codecogs.com/png.image?\dpi{110}u_{D}^{\rm%20sc}) is approximated by 
<div align="center">
  
![\tilde{u}_{D}^{\rm%20sc}(x,y)=\frac{e^{\frac{\mathbf{i}\pi}{4}}}{\sqrt{8\pi%20k}}\frac{e^{\mathbf{i}k|(x,y)|}}{|(x,y)|^{\frac{1}{2}}}\tilde{u}_{D}^{\infty}(\frac{(x,y)}{|(x,y)|})](https://latex.codecogs.com/png.image?\dpi{110}\tilde{u}_{D}^{\rm%20sc}(x,y)=\frac{e^{\frac{\mathbf{i}\pi}{4}}}{\sqrt{8\pi%20k}}\frac{e^{\mathbf{i}k|(x,y)|}}{|(x,y)|^{\frac{1}{2}}}\tilde{u}_{D}^{\infty}\left(\frac{(x,y)}{|(x,y)|}\right)) on ![|(x,y)|=3](https://latex.codecogs.com/png.image?\dpi{110}|(x,y)|=3). 
</div>

Here, we choose ![|(x,y)|=3](https://latex.codecogs.com/png.image?\dpi{110}|(x,y)|=3) rather than ![|(x,y)|=4](https://latex.codecogs.com/png.image?\dpi{110}|(x,y)|=4) to reduce the effect of the reflected wave from the boudnary. Finally, we approximate the far-field operator ![F_{D}](https://latex.codecogs.com/png.image?\dpi{110}F_{D}) as an ![M\times%20M](https://latex.codecogs.com/png.image?\dpi{110}M\times%20M) matrix: 
</div>
<div align="center">
  
![\tilde{F}_{D}=(F_{j\ell})_{j,\ell=1}^{M}](https://latex.codecogs.com/png.image?\dpi{110}\tilde{F}_{D}=(F_{j\ell})_{j,\ell=1}^{M}) 
</div>

with the far-field pattern ![\tilde{u}_{D}^{\infty}(\theta_{j},\theta_{\ell})](https://latex.codecogs.com/png.image?\dpi{110}\tilde{u}_{D}^{\infty}(\theta_{j},\theta_{\ell})) for ![j,\ell=1,\cdots,M](https://latex.codecogs.com/png.image?\dpi{110}j,\ell=1,\cdots,M) at equidistantly distributed directions ![\theta_{j}=2\pi%20j/M](https://latex.codecogs.com/png.image?\dpi{110}\theta_{j}=2\pi%20j/M). In our numerical simulation, we choose ![M=24](https://latex.codecogs.com/png.image?\dpi{110}M=24). 

# Numerical reconstriction # 




[comment]: <> (https://docs.github.com/en/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)
