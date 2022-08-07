```@meta
CurrentModule = VectorSpinVlasovMaxwell1D1V
```

# VectorSpinVlasovMaxwell1D1V

## Scalar spin laser plasma model

Particle distribution function $f(x, p, {\mathbf s}, t)$, $x\in [0,L], p\in \mathbb{R} $ are scalars, ${\mathbf s}=(s_1,s_2,s_3) \in \mathbb{R}^3$, ${\mathbf E} = (E_x, {\mathbf E}_\perp) = (E_x, E_y, E_z)$, ${\mathbf A} = (A_x, {\mathbf A}_\perp) = (0, A_y, A_z)$ and ${\mathbf B} =\nabla\times{\mathbf  A} = (0,- \partial_xA_z,  \partial_xA_y)$.

The scalar spin Vlasov--Maxwell  system introduced in \cite{crouseilles2021geometric} is given by 

```math
\begin{equation}\label{eq:reduced}
\left\{
\begin{aligned}
&\frac{\partial f}{\partial t} + p \frac{\partial f}{\partial x} + [ E_x - \mathfrak{h} s_2 \frac{\partial^2 A_z}{\partial x^2} + \mathfrak{h} s_3 \frac{\partial^2 A_y}{\partial x^2}  - {\mathbf A}_\perp \cdot \frac{\partial {\mathbf A}_\perp}{\partial x} ]\frac{\partial f}{\partial p}  \\ 
& \hspace{3cm}+ [s_3 \frac{\partial A_z}{\partial x} + s_2 \frac{\partial A_y}{\partial x}, -s_1 \frac{\partial A_y}{\partial x}, -s_1 \frac{\partial A_z}{\partial x} ] \cdot \frac{\partial f}{\partial {\mathbf s}} = 0,\\
&\frac{\partial E_x}{\partial t} = -\int_{\mathbb{R}^4} p f  \mathrm{d}{p}\mathrm{d}\mathrm{\mathbf s},\\
&\frac{\partial E_y}{\partial t} = - \frac{\partial^2 A_y}{\partial x^2} + A_y \int_{\mathbb{R}^4}  f  \mathrm{d}{p}\mathrm{d}\mathrm{\mathbf s} + \mathfrak{h}\int_{\mathbb{R}^4} s_3 \frac{\partial f}{\partial x}\mathrm{d}{p}\mathrm{d}\mathrm{\mathbf s},\\
&\frac{\partial E_z}{\partial t} = - \frac{\partial^2 A_z}{\partial x^2} + A_z \int_{\mathbb{R}^4}  f  \mathrm{d}{p}\mathrm{d}\mathrm{\mathbf s} - \mathfrak{h}\int_{\mathbb{R}^4} s_2 \frac{\partial f}{\partial x}\mathrm{d}{p}\mathrm{d}\mathrm{\mathbf s},\\
& \frac{\partial {\mathbf A}_\perp}{\partial t} = - {\mathbf E}_\perp,\\
&\frac{\partial E_x}{\partial x} = \int_{\mathbb{R}^4} f \mathrm{d}{p}\mathrm{d}\mathrm{\mathbf s} - 1. \ \text{(Poisson equation)}
\end{aligned}
\right.
\end{equation}
```
