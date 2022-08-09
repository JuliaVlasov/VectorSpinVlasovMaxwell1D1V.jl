var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"CurrentModule = VectorSpinVlasovMaxwell1D1V","category":"page"},{"location":"api/#VectorSpinVlasovMaxwell1D1V","page":"API","title":"VectorSpinVlasovMaxwell1D1V","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Documentation for VectorSpinVlasovMaxwell1D1V.","category":"page"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [VectorSpinVlasovMaxwell1D1V]","category":"page"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.Mesh","page":"API","title":"VectorSpinVlasovMaxwell1D1V.Mesh","text":"struct Mesh\n\nMesh type to store domain parameters\n\nN::Int64\nNumber of points in v\nM::Int64\nNumber of points in x\nH::Float64\nDomain size v ∈ ]-H,+H[\nL::Float64\nDomain size x ∈ [0,L]\nk::Vector{Float64}\nWave number vector to compute derivative with FFTs\ndx::Float64\nSize step along x\ndv::Float64\nSize step along v\nx::Vector{Float64}\npoints along x direction\nv::Vector{Float64}\npoints along v direction\n\n\n\n\n\n","category":"type"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H1f!-NTuple{8, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H1f!","text":"H1f!(f0, f1, f2, f3, E1, t, L, H)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H2fh!-NTuple{10, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H2fh!","text":"compute the subsystem H2\n\nH2fh!(f0, f1, f2, f3, E3, A3, t, L, H, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.H3fh!-NTuple{10, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.H3fh!","text":"H3fh!(f0, f1, f2, f3, E2, A2, t, L, H, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.HAA!-NTuple{11, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.HAA!","text":"HAA!(f0, f1, f2, f3, E2, E3, A2, A3, t, L, H)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.He!-NTuple{11, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.He!","text":"He!(f0, f1, f2, f3, E1, E2, E3, A2, A3, t, H)\n\n\ncompute the first subsystem He()\n\nM is odd number\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.numeint-Tuple{Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.numeint","text":"computation of integral average in each cell using newton-cotes formula\n\nnumeint(value, N)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-NTuple{9, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E2, A2, op, t, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, Any, Any, Any, Any, HeOperator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, op, t)\n\n\ncompute the first subsystem He()\n\nM is odd number\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, Any, Any, Any, HAAOperator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E2, E3, A2, A3, op, t)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, Any, H2fhOperator, Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"compute the subsystem H2\n\nstep!(f0, f1, f2, f3, E3, A3, op, t, h_int)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.step!-Tuple{Any, Any, Any, Any, Any, H1fOperator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.step!","text":"step!(f0, f1, f2, f3, E1, op, t)\n\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.translation!-Tuple{Any, Any, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.translation!","text":"translation!(df, delta, H)\n\n\ninterpolate df(x - delta)\n\n\n\n\n\n","category":"method"},{"location":"api/#VectorSpinVlasovMaxwell1D1V.translation!-Tuple{Any, Translator, Any}","page":"API","title":"VectorSpinVlasovMaxwell1D1V.translation!","text":"translation!(df, adv, delta)\n\n\ninterpolate df(x - delta) with Parabolic Spline Method (PSM) \n\nWe consider a linear advection problem in p direction\n\nfracpartial fpartial t + a fracpartial fpartial x =0\n\nFrom the conservation of the volume, we have the following identity\n\nf_jell(t)=frac1Delta p int_p_ell-12 ^p_ell+12 f(x_jpt)mathrmdp =\frac1Delta p int_p_ell-12-at ^p_ell+12-at f(x_jp0)mathrmdp\n\nFor simplicity, denote by qin 1 M the index such that p_ell+12-at in p_q-12p_q+12 i.e. p_ell+12-at in C_q, then we have\n\nf_jell(t) =frac1Delta p int_p_q-12-at ^p_q-12 f(x_jp0)mathrmdp+f_jq(0)-frac1Delta p int_p_q+12 ^p_q+12-at f(x_jp0)mathrmdp\n\nHere we need to reconstruct a polynomial function f(x_jp0) using the averages f_jl(0) using the PSM approach. \n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = VectorSpinVlasovMaxwell1D1V","category":"page"},{"location":"#VectorSpinVlasovMaxwell1D1V","page":"Home","title":"VectorSpinVlasovMaxwell1D1V","text":"","category":"section"},{"location":"#Scalar-spin-laser-plasma-model","page":"Home","title":"Scalar spin laser plasma model","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Particle distribution function f(x p mathbf s t), ","category":"page"},{"location":"","page":"Home","title":"Home","text":"xin 0L, \np in mathbbR are scalars, \nmathbf s=(s_1s_2s_3) in mathbbR^3, \nmathbf E = (E_x mathbf E_perp) = (E_x E_y E_z), \nmathbf A = (A_x mathbf A_perp) = (0 A_y A_z), \nmathbf B =nablatimesmathbf  A = (0- partial_xA_z  partial_xA_y).","category":"page"},{"location":"","page":"Home","title":"Home","text":"The scalar spin Vlasov–Maxwell  system is:","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial fpartial t + p fracpartial fpartial x +  E_x - mathfrakh s_2 fracpartial^2 A_zpartial x^2 + mathfrakh s_3 fracpartial^2 A_ypartial x^2  - mathbf A_perp cdot fracpartial mathbf A_perppartial x fracpartial fpartial p   \n hspace3cm+ s_3 fracpartial A_zpartial x + s_2 fracpartial A_ypartial x -s_1 fracpartial A_ypartial x -s_1 fracpartial A_zpartial x  cdot fracpartial fpartial mathbf s = 0\nfracpartial E_xpartial t = -int_mathbbR^4 p f  mathrmdpmathrmdmathrmmathbf s\nfracpartial E_ypartial t = - fracpartial^2 A_ypartial x^2 + A_y int_mathbbR^4  f  mathrmdpmathrmdmathrmmathbf s + mathfrakhint_mathbbR^4 s_3 fracpartial fpartial xmathrmdpmathrmdmathrmmathbf s\nfracpartial E_zpartial t = - fracpartial^2 A_zpartial x^2 + A_z int_mathbbR^4  f  mathrmdpmathrmdmathrmmathbf s - mathfrakhint_mathbbR^4 s_2 fracpartial fpartial xmathrmdpmathrmdmathrmmathbf s\n fracpartial mathbf A_perppartial t = - mathbf E_perp\nfracpartial E_xpartial x = int_mathbbR^4 f mathrmdpmathrmdmathrmmathbf s - 1  text(Poisson equation)\nendaligned\nright","category":"page"},{"location":"","page":"Home","title":"Home","text":"The system numerically solve is the vector model:","category":"page"},{"location":"","page":"Home","title":"Home","text":"f(t xpmathbfs)=frac14pi(f_0(t xp)+3s_1f_1(t xp)+3s_2f_2(t xp)+3s_3f_3(t xp))","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial f_0partial t + p fracpartial f_0partial x + left(E_x - mathbf A_perp cdot fracpartial mathbf A_perppartial x  right) fracpartial f_0partial p - mathfrakhfracpartial^2 A_zpartial x^2fracpartial f_2partial p +  mathfrakhfracpartial^2 A_ypartial x^2 fracpartial f_3partial p  = 0\nfracpartial f_1partial t + p fracpartial f_1partial x + left(E_x - mathbf A_perp cdot fracpartial mathbf A_perppartial x  right) fracpartial f_1partial p\n - fracpartial A_z partial x  f_3  -  fracpartial A_y partial x f_2 = 0\n fracpartial f_2partial t + p fracpartial f_2partial x + left(E_x - mathbf A_perp cdot fracpartial mathbf A_perppartial x  right) fracpartial f_2partial p - fracmathfrakh3 fracpartial^2 A_zpartial x^2fracpartial f_0partial p\n  +  fracpartial A_y partial x f_1 = 0\n  fracpartial f_3partial t + p fracpartial f_3partial x + left(E_x - mathbf A_perp cdot fracpartial mathbf A_perppartial x  right) fracpartial f_3partial p + fracmathfrakh3  fracpartial^2 A_ypartial x^2fracpartial f_0partial p\n  +  fracpartial A_z partial x f_1 = 0\nfracpartial E_xpartial t = -int_mathbbR p f_0  mathrmdmathrmp\nfracpartial E_ypartial t = - fracpartial^2 A_ypartial x^2 + A_y int_mathbbR  f_0  mathrmdmathrmp + mathfrakhint_mathbbR  fracpartial f_3partial xmathrmdmathrmp\nfracpartial E_zpartial t = - fracpartial^2 A_zpartial x^2 + A_z int_mathbbR  f_0  mathrmdmathrmp -mathfrakh int_mathbbR  fracpartial f_2partial xmathrmdmathrmp\n fracpartial mathbf A_perppartial t = - mathbf E_perp\nfracpartial E_xpartial x = int_mathbbR f_0 mathrmdmathrmp - 1 text(Poisson equation)\nendaligned\nright","category":"page"},{"location":"#Time-discretization:-Hamiltonian-splitting-method","page":"Home","title":"Time discretization: Hamiltonian splitting method","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"mathcalH=mathcalH_p+mathcalH_A+mathcalH_E+mathcalH_2+mathcalH_3","category":"page"},{"location":"","page":"Home","title":"Home","text":"where ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nmathcalH_p  = frac12int p^2 f_0 mathrmd xmathrmdp\nmathcalH_A = frac12int mathbf A_perp^2 f_0 mathrmd xmathrmd p+frac12int leftfracpartial mathbf A_perppartial xright^2 mathrmdx\nmathcalH_E = frac12int mathbf E^2  mathrmd x=frac12int (E_x^2+mathbf E_perp^2 ) mathrmd x \nmathcalH_2 = int_Omega  mathfrakh f_2 fracpartial A_zpartial x mathrmdxmathrmdp\nmathcalH_3 = int_Omega  -mathfrakh f_3 fracpartial A_ypartial x mathrmdxmathrmdp\nendaligned ","category":"page"},{"location":"#Subsystem-for-\\mathcal{H}_p","page":"Home","title":"Subsystem for mathcalH_p","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The subsystem fracpartial mathcalZpartial t =  mathcalZ mathcalH_p  associated to mathcalH_p = frac12int p^2 f_0 mathrmd xmathrmdp is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial f_0partial t = f_0 mathcalH_p  = -pfracpartial f_0partial x \nfracpartial mathbffpartial t = mathbff mathcalH_p = -pfracpartial mathbffpartial x \nfracpartial f_1partial t = f_1 mathcalH_p = -pfracpartial f_1partial x \nfracpartial f_2partial t = f_2 mathcalH_p = -pfracpartial f_2partial x \nfracpartial f_3partial t = f_3 mathcalH_p = -pfracpartial f_3partial x\n fracpartial E_xpartial t =  E_x mathcalH_p  =- int_mathbbR p f_0mathrmdp\n fracpartial mathbf E_perppartial t =fracpartial mathbf A_perppartial t =0 \n fracpartial mathbf E_perppartial t =   mathbf E_perp mathcalH_p   = mathbf 0 fracpartial mathbf A_perppartial t =   mathbf A_perp mathcalH_p   = mathbf 0\nendaligned\nright","category":"page"},{"location":"","page":"Home","title":"Home","text":"We denote the initial value as (f_0^0(xp)mathbff^0(xp) E_x^0(x)mathbf E_perp^0(x)mathbf A_perp^0(x)) at time t=0. The solution at time t of this subsystem can be written explicitly, ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n        f_0(xpt)=f_0^0(x-ptp)   mathbff(xpt)=mathbff^0(x-ptp) f_2(xpt)=f_2^0(x-ptp) f_3(xpt)=f_3^0(x-ptp) \n        E_x(xt)=E_x^0(x)-int_0^tint_mathbbR pf_0(xptau) mathrmdpmathrmdtau=E_x^0(x)-int_0^tint_mathbbR pf_0^0(x-ptaup) mathrmdpmathrmdtau \n        mathbf E_perp(xt)=mathbf E_perp^0(x)  mathbf A_perp(xt)=mathbf A_perp^0(x) \n    endaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"Next, we check that the solution propagates the Poisson equation. To do so, we assume that the Poisson equation holds initially, i.e.","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracpartial E_x^0partial x=int_mathbbR f_0^0mathrmdp-1","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then we have, by differentiating the expression of E_x(t x) with respect to x ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nfracpartial E_x(xt)partial x=fracpartial E_x^0partial x-int_0^tint_mathbbR pfracpartial f_0^0(x-ptaup)partial x mathrmdpmathrmdtau=fracpartial E_x^0partial x+int_0^tint_mathbbR fracpartial f_0^0(x-ptaup)partial tau mathrmdpmathrmdtau\n=fracpartial E_x^0partial x+int_mathbbR  f_0^0(x-ptp)mathrmdp-int_mathbbR  f_0^0(xp)mathrmdp=int_mathbbR  f_0(xpt)mathrmdp-1 \nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"which proves that the Poisson equation is satisfied at time t.","category":"page"},{"location":"#Subsystem-for-\\mathcal{H}_A","page":"Home","title":"Subsystem for mathcalH_A","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The subsystem fracpartial mathcalZpartial t =  mathcalZ mathcalH_A  associated to the sub-Hamiltonian mathcalH_A = frac12int mathbf A_perp^2 f_0 mathrmdmathbf xmathrmdmathbf p+frac12int fracpartial mathbf A_perppartial x^2 mathrmdmathbf x is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial f_0partial t = f_0 mathcalH_A  = mathbf A_perp cdot fracpartial mathbf A_perppartial xfracpartial f_0partial p \nfracpartial mathbf fpartial t = mathbf f mathcalH_A = mathbf A_perp cdot fracpartial mathbf A_perppartial x fracpartial mathbf fpartial p \nfracpartial f_2partial t = f_2 mathcalH_A = mathbf A_perp cdot fracpartial mathbf A_perppartial x fracpartial f_2partial p \nfracpartial f_3partial t = f_3 mathcalH_A = mathbf A_perp cdot fracpartial mathbf A_perppartial x fracpartial f_3partial p\n fracpartial E_xpartial t =  E_x mathcalH_A  = 0\n fracpartial mathbf E_perppartial t =   mathbf E_perp mathcalH_A   = - fracpartial^2 mathbf A_perppartial x^2 + mathbf A_perp int_mathbbR  f_0  mathrmdmathrmp\n fracpartial E_xpartial t =fracpartial mathbf A_perppartial t = 0   mathbf A_perp mathcalH_A   = mathbf 0\nendaligned\nright","category":"page"},{"location":"","page":"Home","title":"Home","text":"We denote by (f_0^0(xp)mathbf f^0(xp) E_x^0(x)mathbf E_perp^0(x)mathbf A_perp^0(x)) the initial value at time t=0. The exact solution at time t is,","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n        f_0(xpt)=f_0^0 left( xp+tmathbf A_perp^0(x) cdot fracpartial mathbf A_perp^0(x)partial x  right) \n   mathbf f(xpt)=mathbf f^0 left( xp+tmathbf A_perp^0(x) cdot fracpartial mathbf A_perp^0(x)partial x  right) \n        mathbf E_perp(xt)=mathbf E_perp^0(x)-tfracpartial^2 mathbf A_perp^0(x)partial x^2 +tmathbf A_perp^0(x)int_mathbbR f_0^0(xp) mathrmdp \n        E_x(xt)=E_x^0(x)  mathbf A_perp(xt)=mathbf A_perp^0(x) \nendaligned","category":"page"},{"location":"#Subsystem-for-\\mathcal{H}_E","page":"Home","title":"Subsystem for mathcalH_E","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The subsystem fracpartial mathcalZpartial t =  mathcalZ mathcalH_E  associated to the sub-Hamiltonian mathcalH_E = frac12int mathbf E^2  mathrmdmathbf x is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial f_0partial t = f_0 mathcalH_E  = -E_x fracpartial f_0partial p \nfracpartial mathbf fpartial t = mathbf f mathcalH_E = -E_x fracpartial mathbf fpartial p \nfracpartial f_2partial t = f_2 mathcalH_E = -E_x fracpartial f_2partial p \nfracpartial f_3partial t = f_3 mathcalH_E =-E_x fracpartial f_3partial p\n fracpartial E_xpartial t =  E_x mathcalH_E  = 0\n fracpartial mathbf E_perppartial t =   mathbf E_perp mathcalH_E   = mathbf 0\n fracpartial E_xpartial t = fracpartial mathbf A_perppartial t = 0   mathbf A_perp mathcalH_E   =  -mathbf E_perp\nendaligned\nright","category":"page"},{"location":"","page":"Home","title":"Home","text":"With the initial value(f_0^0(xp)mathbf f^0(xp)E_x^0(x)mathbf E_perp^0(x)mathbf A_perp^0(x)) at time t=0, the solution at time t is as follows, ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n f_0(xpt)=f_0^0 ( xp-tE_x^0(x) ) \n mathbf f(xpt)=mathbf f^0 ( xp-tE_x^0(x) )  \n E_x(xt)=E_x^0(x) \n mathbf E_perp(xt)=mathbf E_perp^0(x) \n mathbf A_perp(xt)=mathbf A_perp^0(x) -tmathbf E_perp^0(x)\nendaligned","category":"page"},{"location":"#Subsystem-for-\\mathcal{H}_2","page":"Home","title":"Subsystem for mathcalH_2","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The subsystem fracpartial mathcalZpartial t =  mathcalZ mathcalH_2  associated to the sub-Hamiltonian mathcalH_2 = int_Omega  mathfrakh f_2 fracpartial A_zpartial x mathrmdxmathrmdp is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial f_0partial t = f_0 mathcalH_2  = mathfrakhfracpartial^2 A_zpartial x^2fracpartial f_2partial p  \nfracpartial f_1partial t = f_1 mathcalH_2 = fracpartial A_z partial x  f_3 \nfracpartial f_2partial t = f_2 mathcalH_2 = fracmathfrakh3 fracpartial^2 A_zpartial x^2fracpartial f_0partial p \nfracpartial f_3partial t = f_3 mathcalH_2 = -fracpartial A_z partial x f_1\n fracpartial  E_zpartial t =   E_z mathcalH_2   =-mathfrakh int_mathbbR  fracpartial f_2partial xmathrmdp\n fracpartial E_xpartial t =fracpartial  E_ypartial t =fracpartial mathbf A_perppartial t =0 \nendaligned\nright","category":"page"},{"location":"","page":"Home","title":"Home","text":"In this subsystem, we observe some coupling between the distribution functions. To write down the exact solution, we reformulate the equations on (f_0 mathbf f) as, using A_z(x t)=A_z^0(x)","category":"page"},{"location":"","page":"Home","title":"Home","text":"    beginaligned\n     partial_t    beginpmatrix\n            f_1  \n            f_3\n        endpmatrix-fracpartial A_z^0 partial x J beginpmatrix\nf_1  \nf_3\nendpmatrix =0  \n partial_t    beginpmatrix\n    f_0 \n    f_2\nendpmatrix-mathfrakhfracpartial^2 A_z^0partial x^2 beginpmatrix\n    0  1 \n    frac13 0\nendpmatrix partial_p    beginpmatrix\nf_0 \nf_2\nendpmatrix =0 \n    endaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"where J denotes the symplectic matrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"J=beginpmatrix\n    0  1 \n-1  0\nendpmatrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"With the initial value (f_0^0(xp)mathbf f^0(xp)E_x^0(x)mathbf E_perp^0(x)mathbf A_perp^0(x)) at time t=0, the exact solution for the first system is","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginpmatrix\n            f_1  \n            f_3\n        endpmatrix(xpt)=expleft(fracpartial A_z^0(x) partial x J tright)beginpmatrix\nf_1^0(xp)  \nf_3^0(xp)\nendpmatrix  textwith exp(Js)=beginpmatrix\n    cos(s)  sin(s) \n-sin(s)  cos(s)\nendpmatrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"Let us now focus on the second system","category":"page"},{"location":"","page":"Home","title":"Home","text":"partial_t  beginpmatrix\n    f_0 \n    f_2\nendpmatrix-mathfrakhfracpartial^2 A_z^0partial x^2 beginpmatrix\n    0  1 \n    frac13 0\nendpmatrix partial_p    beginpmatrix\nf_0 \nf_2\nendpmatrix =0 ","category":"page"},{"location":"","page":"Home","title":"Home","text":"By the eigen-decomposition ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginpmatrix\n    frac12  fracsqrt32 \n    frac12 -fracsqrt32\nendpmatrix\nbeginpmatrix\n    0  1 \n    frac13 0\nendpmatrix\nbeginpmatrix\n    1  1 \n    frac1sqrt3 frac-1sqrt3\nendpmatrix\n=beginpmatrix\nfrac1sqrt3 0 \n    0 -frac1sqrt3\nendpmatrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"then, one can diagonalize the transport equation to get","category":"page"},{"location":"","page":"Home","title":"Home","text":"partial_t  beginpmatrix\n    frac12f_0+fracsqrt32f_2 \n    frac12f_0-fracsqrt32f_2 \nendpmatrix- \nfracmathfrakhsqrt3fracpartial^2 A_z^0partial x^2\nbeginpmatrix\n    1  0 \n    0 -1\nendpmatrix partial_p    beginpmatrix\n    frac12f_0+fracsqrt32f_2 \n    frac12f_0-fracsqrt32f_2 \nendpmatrix =0\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"Thus, we can solve the transport equation","category":"page"},{"location":"","page":"Home","title":"Home","text":"Big(frac12f_0pmfracsqrt32f_2Big)(xpt)=Big(frac12f_0^0 pmfracsqrt32f_2^0Big)( xppm tfracmathfrakhsqrt3fracpartial^2 A_z^0partial x^2(x))","category":"page"},{"location":"","page":"Home","title":"Home","text":"and compute the exact solution at time t as follows, ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n        f_1(xpt)=cos(t fracpartial A_z^0(x) partial x )f_1^0 ( xp)+sin(t fracpartial A_z^0(x) partial x )f_3^0 ( xp) \n        f_3(xpt)=-sin(t fracpartial A_z^0(x) partial x )f_1^0 ( xp)+cos(t fracpartial A_z^0(x) partial x )f_3^0 ( xp) \n        f_0(xpt)=Big(frac12f_0^0 +fracsqrt32f_2^0Big)( xp+tfracmathfrakhsqrt3fracpartial^2 A_z^0partial x^2(x) )+Big(frac12f_0^0 -fracsqrt32f_2^0Big)( xp-tfracmathfrakhsqrt3fracpartial^2 A_z^0partial x^2(x) )\n        f_2(xpt)=frac1sqrt3Big(frac12f_0^0 +fracsqrt32f_2^0Big)( xp+tfracmathfrakhsqrt3fracpartial^2 A_z^0partial x^2(x) )-frac1sqrt3Big(frac12f_0^0 -fracsqrt32f_2^0Big)( xp-tfracmathfrakhsqrt3fracpartial^2 A_z^0partial x^2(x) )\n        mathbf A_perp(xt)=mathbf A_perp^0(x) E_x(xt)=E_x^0(x) E_y(xt)=E_y^0(x)\n        E_z(xt)=E_z^0(x)-tmathfrakhint_mathbbR fracpartial f_2^0partial xmathrmdp\n    endaligned","category":"page"},{"location":"#Subsystem-for-\\mathcal{H}_3","page":"Home","title":"Subsystem for mathcalH_3","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The subsystem fracpartial mathcalZpartial t =  mathcalZ mathcalH_3  associated to the sub-Hamiltonian mathcalH_3 = -int_Omega mathfrakh f_3 fracpartial A_ypartial x mathrmdxmathrmdp is ","category":"page"},{"location":"","page":"Home","title":"Home","text":"left\nbeginaligned\nfracpartial f_0partial t = f_0 mathcalH_3  = -mathfrakhfracpartial^2 A_ypartial x^2 fracpartial f_3partial p  \nfracpartial f_1partial t = f_1 mathcalH_3 =  fracpartial A_y partial x f_2 \nfracpartial f_2partial t = f_2 mathcalH_3 = - fracpartial A_y partial x f_1\nfracpartial f_3partial t =- f_3 mathcalH_3 = -fracmathfrakh3 fracpartial^2 A_ypartial x^2 fracpartial f_0partial p \n fracpartial  E_ypartial t =   E_y mathcalH_3   = mathfrakh int_mathbbR  fracpartial f_3partial xmathrmdp\n fracpartial E_xpartial t = fracpartial  E_zpartial t = fracpartial mathbf A_perppartial t =0 \nendaligned\nright","category":"page"},{"location":"","page":"Home","title":"Home","text":"This subsystem is very similar to the mathcalH_2 one, hence, as previously, we reformulate the equations on the distribution functions as ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n         partial_t    beginpmatrix\n            f_1  \n            f_2\n        endpmatrix\n        -fracpartial A_y^0 partial x J \n        beginpmatrix\n           0  -1 \n           1  0\n        endpmatrix \n        beginpmatrix\n        f_1  \n        f_2\n    endpmatrix =0 \n         partial_t    beginpmatrix\n            f_0 \n            f_3\n        endpmatrix+mathfrakhfracpartial^2 A_y^0partial x^2 beginpmatrix\n            0  1 \n            frac13 0\n        endpmatrix partial_p    beginpmatrix\n            f_0 \n            f_3\n        endpmatrix =0 \n     fracpartial E_xpartial t= 0 \n fracpartial  E_ypartial t =   mathfrakh int_mathbbR  fracpartial f_3partial xmathrmdp\n     fracpartial  E_zpartial t =   E_z H_3   =  0 \n     fracpartial mathbf A_perppartial t =   mathbf A_perp H_3   =  0\n    endaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"with initial value","category":"page"},{"location":"","page":"Home","title":"Home","text":"(f_0^0(xp)mathbf f^0(xp)E_x^0(x)mathbf E_perp^0(x)mathbf A_perp^0(x)) at time t=0. We derive similar formula with mathcalH_2 for the exact solution at time t ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n        f_1(xpt)=cos(t fracpartial A_y^0(x) partial x )f_1^0 ( xp)+sin(t fracpartial A_y^0(x) partial x )f_2^0 ( xp) \n        f_2(xpt)=-sin(t fracpartial A_y^0(x) partial x )f_1^0 ( xp)+cos(t fracpartial A_y^0(x) partial x )f_2^0 ( xp) \n        f_0(xpt)=Big(frac12f_0^0 +fracsqrt32f_3^0Big)( xp-tfracmathfrakhsqrt3fracpartial^2 A_y^0partial x^2(x) )+Big(frac12f_0^0 -fracsqrt32f_3^0Big)( xp+tfracmathfrakhsqrt3fracpartial^2 A_y^0partial x^2(x) )\n        f_3(xpt)=frac1sqrt3Big(frac12f_0^0 +fracsqrt32f_3^0Big)( xp-tfracmathfrakhsqrt3fracpartial^2 A_y^0partial x^2(x) )-frac1sqrt3Big(frac12f_0^0 -fracsqrt32f_3^0Big)( xp+tfracmathfrakhsqrt3fracpartial^2 A_y^0partial x^2(x) )\n    endaligned ","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\n E_y(xt)=E_y^0(x)+tmathfrakhint_mathbbR fracpartial f_3^0partial xmathrmdp \n mathbf A_perp(xt)=mathbf A_perp^0(x) \n E_x(xt)=E_x^0(x) \n E_z(xt)=E_z^0(x) \nendaligned ","category":"page"},{"location":"","page":"Home","title":"Home","text":"To compute the solution E_y(xt), we use the fact that","category":"page"},{"location":"","page":"Home","title":"Home","text":"int_mathbbR f_3(xpt) mathrmdp =int_mathbbR f_3^0(xp)mathrmd p","category":"page"},{"location":"#Fully-discrete-numerical-scheme","page":"Home","title":"Fully discrete numerical scheme","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We assume that the computational domain is 0Ltimes-P_LP_R for x and v. The system is periodic in x-direction with period L and has compact support on -P_LP_R in v-direction. The mesh is as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"x_j=(j-1)Delta x j=1M Delta x=LM  (M text is odd)","category":"page"},{"location":"","page":"Home","title":"Home","text":"p_ell-1 = (ell-1)Delta p-P_L  ell=1N  Delta p=(P_L+P_L)N","category":"page"},{"location":"","page":"Home","title":"Home","text":"We use the spectral Fourier expansion to approximate E_x as it is periodic in the x-direction,","category":"page"},{"location":"","page":"Home","title":"Home","text":"E_xj=sum_k=-(M-1)2^(M-1)2 hatE_xk(t)e^2pi ijkM  text for  j=1M","category":"page"},{"location":"","page":"Home","title":"Home","text":"For the distribution functions (f_0 mathbf f), we use a spectral Fourier expansion for the x-direction and a finite-volume method for the p-direction. For simplicity, we only present the representation for f_0, the notations for mathbf f are the same. Here f_0jell(t) denotes the average of f_0(x_jpt) over a cell C_ell=p_ell-12 p_l+12 with the midpoint p_ell-12=(ell-12)Delta p-P_L, that is,","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_0jell(t)=frac1Delta p int_C_ell f_0(x_jpt)mathrmdp","category":"page"},{"location":"","page":"Home","title":"Home","text":"and also by Fourier expansion in x-direction, then","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_0jell(t)=sum_k=-(M-1)2^(M-1)2 hatf_0kell(t)e^2pi ijkM  j=1M","category":"page"},{"location":"","page":"Home","title":"Home","text":"To evaluate the value of f_0 off-grid in p-direction, we need to reconstruct a continuous function by using the cell average quantity f_0jell. ","category":"page"},{"location":"srs_without_spin/#SRS-without-spin","page":"Validation","title":"SRS without spin","text":"","category":"section"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"using Plots\nusing FFTW\nusing VectorSpinVlasovMaxwell1D1V","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"left\nbeginaligned\nfracpartial fpartial t + p fracpartial fpartial x +  E_x  - mathbf A_perp cdot fracpartial mathbf A_perppartial x fracpartial fpartial p = 0\nfracpartial E_xpartial t = -int_mathbbR p f  mathrmdp\nfracpartial E_ypartial t = - fracpartial^2 A_ypartial x^2 + A_y int_mathbbR  f  mathrmdp\nfracpartial E_zpartial t = - fracpartial^2 A_zpartial x^2 + A_z int_mathbbR  f  mathrmdp\nfracpartial mathbf A_perppartial t = - mathbf E_perp\nfracpartial E_xpartial x = int_mathbbR f mathrmdp - 1\nendaligned\nright","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"We consider the periodic condition with spatial period L=4pik_e, also take H=5 for the computational domain in v-direction.","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"Mathematical domain parameters are taken as N_x=129 N_v=129 Delta t =005","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"We take the following values for physical parameters:","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"alpha=002 k_e=12231 k_0=2k_e v_th=017","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"w_0=26428 k_s=k_e w_s=15799 w_e=10629 ","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"We use a perturbed Maxwellian as an initial condition for f","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"f(t=0xp)=(1+alpha cos(k_e x))frac1sqrt2piv_the^-fracp^22v_th^2","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"and the initial longitudinal electric field","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"E_x(t=0x)=(alpha k_e)sin(k_e x) ","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"Here alpha and k_e are the amplitude and the wave number of the perturbation respectively, and the v_th is the electron thermal speed. For the transverse fields, we consider an incident electromagnetic wave with circular polarization:","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"beginaligned\n E_y(t=0x)=E_0 cos(k_0 x) \n E_z(t=0x)=E_0 sin(k_0 x)\n A_y(t=0x)=-fracE_0w_0 sin(k_0 x)  \n A_z(t=0x)=fracE_0w_0 cos(k_0 x)\nendaligned ","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"where the k_0 and w_0 are the wave number and the amplitude of the transverse electric field respectively. ","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"We also take the amplitude of the incident wave E_ref=0325 as a reference value.  be in the range 025E_ref leq E_0  leq 2E_ref","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":"time evolution of the longitudinal electric field norm","category":"page"},{"location":"srs_without_spin/","page":"Validation","title":"Validation","text":" E_x (t) =left(frac12int_0^L E_x^2(tx) mathrmdmathrmxright )^frac12","category":"page"},{"location":"example/#Example","page":"Example","title":"Example","text":"","category":"section"},{"location":"example/","page":"Example","title":"Example","text":"using Plots\nusing FFTW\nusing VectorSpinVlasovMaxwell1D1V\n\nfunction run()\n\n    T = 50 # 4000  # final time\n    M = 65   # partition of x\n    N = 129   # partition of v\n    H = 5.0 / 2   # v domain size()\n    kkk = 1.2231333040331807  #ke\n    L = 4pi / kkk  # x domain size()\n    h = 0.04 #time step size()\n    NUM = floor(Int, T / h + 1.1) # time step number\n    a = 0.02 # 0.001; perturbation coefficient\n    h_int = 0.2 # hbar\n    k0 = 2.0 * kkk\n    ww = sqrt(1.0 + k0^2.0) # w0\n    ata = 0.2\n\n    mesh = Mesh(N, M, H, L)\n    adv  = Translator(mesh)\n    \n    E1, E2, E3, A2, A3 = initialfields( mesh, a, ww, kkk, k0)\n    f0, f1, f2, f3 = initialfunction(mesh, a, kkk, ata)\n\n    Ex_energy = Float64[]\n    E_energy = Float64[]\n    B_energy = Float64[]\n    energy = Float64[]\n    Sz = Float64[]\n    Tvalue = Vector{Float64}[]\n    time = Float64[]\n\n    results = Diagnostics(f0, f2, f3, E1, E2, E3, A2, A3, mesh, h_int)\n\n    H2fh = H2fhOperator(adv)\n    He = HeOperator(adv)\n    HAA = HAAOperator(adv)\n    H3fh = H3fhOperator(adv)\n    H1f = H1fOperator(adv)\n\n    for i = 1:NUM # Loop over time\n\n        step!(f0, f1, f2, f3, E3, A3, H2fh, h/2, h_int)\n        step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h/2)\n        step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h/2)\n        step!(f0, f1, f2, f3, E2, A2, H3fh, h/2, h_int)\n        step!(f0, f1, f2, f3, E1, H1f, h)\n        step!(f0, f1, f2, f3, E2, A2, H3fh, h/2, h_int)\n        step!(f0, f1, f2, f3, E2, E3, A2, A3, HAA, h/2)\n        step!(f0, f1, f2, f3, E1, E2, E3, A2, A3, He, h/2)\n        step!(f0, f1, f2, f3, E3, A3, H2fh, h/2, h_int)\n        \n        save!(results, i*h, f0, f2, f3, E1, E2, E3, A2, A3)\n\n    end\n\n    results\n\nend","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"data = run()","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(data.time, data.Ex_energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(data.time, data.E_energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(data.time, data.B_energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(data.time, data.energy)","category":"page"},{"location":"example/","page":"Example","title":"Example","text":"plot(data.time, data.Sz)","category":"page"}]
}
