using CSV
using DataFrames
using Plots
using FFTW
using MAT
using ProgressMeter
using Random
using SpinGEMPIC
using Statistics
import GEMPIC: OneDGrid, Maxwell1DFEM
import GEMPIC: eval_uniform_periodic_spline_curve
using VectorSpinVlasovMaxwell1D1V

import VectorSpinVlasovMaxwell1D1V: initialfunction
import VectorSpinVlasovMaxwell1D1V: initialfields
import VectorSpinVlasovMaxwell1D1V: diagnostics
import VectorSpinVlasovMaxwell1D1V: H2fh!
import VectorSpinVlasovMaxwell1D1V: He!
import VectorSpinVlasovMaxwell1D1V: HAA!
import VectorSpinVlasovMaxwell1D1V: H3fh!
import VectorSpinVlasovMaxwell1D1V: H1f!

function compute_rho(mesh::Mesh, f::Array{Float64,2})

   dv = mesh.dv
   ρ = dv * sum(real(f), dims=1)
   vec(ρ .- mean(ρ)) # vec squeezes the 2d array returned by sum function
end

function compute_e(mesh::Mesh, ρ::Vector{Float64})
   nx = mesh.nx
   modes = mesh.kx
   modes[1] = 1.0
   ρ = fft(ρ)./modes
   ρk .*= -1im
   ifft!(ρk)
   real(ρk)
end


T = 1000 # 4000  # final time
nx = 128  # partition of x
nv = 128   # partition of v
vmin, vmax = -6, 6   # v domain size()
ke = 1.22 #ke
xmin, xmax = 0, 4pi / ke  # x domain size()
h = 0.05 #time step size()
nsteps = floor(Int, T / h + 1.1) # time step number
a = 0.02 # 0.001; perturbation coefficient
h_int = 0.2 # hbar
k0 = 2.0 * ke
ww = 2.63 #sqrt(1.0 + k0^2.0) # w0
ata = 0.2

mesh = Mesh(xmin, xmax, nx, vmin, vmax, nv)
#adv = BSplineAdvection(mesh)
adv = PSMAdvection(mesh)

E1, E2, E3, A2, A3 = initialfields( mesh, a, ww, ke, k0)
f0, f1, f2, f3 = initialfunction(mesh, a, ke, ata)

domain = [xmin, xmax, xmax - xmin]
n_particles = 1000000
mesh_1d = OneDGrid( xmin, xmax, nx)
spline_degree = 3
σ, μ = 0.17, 0.0
kx, α = ke, a

df = CosGaussian(kx, a, σ, μ)

rng = MersenneTwister(123)
mass, charge = 1.0, 1.0

particle_group = ParticleGroup( n_particles, mass, charge, 1)   
sample!(rng, particle_group, df, mesh_1d)
set_common_weight(particle_group, (1.0/n_particles))

kernel_smoother2 = ParticleMeshCoupling( mesh_1d, n_particles, spline_degree-2) 
kernel_smoother1 = ParticleMeshCoupling( mesh_1d, n_particles, spline_degree-1)    
kernel_smoother0 = ParticleMeshCoupling( mesh_1d, n_particles, spline_degree)

maxwell_solver = Maxwell1DFEM(mesh_1d, spline_degree)

ρ_scalar = zeros(nx)
efield_poisson = zeros(nx)

solve_poisson!( efield_poisson, particle_group, kernel_smoother0, maxwell_solver, ρ_scalar )

ρ_vector = compute_rho(mesh, f0)

p = plot(layout=(1,2))
xg = LinRange(xmin, xmax, nx)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, ρ_scalar .- mean(ρ_scalar))
plot!(p[1,1], xg, sval, title=:ρ, label="scalar")
plot!(p[1,1], xg, ρ_vector, title=:ρ, label="vector")
sval = eval_uniform_periodic_spline_curve(spline_degree-1, efield_poisson)
plot!(p[1,2], xg, sval, title=:ex, label="scalar" )
plot!(p[1,2], xg, real(ifft(E1)), title=:ex, label="vector" )

xp = view(particle_group.array, 1, :)
vp = view(particle_group.array, 2, :)
s1 = view(particle_group.array, 3, :)
s2 = view(particle_group.array, 4, :)
s3 = view(particle_group.array, 5, :)
wp = view(particle_group.array, 6, :)

p = plot(layout=(3,1))
histogram!(p[1,1], xp, weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[1,1], x-> (1+α*cos(kx*x))/(sqrt(2π)/kx), 0., 2π/kx, lab="")
histogram!(p[2,1], vp, weights=wp, normalize=true, bins = 100, lab = "")
plot!(p[2,1], v-> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")

