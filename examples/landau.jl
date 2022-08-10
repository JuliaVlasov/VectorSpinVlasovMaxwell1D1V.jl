using FFTW
using LinearAlgebra
using Plots
using ProgressMeter
using Statistics
using VectorSpinVlasovMaxwell1D1V
import VectorSpinVlasovMaxwell1D1V: bspline

"""
    UniformMesh(xmin, xmax, nx)

1D uniform mesh data
"""
struct UniformMesh
   xmin  :: Float64
   xmax  :: Float64
   nx    :: Int
   dx    :: Float64
   x     :: Vector{Float64}
   function UniformMesh(xmin, xmax, nx)
      dx = (xmax - xmin) / nx
      x  = LinRange(xmin, xmax, nx+1)[begin:end-1]     
      new( xmin, xmax, nx, dx, x)
   end
end

"""
    compute_rho(meshv, f)

Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho(mesh::Mesh, f::Array{Float64,2})

   dv = mesh.dv
   rho = dv * sum(real(f), dims=1)
   vec(rho .- mean(rho)) # vec squeezes the 2d array returned by sum function
end

"""
    compute_e(meshx, rho)
compute Ex using that -ik*Ex = rho
"""
function compute_e(mesh::Mesh, rho::Vector{Float64})
   nx = mesh.nx
   modes = mesh.kx
   modes[1] = 1.0
   rhok = fft(rho)./modes
   rhok .*= -1im
   ifft!(rhok)
   real(rhok)
end

"""
    Advection(f, p, mesh, v, nv, dt)

Advection type

"""
struct Advection
    
    p        :: Int64 
    step     :: Float64
    modes    :: Vector{Float64}
    eig_bspl :: Vector{Float64}
    eigalpha :: Vector{Complex{Float64}}
    
    function Advection( p, mesh :: UniformMesh )
        nx        = mesh.nx
        step      = mesh.dx
        modes     = zeros(Float64, nx)
        modes    .= [2π * i / nx for i in 0:nx-1]
        eig_bspl  = zeros(Float64, nx)
        eig_bspl .= bspline(p, -div(p+1,2), 0.0)
        for i in 1:div(p+1,2)-1
            eig_bspl .+= bspline(p, i-(p+1)÷2, 0.0) * 2 .* cos.(i * modes)
        end
        eigalpha  = zeros(Complex{Float64}, nx)
        new( p, step, modes, eig_bspl, eigalpha )
    end

    function Advection( p, mesh :: Mesh )
        nv        = mesh.nv
        step      = mesh.dv
        modes     = zeros(Float64, nv)
        modes    .= [2π * i / nv for i in 0:nx-1]
        eig_bspl  = zeros(Float64, nv)
        eig_bspl .= bspline(p, -div(p+1,2), 0.0)
        for i in 1:div(p+1,2)-1
            eig_bspl .+= bspline(p, i-(p+1)÷2, 0.0) * 2 .* cos.(i * modes)
        end
        eigalpha  = zeros(Complex{Float64}, nv)
        new( p, step, modes, eig_bspl, eigalpha )
    end
    
end

# + slideshow={"slide_type": "slide"}
function (adv :: Advection)(f    :: Array{Float64,2}, 
                            v    :: Vector{Float64}, 
                            dt   :: Float64)
    
   f̂ = fft(f,1)
    
   @inbounds for j in eachindex(v)
      alpha = dt * v[j] / adv.step
      # compute eigenvalues of cubic splines evaluated at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(adv.eigalpha,0.0im)
      for i in -div(adv.p-1,2):div(adv.p+1,2)
         adv.eigalpha .+= (bspline(adv.p, i-div(adv.p+1,2), beta) 
                        .* exp.((ishift+i) * 1im .* adv.modes))
      end
          
      # compute interpolating spline using fft and properties of circulant matrices
      
      f̂[:,j] .*= adv.eigalpha ./ adv.eig_bspl
        
   end
        
   f .= real(ifft(f̂,1))
    
end            

function landau( ϵ, kx, mesh)
    nx = mesh.nx
    nv = mesh.nv
    x  = mesh.x
    v  = mesh.v
    f  = zeros(nv,nx)
    f .= transpose(1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* exp.(-0.5*v.^2)
    f
end

function landau_damping(tf::Float64, nt::Int64)
    
  # Set grid
  p = 3
  nx, nv = 128, 256
  xmin, xmax = 0.0, 4π
  vmin, vmax = -6., 6.
  meshx = UniformMesh(xmin, xmax, nx)
  meshv = UniformMesh(vmin, vmax, nv)
  mesh = Mesh(xmin, xmax, nx, vmin, vmax, nv)
  mesh.v .= meshv.x
  x, v = mesh.x, mesh.v    
  dx = mesh.dx
  
  # Set distribution function for Landau damping
  ϵ, kx = 0.001, 0.5
  f = landau( ϵ, kx, mesh)
  fᵗ = zeros(nx,nv)
    
  # Instantiate advection functions
  advection_x! = Advection(p, meshx)
  advection_v! = Advection(p, meshv)
  
  # Set time step
  dt = tf / nt
  
  # Run simulation
  ℰ = Float64[]
  
  @showprogress 1 for it in 1:nt
        
       transpose!(fᵗ, f)
       advection_x!(fᵗ, v, 0.5dt)
       transpose!(f, fᵗ)

       ρ = compute_rho(mesh, f)
       e = compute_e(mesh, ρ)
        
       push!(ℰ, 0.5*log(sum(e.*e)*dx))
        
       advection_v!(f, e, dt)
    
       transpose!(fᵗ, f)
       advection_x!(fᵗ, v, 0.5dt)
       transpose!(f, fᵗ)
        
  end
                  
  ℰ

end

nt = 1000
tf = 100.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau_damping(tf, nt);

plot( t, nrj; label = "E")
plot!(t, -0.1533*t.-5.50; label="-0.1533t.-5.5")

