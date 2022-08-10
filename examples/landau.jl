using FFTW
using LinearAlgebra
using Plots
using ProgressMeter
using Statistics
using VectorSpinVlasovMaxwell1D1V
import VectorSpinVlasovMaxwell1D1V: bspline

"""
    compute_rho(f)
Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho(mesh :: Mesh, f)
    
   dv =  mesh.dv
   rho = dv * sum(real(f), dims=2)
   vec(rho .- mean(rho)) 

end

"""
    compute_e(f)
compute Ex using that -ik*Ex = rho 
"""
function compute_e( mesh :: Mesh, f )

   rho = compute_rho( mesh, f )
   nx = mesh.nx
   xmin = mesh.xmin
   xmax = mesh.xmax
   k =  2π / (xmax - xmin)
   modes = k * vcat(0:div(nx,2)-1, -div(nx,2):-1)
   modes[1] = 1.0
   rhok = fft(rho) ./ modes .* (-1im)
   real(ifft!(rhok))

end


function landau()

    total_time = 50 
    nv = 128   
    nx = 128  
    vmin, vmax = -6., 6  
    ke =  0.5
    xmin, xmax = 0, 2pi / ke 
    dt = 0.1 
    nsteps = floor(Int, total_time / dt + 1.1) 
    α = 0.01 

    mesh = Mesh(xmin, xmax, nx, vmin, vmax, nv)
    mesh.v .= LinRange(-vmin, vmax, nv+1)[1:end-1]
    adv = BSplineAdvection(mesh)

    f = zeros(mesh.nv, mesh.nx)
    f .= (1 .+ α*cos.(ke*mesh.x'))/sqrt(2π) .* exp.(-0.5*mesh.v.^2)

    e = compute_e(mesh, f)
    energy = [0.5 * log(sum(e.^2))]
    time = [0.0]


    @showprogress 1 for i = 1:nsteps # Loop over time

        advection!(f', adv, mesh.v .* 0.5dt )
        e .= compute_e(mesh, f)
        advection!(f, adv, e .* dt )
        advection!(f', adv, mesh.v .* 0.5dt )
        push!(energy, 0.5 * log(sum(e.^2)))
        push!(time, i * dt)

    end

    time, energy

end

#time, energy = landau()
#
#plot(time, energy, label="electric energy")


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
      x  = range(xmin, stop=xmax, length=nx+1)[1:end-1]     
      new( xmin, xmax, nx, dx, x)
   end
end

"""
    compute_rho(meshv, f)

Compute charge density
ρ(x,t) = ∫ f(x,v,t) dv
"""
function compute_rho(meshv::UniformMesh,
        f::Array{Complex{Float64},2})

   dv = meshv.dx
   rho = dv * sum(real(f), dims=2)
   vec(rho .- mean(rho)) # vec squeezes the 2d array returned by sum function
end

"""
    compute_e(meshx, rho)
compute Ex using that -ik*Ex = rho
"""
function compute_e(meshx::UniformMesh, rho::Vector{Float64})
   nx = meshx.nx
   k =  2π / (meshx.xmax - meshx.xmin)
   modes = zeros(Float64, nx)
   modes .= k * vcat(0:div(nx,2)-1,-div(nx,2):-1)
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
    mesh     :: UniformMesh
    modes    :: Vector{Float64}
    eig_bspl :: Vector{Float64}
    eigalpha :: Vector{Complex{Float64}}
    
    function Advection( p, mesh )
        nx        = mesh.nx
        modes     = zeros(Float64, nx)
        modes    .= [2π * i / nx for i in 0:nx-1]
        eig_bspl  = zeros(Float64, nx)
        eig_bspl  = zeros(Float64, nx)
        eig_bspl .= bspline(p, -div(p+1,2), 0.0)
        for i in 1:div(p+1,2)-1
            eig_bspl .+= bspline(p, i-(p+1)÷2, 0.0) * 2 .* cos.(i * modes)
        end
        eigalpha  = zeros(Complex{Float64}, nx)
        new( p, mesh, modes, eig_bspl, eigalpha )
    end
    
end

# + slideshow={"slide_type": "slide"}
function (adv :: Advection)(f    :: Array{Complex{Float64},2}, 
                            v    :: Vector{Float64}, 
                            dt   :: Float64)
    
   nx = adv.mesh.nx
   nv = length(v)
   dx = adv.mesh.dx
    
   fft!(f,1)
    
   @inbounds for j in 1:nv
      alpha = dt * v[j] / dx
      # compute eigenvalues of cubic splines evaluated at displaced points
      ishift = floor(-alpha)
      beta   = -ishift - alpha
      fill!(adv.eigalpha,0.0im)
      for i in -div(adv.p-1,2):div(adv.p+1,2)
         adv.eigalpha .+= (bspline(adv.p, i-div(adv.p+1,2), beta) 
                        .* exp.((ishift+i) * 1im .* adv.modes))
      end
          
      # compute interpolating spline using fft and properties of circulant matrices
      
      f[:,j] .*= adv.eigalpha ./ adv.eig_bspl
        
   end
        
   ifft!(f,1)
    
end            

function landau( ϵ, kx, meshx, meshv)
    nx = meshx.nx
    nv = meshv.nx
    x  = meshx.x
    v  = meshv.x
    f  = zeros(Complex{Float64},(nx,nv))
    f .= (1.0.+ϵ*cos.(kx*x))/sqrt(2π) .* transpose(exp.(-0.5*v.^2))
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
  x, v = meshx.x, meshv.x    
  dx = meshx.dx
  
  # Set distribution function for Landau damping
  ϵ, kx = 0.001, 0.5
  f = landau( ϵ, kx, meshx, meshv)
  fᵗ = zeros(Complex{Float64},(nv,nx))
    
  # Instantiate advection functions
  advection_x! = Advection(p, meshx)
  advection_v! = Advection(p, meshv)
  
  # Set time step
  dt = tf / nt
  
  # Run simulation
  ℰ = Float64[]
  
  for it in 1:nt
        
       advection_x!(f, v, 0.5dt)

       ρ = compute_rho(meshv, f)
       e = compute_e(meshx, ρ)
        
       push!(ℰ, 0.5*log(sum(e.*e)*dx))
        
       transpose!(fᵗ, f)
       advection_v!(fᵗ, e, dt)
       transpose!(f, fᵗ)
    
       advection_x!(f, v, 0.5dt)
        
  end
                  
  ℰ

end

nt = 1000
tf = 100.0
t  = range(0.0, stop=tf, length=nt)
@time nrj = landau_damping(tf, nt);

plot( t, nrj; label = "E")
plot!(t, -0.1533*t.-5.50; label="-0.1533t.-5.5")

