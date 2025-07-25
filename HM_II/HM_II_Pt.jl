using SpecialFunctions, LinearAlgebra
using Plots, Printf, PlotThemes

#=
This code implements Hydro-Mechanical model II : Nonlinear diffusion equation for Porosity from Schmalholz et al. 2024
Equation of porosity diffusion with porosity-dependent permeability 
Without (de)hydratation reactions
Including poroelastic effects

Equations used to find the PDE of porosity in this case : 
1- Fluid mass conservation (per unit volume)
  ∂/∂t(ϕρf) = ∂/∂x(vfϕρf)       [1]
2- Dracy's law
  ϕvf = k/ηf ∂Pf/dx             [2]
3- Pore compressibility equation (poroelasticity, Yarushina & Podladchikov, 2015)
  ∂ϕ/∂t = βϕ ∂Pf/∂t             [3]
4 - Permeability–porosity relationship (Kozeny–Carman formulation)
  k = kref (ϕ/ϕref)ⁿ            [7]

Resulting nonlinear porosity diffusion equation:
∂ϕ/∂t = ∂/∂x (kref/(Βϕ*ηf) (ϕ/ϕref)^n ∂ϕ/dx)   [8]
=#

v2c1D(A)     = 0.5*(A[2:end]+A[1:end-1])                   # vertice to centroid in x 

v2cx2D(A)    = 0.5*(A[2:end,:]+A[1:end-1,:])              # vertice to centroid in x 
v2cy2D(A)    = 0.5*(A[:,2:end]+A[:,1:end-1])              # vertice to centroid in y 
v2c2D(A)     = 0.5*(v2cx(A[:,2:end])+v2cx(A[:,1:end-1]))  # vertice to centroid

v2cx3D(A)    = 0.5*(A[2:end,:,:]+A[1:end-1,:,:])              # vertice to centroid in x 
v2cy3D(A)    = 0.5*(A[:,2:end,:]+A[:,1:end-1,:])              # vertice to centroid in y 
v2cz3D(A)    = 0.5*(A[:,:,2:end]+A[:,:,1:end-1])              # vertice to centroid in y 
v2c3D(A)     = 0.5*(v2cx3D(v2cy3D(A[:,:,2:end]))+v2cx3D(v2cy3D(A[:,:,1:end-1])))  # vertice to centroid

function HM_II()

    # Which dimension do you want ? 
    # -----------------------------
    Dim     = 3
    
    # [For 2D or 3D] Which anomaly do you want ? 
    # ------------------------------------------
    layer       = 0 
    gaussienne  = 1
    ellipse     = 0


    # I - Initialisalisation
    # ----------------------
 
    # Physics
    ηf_phys     = 1000                              # Viscosity of fluid            | [Pa.s]
    kref_phys   = 0.2                               # Reference permeability        | [m2]
    βϕ_phys     = 0.1                               # Pore compressibility          | [Pa-1]
    ϕ_amb       = 0.1                               # Ambiant porosity 
    ϕ_max       = 0.5                               # Maximum of porosity anomaly         -> For initial condition
    ϕref        = 0.2                               # reference porosity
    n_darcy     = 3                                 # Carman Kozeny exponent
    # Spacial
    L_phys      = 1                                 # Length of the box             | [m]
    xo          = 0.0                               # localisation of the center of the anomaly on x
    yo          = 0.0                               # localisation of the center of the anomaly on y
    zo          = 0.0                               # localisation of the center of the anomaly on z
    dϕ_phys     = 0.1                               # Magnitude of ϕ perturbation   | [m] -> For initial condition

    # Numerics 
    Nx          = 51
    Ny          = 51
    Nz          = 51
    nt          = 10
    tol         = 1e-5
    niter       = 1e5
    nout        = 1000

    #---------------------------------------
    
    car = (
        L  = L_phys,                         # | [m] 
        η  = ηf_phys,                        # | [Pa.s]
        βϕ = βϕ_phys,                        # | [Pa-1] 
        t  = ηf_phys*βϕ_phys                 # | [s]               
    )

    physics = (
        L       = (x = L_phys / car.L, y = L_phys / car.L, z = L_phys / car.L),
        ϕ_amb   = ϕ_amb,
        ϕ_max   = ϕ_max,
        ϕref    = ϕref,
        n_darcy = n_darcy,
        ηf      = ηf_phys   / car.η,
        βϕ      = βϕ_phys   / car.βϕ,
        kref    = kref_phys / car.L^2, 
        dϕ      = dϕ_phys   / car.L,
        xo      = xo,
        yo      = yo,
        zo      = zo
    )

    κ = physics.kref / physics.βϕ / physics.ηf                                   # Diffusivity coefficient

    # Numerics 
    n = (
        v = (x = Nx, y = Ny, z = Nz),           # nb of vertices
        c = (x = Nx-1, y = Ny-1, z = Nz-1),      # nb of centroids
        iter = niter,
        out = nout,
        t  = nt
    )

    Δ = (
        x = physics.L.x/n.c.x,
        y = physics.L.y/n.c.y,
        z = physics.L.z/n.c.z,
    )    
    
    Δt = minimum(Δ)^2/(2.1*κ)

    mesh = (                                                                           
        v     = (x=LinRange(-physics.L.x/2, physics.L.x/2, n.v.x), y=LinRange(-physics.L.y/2, physics.L.y/2, n.v.y), z=LinRange(-physics.L.z/2, physics.L.z/2, n.v.z)),                               # vertice array
        c     = (x=LinRange(-physics.L.x/2+Δ.x/2, physics.L.x/2-Δ.x/2, n.c.x), y=LinRange(-physics.L.y/2+Δ.y/2, physics.L.y/2-Δ.y/2, n.c.y), z=LinRange(-physics.L.z/2+Δ.z/2, physics.L.z/2-Δ.z/2, n.c.z))               # centroid array
    )

    num = (
        tol   = tol,
        θ     = 4.3 * (minimum(physics.L))/maximum(n.c),
        nout  = nout
    )

    ano_init    = (
        layer = layer, 
        gaussienne = gaussienne, 
        ellipse = ellipse
    )

    Dim == 0 ? error("0 dimention = Nan diffusion ") : nothing
    Dim > 3 ? error("On s'arrete à la troisieme dimention d'espace ici ! XD") : nothing

    theme(:dark::Symbol;)

    # Calcul in function of the dimension
    Dim == 1 ? Dim1(;physics,κ,n,Δ,Δt,mesh,num) : nothing
    Dim == 2 ? Dim2(;physics,κ,n,Δ,Δt,mesh,num, ano_init) : nothing
    Dim == 3 ? Dim3(;physics,κ,n,Δ,Δt,mesh,num, ano_init) : nothing

end

function Dim1(;physics,κ,n,Δ,Δt,mesh,num)

    # Allocate array 
    ϕ       = zeros(n.c.x)
    ϕ0      = zero(ϕ)
    ϕ_init  = zero(ϕ)
    ∇_qϕ    = zero(ϕ)
    R_ϕ     = zero(ϕ)
    ∂ϕ∂τ    = zero(ϕ)

    qϕ     = zeros(n.v.x)  
    coeff_diff  = zeros(n.c.x-1)


    # II - Initialization
    # -------------------
    
    # We make a gaussienne 
    ϕ       = physics.ϕ_amb * ones(n.c.x)
    ϕ      .= ϕ + physics.ϕ_max * exp.(-mesh.c.x.^2/2/physics.dϕ^2)
    ϕ_init  = copy(ϕ)

    Δτ = 0

    # III - Computation
    # -----------------

    # Time loop 
    for it = 1: n.t

        @printf("Time step = %06d\n", it)

        # update
        ϕ0 .= ϕ

        # Pseudo-time step calculation
        # Δτ     = 5e2 * Δt * Δ.x^2 / ( maximum(coeff_diff) * 4)      # Don't ask me why we use 5e2...🥲

        #Pseudo-time step 
        for iter=1:n.iter

            # Limits -> flux = 0
            qϕ[1]        = 0
            qϕ[end]      = 0

            # Diffusion factor calculation
            # We make interpolation to have a good indice number and to be in phase with qϕ
            coeff_diff .= κ .* (v2c1D(ϕ) ./ physics.ϕref).^ physics.n_darcy

            iter == 1 ? Δτ = 5e2 * Δt * Δ.x^2 / ( maximum(coeff_diff) * 4) : nothing
            
            # Flux
            qϕ[2:end-1] .= - coeff_diff .* diff(ϕ) ./Δ.x

            # ∇q
            ∇_qϕ        .= diff(qϕ) ./ Δ.x 

            # Residual
            R_ϕ         .= - ((ϕ  .- ϕ0)  ./ Δt + ∇_qϕ)

            # Pseudo rate update
            ∂ϕ∂τ        .= R_ϕ .+ (1.0-num.θ).*∂ϕ∂τ   

            # Solution update 
            ϕ          .+= Δτ  .* ∂ϕ∂τ 
            
            # Error calculation
            errϕ          = norm(R_ϕ) / length(R_ϕ)
            if mod(iter, 100)==0
                @printf("iter = %06d:  Err. ϕ = %1.6e\n", iter, errϕ)
            end 
            iter == n.iter ? error("That don't converge !") : nothing
            isnan(errϕ) ? error("Residual err = NaN!") : nothing
            errϕ<num.tol ? break : nothing

        end

    
        # IV - Visualization 
        # ------------------

        # Plots
        p1 = plot(mesh.c.x, ϕ_init,label="ϕ init", title = "ϕ", xlabel="x", ylabel="ϕ")
        p1 = scatter!(mesh.c.x, ϕ0, label="ϕ0", markercolor=:white, markerstrokecolor=:black)
        p1 = scatter!(mesh.c.x, ϕ, label="ϕ",  markercolor=:white, markerstrokecolor=:blue)
        display(plot(p1, ratio_aspect=1))
        sleep(0.1)
    end 
end

function Dim2(;physics,κ,n,Δ,Δt,mesh,num, ano_init)

    # Allocate array 
    ϕ       = zeros(n.c.x,n.c.y)
    ϕ0      = zero(ϕ)
    ϕ_init  = zero(ϕ)
    ∇_qϕ    = zero(ϕ)
    R_ϕ     = zero(ϕ)
    ∂ϕ∂τ    = zero(ϕ)  
    Δτ      = 0

    qϕ      = (x = zeros(n.v.x,n.c.y), y = zeros(n.c.x,n.v.y))
    coeff_diff  = (x = zeros(n.c.x-1,n.c.y), y = zeros(n.c.x,n.c.y-1))

    # II - Initialization
    # -------------------
    
    # We make a gaussienne 
    if ano_init.gaussienne == 1
        ϕ      = physics.ϕ_amb * ones(n.c.x, n.c.y)
        for i=1:n.c.x
            for j=1:n.c.y
                ϕ[i,j]   = ϕ[i,j] + physics.ϕ_max*exp(-((mesh.c.y[j]-physics.yo)^2 + (mesh.c.x[i]-physics.xo)^2)/(physics.dϕ^2)) 
            end
        end
    end 

    ϕ_init  = copy(ϕ)


    # III - Computation
    # -----------------

    # Time loop 
    for it = 1:n.t

        @printf("Time step = %06d\n", it)

        # update
        ϕ0 .= ϕ

        #Pseudo-time step 
        for iter=1:n.iter

            # Limits -> flux = 0
            qϕ.x[1,:] .= 0 ; qϕ.x[end,:] .= 0
            qϕ.y[:,1] .= 0 ; qϕ.y[:,end] .= 0

            # Diffusion factor calculation
            # We make interpolation to have a good indice number and to be in phase with qϕ
            coeff_diff.x .= κ .* (v2cx2D(ϕ) ./ physics.ϕref).^ physics.n_darcy
            coeff_diff.y .= κ .* (v2cy2D(ϕ) ./ physics.ϕref).^ physics.n_darcy

            iter == 1 ? Δτ = 1e1 * Δt * minimum(Δ)^2 / ( maximum(physics.kref/(physics.βϕ*physics.ηf) * (ϕ./physics.ϕref).^physics.n_darcy) * 4) : nothing

            # Flux
            @. qϕ.x[2:end-1,:] = - coeff_diff.x * (ϕ[2:end,:] - ϕ[1:end-1,:]) /Δ.x
            @. qϕ.y[:,2:end-1] = - coeff_diff.y * (ϕ[:,2:end] - ϕ[:,1:end-1]) /Δ.y

            # ∇q
            @. ∇_qϕ      = (qϕ.x[2:end,:] - qϕ.x[1:end-1,:]) / Δ.x + (qϕ.y[:,2:end] - qϕ.y[:,1:end-1]) / Δ.y

            # Residual
            @. R_ϕ       = - ((ϕ - ϕ0) / Δt + ∇_qϕ)

            # Pseudo rate update
            @. ∂ϕ∂τ      = R_ϕ + (1.0-num.θ) * ∂ϕ∂τ             
            
            # Solution update
            @. ϕ        += Δτ  * ∂ϕ∂τ  
            
            # Error calculation
            errϕ      = norm(R_ϕ) / length(R_ϕ)
            if mod(iter, num.nout)==0 || errϕ<num.tol
                @printf("iter = %06d:  Err. ϕ = %2.1e\n", iter, errϕ)
            end

            if isnan(errϕ)
                error("Residual err = NaN!")
            end
            errϕ<num.tol ? break : nothing
        end
    

        # IV - Visualization 
        # ------------------

        # Plots
        p1 = heatmap(mesh.c.x, mesh.c.y, ϕ', clims=(0.1,0.5), title="Step $it")
        # p1 = plot!(mesh.c.x, mesh.c.y[Int64(round(n.c.y/2))], lc=:red, label="section on x")
        # p1 = plot!(mesh.c.x[Int64(round(n.c.x/2))], mesh.c.y, lc=:green, label="section on y")
        display(plot(p1, aspect_ratio=1))

        p3 = plot(mesh.c.x, ϕ[Int64(round(n.c.x/2)),:],label="ϕ", title = "horizontally", xlabel="x", ylabel="ϕ", lc=:red, ylims=(0.1,0.5))
        p3 = scatter!(mesh.c.x, ϕ0[Int64(round(n.c.x/2)),:] , label="ϕ0", ms=2.5, marker=:circle, markercolor=:white , markeralpha =:0.5, markerstrokecolor=:orange)

        p4 = plot(mesh.c.y, ϕ[:,Int64(round(n.c.x/2))], label="ϕ", title = "vertically", xlabel="y", ylabel="ϕ", lc=:green, ylims=(0.1,0.5))
        p4 = scatter!(mesh.c.y, ϕ0[:,Int64(round(n.c.x/2))], label="ϕ0", ms=2.5, marker=:circle, markercolor=:white , markeralpha =:0.5, markerstrokecolor=:blue)

        if it == n.t
            display(plot(p3,p4, layout=(1,2)))
        end
    end 
end

function Dim3(;physics,κ,n,Δ,Δt,mesh,num, ano_init)

    # Allocate array 
    ϕ      = zeros(n.c...)
    ϕ0     = zero(ϕ)
    ϕ_init = zero(ϕ)
    ∇_qϕ   = zero(ϕ)
    R_ϕ    = zero(ϕ)
    ∂ϕ∂τ   = zero(ϕ)  
    Δτ     = 0

    qϕ     = (x = zeros(n.v.x,n.c.y,n.c.z), y = zeros(n.c.x,n.v.y,n.c.z), z = zeros(n.c.x,n.c.y,n.v.z))
    coeff_diff  = (x = zeros(n.c.x-1,n.c.y,n.c.z), y = zeros(n.c.x,n.c.y-1,n.c.z), z = zeros(n.c.x,n.c.y,n.c.z-1))


    # II - Initialization
    # -------------------
    
    # We make a gaussienne 
    if ano_init.gaussienne == 1
        ϕ      = physics.ϕ_amb * ones(n.c...)
        for i=1:n.c.x
            for j=1:n.c.y
                for k=1:n.c.z
                    ϕ[i,j,k]   = ϕ[i,j,k] + physics.ϕ_max * exp(-((mesh.c.z[k]-physics.zo)^2 + (mesh.c.y[j]-physics.yo)^2 + (mesh.c.x[i]-physics.xo)^2) / (physics.dϕ^2)) 
                end    
            end  
        end 
    end

    ϕ_init  = copy(ϕ)


    # III - Computation
    # -----------------

    # Time loop 
    for it = 1:n.t

        @printf("Time step = %06d\n", it)

        # Update
        ϕ0 .= ϕ

        #Pseudo-time step 
        for iter=1:n.iter

            # Limits -> flux = 0
            qϕ.x[1,:,:] .= 0 ; qϕ.x[end,:,:] .= 0
            qϕ.y[:,1,:] .= 0 ; qϕ.y[:,end,:] .= 0
            qϕ.z[:,:,1] .= 0 ; qϕ.y[:,:,end] .= 0

            # Diffusion factor calculation
            # We make interpolation to have a good indice number and to be in phase with qϕ
            coeff_diff.x .= κ .* (v2cx3D(ϕ) ./ physics.ϕref).^ physics.n_darcy
            coeff_diff.y .= κ .* (v2cy3D(ϕ) ./ physics.ϕref).^ physics.n_darcy
            coeff_diff.z .= κ .* (v2cz3D(ϕ) ./ physics.ϕref).^ physics.n_darcy

            iter == 1 ? Δτ = 1e1 * Δt * minimum(Δ)^2 / ( maximum(physics.kref/(physics.βϕ*physics.ηf) * (ϕ./physics.ϕref).^physics.n_darcy) * 4) : nothing

            # Flux
            @. qϕ.x[2:end-1,:,:] = - coeff_diff.x * (ϕ[2:end,:,:] - ϕ[1:end-1,:,:]) /Δ.x
            @. qϕ.y[:,2:end-1,:] = - coeff_diff.y * (ϕ[:,2:end,:] - ϕ[:,1:end-1,:]) /Δ.y
            @. qϕ.z[:,:,2:end-1] = - coeff_diff.z * (ϕ[:,:,2:end] - ϕ[:,:,1:end-1]) /Δ.z

            # ∇q
            @. ∇_qϕ      = (qϕ.x[2:end,:,:] - qϕ.x[1:end-1,:,:]) / Δ.x + (qϕ.y[:,2:end,:] - qϕ.y[:,1:end-1,:]) / Δ.y + (qϕ.z[:,:,2:end] - qϕ.z[:,:,1:end-1]) / Δ.z

            # Residual
            @. R_ϕ       = - ((ϕ - ϕ0) / Δt + ∇_qϕ)

            # Pseudo rate update
            @. ∂ϕ∂τ      = R_ϕ + (1.0-num.θ) * ∂ϕ∂τ             
            
            # Solution update
            @. ϕ        += Δτ  * ∂ϕ∂τ  
            
            # Error calculation
            errϕ         = norm(R_ϕ) / length(R_ϕ)
            if mod(iter, 100)==0 || errϕ<num.tol
                @printf("iter = %06d:  Err. ϕ = %2.1e\n", iter, errϕ)
            end

            if isnan(errϕ)
                error("Residual err = NaN!")
            end
            errϕ<num.tol ? break : nothing
        end
    

        # IV - Visualization 
        # ------------------

        # Plots
        p01 = heatmap(mesh.c.x, mesh.c.y, ϕ_init[:,:,Int64(round(n.c.z/2))]', xlabel = 'x', ylabel='y',clims=(0.1,0.5), title="initial condition")
        p02 = heatmap(mesh.c.y, mesh.c.z, ϕ_init[Int64(round(n.c.x/2)),:,:]', xlabel = 'y', ylabel='z',clims=(0.1,0.5))
        p03 = heatmap(mesh.c.x, mesh.c.z, ϕ_init[:,Int64(round(n.c.y/2)),:]', xlabel = 'x', ylabel='z',clims=(0.1,0.5))    

        p11 = heatmap(mesh.c.x, mesh.c.y, ϕ[:,:,Int64(round(n.c.z/2))]',  xlabel = 'x', ylabel='y', clims=(0.1,0.5), title="Step $it")
        p12 = heatmap(mesh.c.y, mesh.c.z, ϕ[Int64(round(n.c.x/2)),:,:]', xlabel = 'y', ylabel='z', clims=(0.1,0.5))
        p13 = heatmap(mesh.c.x, mesh.c.z, ϕ[:,Int64(round(n.c.y/2)),:]', xlabel = 'x', ylabel='z', clims=(0.1,0.5))    
        
        if it ==1 
            display(plot(p01, p02, p03))
        end

        display(plot(p11, p12, p13))

        p2 = plot(mesh.c.x, ϕ[Int64(round(n.c.x/2)),:,Int64(round(n.c.z/2))],label="ϕ", title = "horizontally", xlabel="x", ylabel="ϕ", lc=:red, ylims=(0.1,0.5))
        p2 = scatter!(mesh.c.x, ϕ0[Int64(round(n.c.x/2)),:,Int64(round(n.c.z/2))] , label="ϕ0", ms=2.5, marker=:circle, markercolor=:white , markeralpha =:0.5, markerstrokecolor=:orange)
        
        if it == n.t
            display(plot(p2))
        end
    end 
end

HM_II()