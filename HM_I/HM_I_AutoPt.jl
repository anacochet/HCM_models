using SpecialFunctions, LinearAlgebra
using Plots, Printf, PlotThemes

#= 
This code solves the PDEs of Model I from Schmalholz et al., 2024:
Equations for porosity (ϕ) and fluid pressure (Pf) diffusion
Without (de)hydratation reactions
Including poroelastic effects

PDEs:
∂Pf/∂t = ∂/∂x (k / (Bϕ * ηf) * ∂Pf/∂x)   [1] Eq. 5 in Schmalholz et al., 2024
=#

function HM_I()

    # Which dimension do you want ? 
    # -----------------------------
    Dim     = 3
    
    # [For 2D or 3D] Which anomaly do you want ? 
    # ------------------------------------------
    layer       = 0 
    gaussienne  = 1
    ellipse     = 0

    # Calcul of error with the analytical solution for a gaussienne
    ana_sol     = 1


    # I - Initialisalisation
    # ----------------------
 
    # Physics
    ηf_phys     = 1000                              # Viscosity of fluid            | [Pa.s]
    k_phys      = 0.5                               # Permeability                  | [m2]
    βϕ_phys     = 0.1                               # Pore compressibility          | [Pa-1]
    P_BG_phys   = 1.6e9                             # Pressure background           | [Pa]
    Pf_max_phys = 0.8e9                             # Maximum of Pf anomaly         | [Pa]
    dPf_phys    = 0.1                               # Magnitude of Pf pertubation = standard deviation of the gaussienne | [m]
    # Spacial
    L_phys      = 1                                 # Length of the box             | [m]
    xo          = 0.0                               # localisation of the center of the anomaly on x
    yo          = 0.0                               # localisation of the center of the anomaly on y
    zo          = 0.0                               # localisation of anomaly on z

    # Numerics 
    Nx          = 101
    Ny          = 101
    Nz          = 101
    nt          = 10
    niter       = 1e5
    tol         = 1e-5
    nout        = 100
    CFL         = 0.98
    c_fact      = 0.5


    #---------------------------------------
    
    car = (
        Lc = L_phys, 
        Pc = P_BG_phys, 
        ηc = ηf_phys, 
        tc = ηf_phys/P_BG_phys
    )

    physics = (
        L      = (x = L_phys / car.Lc, y = L_phys / car.Lc, z = L_phys / car.Lc),
        P_BG   = P_BG_phys / car.Pc,
        ηf     = ηf_phys   / car.ηc,
        Βϕ     = βϕ_phys   * car.Pc,
        k      = k_phys    / car.Lc^2, 
        dPf    = dPf_phys  / car.Lc, 
        Pf_max = Pf_max_phys / car.Pc,
        xo     = xo,
        yo     = yo,
        zo     = zo
    )

    κ = physics.k / physics.Βϕ / physics.ηf                                   # Diffusivity coefficient

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

    ano_init    = (
        layer = layer, 
        gaussienne = gaussienne, 
        ellipse = ellipse
    )

    num = (
        CFL = CFL, 
        c_fact = c_fact,
        nout = nout,
        tol = tol 
    )


    Dim == 0 ? error("0 dimention = Nan diffusion ") : nothing
    Dim > 3 ? error("On s'arrete à la troisieme dimention d'espace ici ! XD") : nothing
    (ana_sol == 1) && (ana_sol + gaussienne != 2) ? error("La solution analytique doit être initiée pour une guassienne") : nothing

    theme(:dark::Symbol;)


    # Calcul in function of the dimension
    Dim == 1 ? Dim1(;physics,κ,n,Δ,Δt,mesh,num) : nothing
    Dim == 2 ? Dim2(;physics,κ,n,Δ,Δt,mesh,num, ano_init, ana_sol) : nothing
    Dim == 3 ? Dim3(;physics,κ,n,Δ,Δt,mesh,num, ano_init, ana_sol) : nothing

end


function Dim1(;physics,κ,n,Δ,Δt,mesh,num)

    # Allocate array 
    Pf      = zeros(n.c.x)
    Pf0     = zero(Pf)
    Pf_init = zero(Pf)
    ∇_qPf   = zero(Pf)
    R_Pf    = zero(Pf)
    R_Pf0   = zero(Pf)
    ∂Pf∂τ   = zero(Pf)

    qPf     = zeros(n.v.x)
    errPf   = 0


    # II - Initialization
    # -------------------
    
    # We make a gaussienne 
    Pf      = physics.P_BG * ones(n.c.x)
    Pf     .= Pf + physics.Pf_max * exp.(-mesh.c.x.^2/2/physics.dPf^2)
    Pf_init = copy(Pf)

    # Calcul du coeff lambda
    λmax = 2 * abs(κ ./ Δ.x .^ 2) + abs(2 * κ ./ Δ.x .^ 2 + 1 ./ Δt)

    # Pseudo time step is proportional to λmax (e.g., Oakley, 1995)
    Δτ    = 2 / sqrt(maximum(λmax)) * num.CFL


    # III - Computation
    # -----------------

    # Time loop 
    for it = 1: n.t

        @printf("Time step = %06d\n", it)

        # update
        Pf0 .= Pf

        # Damping coefficient is proportional to λmin (e.g., Oakley, 1995)
        λmin  = 0.          # wild guess  
        c     = 2*sqrt.(λmin)*num.c_fact
        β     = 2 .* Δτ ./ (2 .+ c.*Δτ)
        α     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)

        #Pseudo-time step 
        for iter=1:n.iter

            # Update old residual 
            R_Pf0   .= R_Pf

            # Limits -> flux = 0
            qPf[1]        = 0
            qPf[end]      = 0

            # Flux
            qPf[2:end-1] .= - physics.k .* diff(Pf) ./Δ.x

            # ∇q
            ∇_qPf        .= diff(qPf) ./ Δ.x 

            # Residual
            R_Pf         .= - ((physics.Βϕ * physics.ηf) * (Pf - Pf0) / Δt + ∇_qPf )

            # Pseudo-rate update
            ∂Pf∂τ       .= β * R_Pf .+ α * ∂Pf∂τ

            # Solution update
            Pf          .+= Δτ/(physics.ηf*physics.Βϕ) * ∂Pf∂τ
            
            if mod(iter, num.nout)==0 || iter==1
                # Error calc 
                iter == 1 ? errPf = norm(R_Pf)/sqrt(length(R_Pf)) : nothing
                Abs_errPf = norm(R_Pf)/sqrt(length(R_Pf))
                Rel_errPf = Abs_errPf / errPf
                
                @printf("iter = %06d:  Abs. Err. Pf = %2.1e Abs. Rel. Err. Pf = %2.1e\n", iter, Abs_errPf, Rel_errPf)

                Rel_errPf < num.tol && Abs_errPf < num.tol ? break : nothing

                if isnan(Abs_errPf)
                    error("Residual errPf = NaN!")
                end

                # Mettre a jour λmax ici si la le coeff est non constant dans le temps

                # Rayleigh quotient (e.g., Joldes et al., 2011)
                λmin  = abs.((sum(Δτ*∂Pf∂τ.*(R_Pf .- R_Pf0)))) / sum( (Δτ*∂Pf∂τ) .* (Δτ*∂Pf∂τ) )

                # Dynamic evaluation of PT iteration parameters
                c     = 2*sqrt.(λmin)*num.c_fact  
                β     = 2 .* Δτ ./ (2 .+ c.*Δτ)
                α     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)
            end      

        end
    
        # IV - Visualization 
        # ------------------

        # Analytical solution 
        t       = it * Δt 
        Pfana   = zeros(n.c.x)
        Pfana   = @. physics.P_BG + physics.Pf_max / sqrt(1 + 2 * t * κ / physics.dPf^2) * exp(- mesh.c.x^2. / (2 * (physics.dPf^2 + 2 * t * κ)) )

        # Plots
        p1 = plot(mesh.c.x, Pf_init,label="Pf init", title = "Pf", xlabel="x", ylabel="Pf")
        p1 = scatter!(mesh.c.x, Pf0, label="Pf0", markercolor=:white, markerstrokecolor=:black)
        p1 = scatter!(mesh.c.x, Pf, label="Pf",  markercolor=:white, markerstrokecolor=:blue)
        p1 = plot!(mesh.c.x, Pfana, label="Pf ana")
        
        display(plot(p1, ratio_aspect=1))
        sleep(0.1)
    end 
end

function Dim2(;physics,κ,n,Δ,Δt,mesh,num,ano_init,ana_sol)

    # Allocate array 
    Pf      = zeros(n.c.x,n.c.y)
    Pf0     = zero(Pf)
    Pf_init = zero(Pf)
    ∇_qPf   = zero(Pf)
    R_Pf    = zero(Pf)
    R_Pf0   = zero(Pf)
    Pfana   = zero(Pf)
    ∂Pf∂τ   = zero(Pf)  

    qPf     = (x = zeros(n.v.x,n.c.y), y = zeros(n.c.x,n.v.y))
    errPf   = 0

    # II - Initialization
    # -------------------
    
    # We make a gaussienne 
    if ano_init.gaussienne == 1
        Pf      = physics.P_BG * ones(n.c.x, n.c.y)
        for i=1:n.c.x
            for j=1:n.c.y
                Pf[i,j]   = Pf[i,j] + physics.Pf_max*exp(-((mesh.c.y[j]-physics.yo)^2 + (mesh.c.x[i]-physics.xo)^2)/(physics.dPf^2)) 
            end
        end
    end 

    Pf_init  = copy(Pf)

    # Calcul du coeff lambda
    λmax = 2 * abs(κ ./ Δ.x .^ 2) + 2 * abs(κ ./ Δ.y .^ 2) + abs(2 * κ ./ Δ.y .^ 2 + 2 * κ ./ Δ.x .^ 2 + 1 ./ Δt)

    # Pseudo time step is proportional to λmax (e.g., Oakley, 1995)
    Δτ    = 2 / sqrt(maximum(λmax)) * num.CFL


    # III - Computation
    # -----------------

    # Time loop 
    for it = 1:n.t

        @printf("Time step = %06d\n", it)

        # update
        Pf0 .= Pf

        # Damping coefficient is proportional to λmin (e.g., Oakley, 1995)
        λmin  = 0.          # wild guess  
        c     = 2*sqrt.(λmin)*num.c_fact
        β     = 2 .* Δτ ./ (2 .+ c.*Δτ)
        α     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)

        #Pseudo-time step 
        for iter=1:n.iter

            # Update old residual 
            R_Pf0   .= R_Pf

            # Limits -> flux = 0
            qPf.x[1,:] .= 0 ; qPf.x[end,:] .= 0
            qPf.y[:,1] .= 0 ; qPf.y[:,end] .= 0

            # Flux
            @. qPf.x[2:end-1,:] = - physics.k * (Pf[2:end,:] - Pf[1:end-1,:]) /Δ.x
            @. qPf.y[:,2:end-1] = - physics.k * (Pf[:,2:end] - Pf[:,1:end-1]) /Δ.y

            # ∇q
            @. ∇_qPf      = (qPf.x[2:end,:] - qPf.x[1:end-1,:]) / Δ.x + (qPf.y[:,2:end] - qPf.y[:,1:end-1]) / Δ.y

            # Residual
            @. R_Pf       = - ((physics.Βϕ * physics.ηf) * (Pf - Pf0) / Δt + ∇_qPf)

            # Pseudo rate update
            @. ∂Pf∂τ      = β * R_Pf .+ α * ∂Pf∂τ            
            
            # Solution update
            @. Pf        += Δτ/(physics.ηf*physics.Βϕ) * ∂Pf∂τ
            
            if mod(iter, num.nout)==0 || iter==1
                # Error calc 
                iter == 1 ? errPf = norm(R_Pf)/sqrt(length(R_Pf)) : nothing
                Abs_errPf = norm(R_Pf)/sqrt(length(R_Pf))
                Rel_errPf = Abs_errPf / errPf
                
                @printf("iter = %06d:  Abs. Err. Pf = %2.1e Abs. Rel. Err. Pf = %2.1e\n", iter, Abs_errPf, Rel_errPf)

                Rel_errPf < num.tol && Abs_errPf < num.tol ? break : nothing

                if isnan(Abs_errPf)
                    error("Residual errVsx = NaN!")
                end

                # Mettre a jour λmax ici si la le coeff est non constant dans le temps

                # Rayleigh quotient (e.g., Joldes et al., 2011)
                λmin  = abs.((sum(Δτ*∂Pf∂τ.*(R_Pf .- R_Pf0)))) / sum( (Δτ*∂Pf∂τ) .* (Δτ*∂Pf∂τ) )

                # Dynamic evaluation of PT iteration parameters
                c     = 2*sqrt.(λmin)*num.c_fact  
                β     = 2 .* Δτ ./ (2 .+ c.*Δτ)
                α     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)
            end       

        end
    

        # IV - Visualization 
        # ------------------

        if ana_sol == 1
            # Analytical solution 
            t       = it * Δt 
            Pfana   = physics.P_BG*ones(n.c.x, n.c.y) 
            for i=1:n.c.x 
                for j=1:n.c.y
                    Pfana[i,j] = Pfana[i,j] + physics.Pf_max / (1 + 4 * t * κ / ((physics.dPf)^2)) * exp(-(mesh.c.x[i]^2 + mesh.c.y[j]^2) / ((physics.dPf)^2 + 4 * t * κ)) 
                end
            end
            erreur_max = round(maximum(Pf.-Pfana), digits = 3)
        end


        # Plots
        p1 = heatmap(mesh.c.x, mesh.c.y, Pf', clims=(1,1.5), title="Step $it")
        # p1 = plot!(mesh.c.x, mesh.c.y[Int64(round(n.c.y/2))], lc=:red, label="section on x")
        # p1 = plot!(mesh.c.x[Int64(round(n.c.x/2))], mesh.c.y, lc=:green, label="section on y")
        if ana_sol == 1
            p2 = heatmap(mesh.c.x, mesh.c.y, (Pf.-Pfana)', title="max err. $erreur_max", colormap = :roma)
            display(plot(p1,p2, layout=(1,2), aspect_ratio=1))
        else 
            display(plot(p1, aspect_ratio=1))
        end

        if ana_sol == 1 
            p3 = plot(mesh.c.x, Pf[Int64(round(n.c.x/2)),:],label="Pf", title = "horizontally", xlabel="x", ylabel="Pf", lc=:red, ylims=(1.0,1.5))
            p3 = scatter!(mesh.c.x, Pf0[Int64(round(n.c.x/2)),:] , label="Pf0", ms=2.5, marker=:circle, markercolor=:white , markeralpha =:0.5, markerstrokecolor=:orange)
            p3 = plot!(mesh.c.x, Pfana[Int64(round(n.c.x/2)),:], label="Pf ana", lc=:white)

            p4 = plot(mesh.c.y, Pf[:,Int64(round(n.c.x/2))], label="Pf", title = "vertically", xlabel="y", ylabel="Pf", lc=:green, ylims=(1.0,1.5))
            p4 = scatter!(mesh.c.y, Pf0[:,Int64(round(n.c.x/2))], label="Pf0", ms=2.5, marker=:circle, markercolor=:white , markeralpha =:0.5, markerstrokecolor=:blue)
            p4 = plot!(mesh.c.y, Pfana[:,Int64(round(n.c.x/2))], label="Pf ana", lc=:white)

            if it == n.t
                display(plot(p3,p4, layout=(1,2)))
            end
        end
        
    end 
end

function Dim3(;physics,κ,n,Δ,Δt,mesh,num,ano_init,ana_sol)

    # Allocate array 
    Pf      = zeros(n.c...)
    Pf0     = zero(Pf)
    Pf_init = zero(Pf)
    ∇_qPf   = zero(Pf)
    R_Pf    = zero(Pf)
    R_Pf0   = zero(Pf)
    Pfana   = zero(Pf)
    ∂Pf∂τ   = zero(Pf)  

    qPf     = (x = zeros(n.v.x,n.c.y,n.c.z), y = zeros(n.c.x,n.v.y,n.c.z), z = zeros(n.c.x,n.c.y,n.v.z))
    errPf   = 0


    # II - Initialization
    # -------------------
    
    # We make a gaussienne 
    if ano_init.gaussienne == 1
        Pf      = physics.P_BG * ones(n.c...)
        for i=1:n.c.x
            for j=1:n.c.y
                for k=1:n.c.z
                    Pf[i,j,k]   = Pf[i,j,k] + physics.Pf_max * exp(-((mesh.c.z[k]-physics.zo)^2 + (mesh.c.y[j]-physics.yo)^2 + (mesh.c.x[i]-physics.xo)^2) / (physics.dPf^2)) 
                end    
            end  
        end 
    end

    Pf_init  = copy(Pf)

    # Calcul du coeff lambda
    λmax = 2 * abs(κ ./ Δ.x .^ 2) + 2 * abs(κ ./ Δ.y .^ 2) + 2 * abs(κ ./ Δ.z .^ 2) + abs(2 * κ ./ Δ.z .^ 2 + 2 * κ ./ Δ.y .^ 2 + 2 * κ ./ Δ.x .^ 2 + 1 ./ Δt)

    # Pseudo time step is proportional to λmax (e.g., Oakley, 1995)
    Δτ    = 2 / sqrt(maximum(λmax)) * num.CFL


    # III - Computation
    # -----------------

    # Time loop 
    for it = 1:n.t

        @printf("Time step = %06d\n", it)

        # Update
        Pf0 .= Pf

        # Damping coefficient is proportional to λmin (e.g., Oakley, 1995)
        λmin  = 0.          # wild guess  
        c     = 2*sqrt.(λmin)*num.c_fact
        β     = 2 .* Δτ ./ (2 .+ c.*Δτ)
        α     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)

        #Pseudo-time step 
        for iter=1:n.iter

            # Update old residual 
            R_Pf0   .= R_Pf

            # Limits -> flux = 0
            qPf.x[1,:,:] .= 0 ; qPf.x[end,:,:] .= 0
            qPf.y[:,1,:] .= 0 ; qPf.y[:,end,:] .= 0
            qPf.z[:,:,1] .= 0 ; qPf.y[:,:,end] .= 0

            # Flux
            @. qPf.x[2:end-1,:,:] = - physics.k * (Pf[2:end,:,:] - Pf[1:end-1,:,:]) /Δ.x
            @. qPf.y[:,2:end-1,:] = - physics.k * (Pf[:,2:end,:] - Pf[:,1:end-1,:]) /Δ.y
            @. qPf.z[:,:,2:end-1] = - physics.k * (Pf[:,:,2:end] - Pf[:,:,1:end-1]) /Δ.z

            # ∇q
            @. ∇_qPf      = (qPf.x[2:end,:,:] - qPf.x[1:end-1,:,:]) / Δ.x + (qPf.y[:,2:end,:] - qPf.y[:,1:end-1,:]) / Δ.y + (qPf.z[:,:,2:end] - qPf.z[:,:,1:end-1]) / Δ.z

            # Residual
            @. R_Pf       = - ((physics.Βϕ * physics.ηf) * (Pf - Pf0) / Δt + ∇_qPf)

            # Pseudo rate update
            @. ∂Pf∂τ      =β * R_Pf .+ α * ∂Pf∂τ       
            
            # Solution update
            @. Pf        += Δτ/(physics.ηf*physics.Βϕ) * ∂Pf∂τ
            
            if mod(iter, num.nout)==0 || iter==1
                # Error calc 
                iter == 1 ? errPf = norm(R_Pf)/sqrt(length(R_Pf)) : nothing
                Abs_errPf = norm(R_Pf)/sqrt(length(R_Pf))
                Rel_errPf = Abs_errPf / errPf
                
                @printf("iter = %06d:  Abs. Err. Pf = %2.1e Abs. Rel. Err. Pf = %2.1e\n", iter, Abs_errPf, Rel_errPf)

                Rel_errPf < num.tol && Abs_errPf < num.tol ? break : nothing

                if isnan(Abs_errPf)
                    error("Residual errVsx = NaN!")
                end

                # Mettre a jour λmax ici si la le coeff est non constant dans le temps

                # Rayleigh quotient (e.g., Joldes et al., 2011)
                λmin  = abs.((sum(Δτ*∂Pf∂τ.*(R_Pf .- R_Pf0)))) / sum( (Δτ*∂Pf∂τ) .* (Δτ*∂Pf∂τ) )

                # Dynamic evaluation of PT iteration parameters
                c     = 2*sqrt.(λmin)*num.c_fact  
                β     = 2 .* Δτ ./ (2 .+ c.*Δτ)
                α     = (2 .- c.*Δτ) ./ (2 .+ c.*Δτ)
            end       

        end
    

        # IV - Visualization 
        # ------------------

        # Analytical solution 
        if ana_sol == 1 
            t       = it * Δt 
            Pfana   = physics.P_BG*ones(n.c.x, n.c.y, n.c.z) 
            for i=1:n.c.x 
                for j=1:n.c.y
                    for k=1:n.c.z
                        Pfana[i,j,k] = Pfana[i,j,k] + physics.Pf_max / (1 + 4 * t * κ / ((physics.dPf)^2))^(3/2) * exp(-(mesh.c.x[i]^2 + mesh.c.y[j]^2 + mesh.c.z[k]^2) / ((physics.dPf)^2 + 4 * t * κ)) 
                    end 
                end 
            end

            erreur_max_x = round(maximum(Pf[:,:,Int64(round(n.c.z/2))].-Pfana[:,:,Int64(round(n.c.z/2))]), digits = 4)
            erreur_max_y = round(maximum(Pf[Int64(round(n.c.x/2)),:,:].-Pfana[Int64(round(n.c.x/2)),:,:]), digits = 4)
            erreur_max_z = round(maximum(Pf[:,Int64(round(n.c.y/2)),:].-Pfana[:,Int64(round(n.c.y/2)),:]), digits = 4)

        end 

        # Plots
        p01 = heatmap(mesh.c.x, mesh.c.y, Pf_init[:,:,Int64(round(n.c.z/2))]', xlabel = 'x', ylabel='y', title="initial condition")
        p02 = heatmap(mesh.c.y, mesh.c.z, Pf_init[Int64(round(n.c.x/2)),:,:]', xlabel = 'y', ylabel='z')
        p03 = heatmap(mesh.c.x, mesh.c.z, Pf_init[:,Int64(round(n.c.y/2)),:]', xlabel = 'x', ylabel='z')    

        p11 = heatmap(mesh.c.x, mesh.c.y, Pf[:,:,Int64(round(n.c.z/2))]',  xlabel = 'x', ylabel='y', clims=(1,1.5), title="Step $it")
        p12 = heatmap(mesh.c.y, mesh.c.z, Pf[Int64(round(n.c.x/2)),:,:]', xlabel = 'y', ylabel='z', clims=(1,1.5))
        p13 = heatmap(mesh.c.x, mesh.c.z, Pf[:,Int64(round(n.c.y/2)),:]', xlabel = 'x', ylabel='z', clims=(1,1.5))    
        
        if it ==1 
            display(plot(p01, p02, p03))
        end

        display(plot(p11, p12, p13))

        if ana_sol == 1 
            p21 = heatmap(mesh.c.x, mesh.c.y, (Pf[:,:,Int64(round(n.c.z/2))].- Pfana[:,:,Int64(round(n.c.z/2))])', title="err. $erreur_max_z", colormap = :roma)
            p22 = heatmap(mesh.c.y, mesh.c.z, (Pf[Int64(round(n.c.x/2)),:,:].- Pfana[Int64(round(n.c.x/2)),:,:])', title="err. $erreur_max_x", colormap = :roma)
            p23 = heatmap(mesh.c.x, mesh.c.z, (Pf[:,Int64(round(n.c.y/2)),:].- Pfana[:,Int64(round(n.c.y/2)),:])', title="err. $erreur_max_y", colormap = :roma)
        
            p3 = plot(mesh.c.x, Pf[Int64(round(n.c.x/2)),:,Int64(round(n.c.z/2))],label="Pf", title = "horizontally", xlabel="x", ylabel="Pf", lc=:red, ylims=(1.0,1.5))
            p3 = scatter!(mesh.c.x, Pf0[Int64(round(n.c.x/2)),:,Int64(round(n.c.z/2))] , label="Pf0", ms=2.5, marker=:circle, markercolor=:white , markeralpha =:0.5, markerstrokecolor=:orange)
            p3 = plot!(mesh.c.x, Pfana[Int64(round(n.c.x/2)),:,Int64(round(n.c.z/2))], label="Pf ana", lc=:white)

            
            display(plot(p21, p22, p23))
        
            if it == n.t
                display(plot(p3))
            end
    
        end     
    
    end 
end

HM_I()