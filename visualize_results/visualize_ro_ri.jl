using LESStudySetup
using CairoMakie
using SixelTerm
using Printf, Dates
using Oceananigans: compute!

using LESStudySetup.Diagnostics
using LESStudySetup.Diagnostics: load_snapshots,ζ,N²,M²
set_theme!(Theme(fontsize = 12))

function visualize_Ro_Ri(cooling, wind, dTf, idxes, ks)
    cooling = @sprintf("%03d", cooling)
    wind = replace("$(wind)","." => "" )
    if dTf < 0
        fileparams = "four_vortices_cooling_$(cooling)_wind_$(wind)"
    else
        if length(wind) < 2
            wind = "0" * wind
        end
        fileparams = "free_surface_short_test_$(cooling)_wind_$(wind)_dTf_$(dTf)"
    end
    filehead = "/orcd/data/abodner/001/simone/LESStudySetup.jl/experiments/freely_evolving_experiments/"
    filename = filehead * "hydrostatic_snapshots_" * fileparams * ".jld2"
    metadata = filehead * "experiment_" * fileparams * "_metadata.jld2"
    filesave = "/orcd/data/abodner/001/shirui/LESStudySetup/results/"

    # load all the data!!
    snapshots = load_snapshots(filename; metadata)

    # Let's pick the last snapshot!
    times = snapshots[:T].times
    snapshot_number = length(times) ÷ 2
    nday = @sprintf("%.0f", (times[snapshot_number])/60^2/24)
    println("Plotting timeseries $snapshot_number with $(nday) days of simulation...")
    t0 = now()

    # Plot the Ro and Ri
    ρ₀ = parameters.ρ₀
    f = parameters.f
    τw = parameters.τw
    θ = parameters.θ
    H0 = parameters.m₀
    uf2 = τw/ρ₀
    r = 16/7.5e-9 # EBF in W/m^2
    
    xT, yT, zT = nodes(snapshots[:T][snapshot_number])
    tns = 1:5:snapshot_number
    Mi = size(idxes,1)
    Ros = zeros(Mi,length(ks),length(tns))
    Ris = zeros(Mi,length(ks),length(tns))
    EBFs = zeros(Mi,length(tns))

    fig = Figure(size = (900, 800))
    fig2 = Figure(size = (900, 800))
    fig3 = Figure(size = (900, 300))
    c5 = [:red, :blue, :green, :purple, :orange]
    label5 = ["(25,75)","(50,75)", "(50,50)", "(50,25)","(75,75)"]

    for (i,n) in enumerate(tns)
        println("Computing Ro and Ri for snapshot $n, wall time: $((now() - t0).value/1e3) seconds.")

        Ro = compute!(Field(ζ(snapshots, n)/f))
        N2 = compute!(Field(N²(snapshots, n)))
        M2 = compute!(Field(M²(snapshots, n)))
        #Ri = compute!(Field(N2*f^2/M2^2))
        #EBF = compute!(Field(uf2*M2*cos(θ)/f*r))
        Ros[:,:,i] = stack([interior(Ro,idxes[m][1],idxes[m][2],ks) for m = 1:Mi],dims=1)
        println("Ros[:,:,i]")
        println(Ros[:,:,i])
        N2i = stack([interior(N2,idxes[m][1],idxes[m][2],ks) for m = 1:Mi],dims=1)
        M2i = stack([interior(M2,idxes[m][1],idxes[m][2],ks) for m = 1:Mi],dims=1)
        
        for m = 1:Mi,k = 1:length(ks)
            if M2i[m,k] != 0
                Ris[m,k,i] = f^2*abs(N2i[m,k])/M2i[m,k]^2
            end
        end
        println("Ris[:,:,i]")
        println(Ris[:,:,i])
        println("Lsi in m")
        println(2π*M2i*H0/(f^2).*sqrt.((1 .+ Ris[:,:,i])*2/5))
        println("τs in days")
        println(sqrt.(54/5*(1 .+ Ris[:,:,i]))/f/60^2/24)
        println("EBF[:,i]")
        EBFs[:,i] = uf2*cos(θ)/f*r*M2i[:,1]
        println(EBFs[:,i])

        # Plot Ros
        for k = 1:length(ks)
            ax = Axis(fig[k, 1]; title = "z = $(zT[ks[k]]) m", ylabel = L"\text{Ro}")
            for m = 1:Mi
                lines!(ax,times[tns]/60^2/24, Ros[m,k,:], color = c5[m], label = label5[m])
            end
            if k == length(ks)
                ax.xlabel = "time (days)"
            else
                hidexdecorations!(ax, ticks = false)
            end
            if k == 1
                axislegend(ax)
            end
        end
        save(filesave * "Ro_" * fileparams * "_d$(nday).pdf", fig)

        # Plot Ris
        for k = 1:length(ks)
            ax = Axis(fig2[k, 1]; title = "z = $(zT[ks[k]]) m", ylabel = L"\text{Ri}")
            for m = 1:Mi
                lines!(ax,times[tns]/60^2/24, Ris[m,k,:], color = c5[m], label = label5[m])
            end
            if k == length(ks)
                ax.xlabel = "time (days)"
            else
                hidexdecorations!(ax, ticks = false)
            end
            if k == 1
                axislegend(ax)
            end
        end
        save(filesave * "Ri_" * fileparams * "_d$(nday).pdf", fig2)

        # Plot EBFs
        ax = Axis(fig3[1, 1]; xlabel = "time (days)", ylabel = "EBF in W/m²")
        for m = 1:Mi
            lines!(ax,times[tns]/60^2/24, EBFs[m,:], color = c5[m], label = label5[m])
        end
        axislegend(ax)
        save(filesave * "EBF_" * fileparams * "_d$(nday).pdf", fig3)
    end

    return
end

cooling, wind, dTf = 50, 0.1, 2
idxes = [[250,750],[500,750],[500,500],[500,250],[750,750]]
ks = [124,113,95]
visualize_Ro_Ri(cooling, wind, dTf, idxes, ks)