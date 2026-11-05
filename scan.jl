#!/usr/bin/env julia
################################################################################
# Scan of clusternoisr.jl over the cluster width and the signal over noise.    #
#                                                                              #
# The cluster width (Gauss sigma, in strips) and the peak strip over noise are #
# scanned; amp_signal is set such that the peak strip has the wanted signal    #
# over noise for each cluster width. Per cell, the core sigma, the mean and    #
# the outlier fraction of the three options are tabulated, and two figures     #
# are saved:                                                                   #
#  - clusternoisr_scan_n<100 n_threshold>.png: core sigma and outlier          #
#    fraction against the cluster width, one column per signal over noise      #
#  - clusternoisr_scurve_SNR<peak/noise>_n<100 n_threshold>.png: the mean      #
#    residual against the hit position within a strip (the "S-curve"), one     #
#    panel per cluster width. A step in it marks where an edge strip enters    #
#    or leaves the cluster.                                                    #
#                                                                              #
# usage:  julia scan.jl                                                        #
# AMP_NOISE and N_THRESHOLD are taken from the environment as in               #
# clusternoisr.jl; the other parameters of clusternoisr.jl are overridden by   #
# the grid below.                                                              #
#                                                                              #
# (C) 2026, Philipp Bielefeldt, Uni Bonn                                       #
# released under the terms of Mozilla Public License (see file)                #
################################################################################

include(joinpath(@__DIR__, "clusternoisr.jl"));


### scan grid ###

# cluster width (Gauss sigma) in strips
const widths = [0.25, 0.5, 1.0, 2.0, 4.0];
# peak strip amplitude over the mean noise rms
const peak_snrs = [4.0, 7.0, 15.0];
# events per cell (the statistical error of sigma is 1/sqrt(2N), 0.7 % here)
const events_per_cell = 10000;
# the S-curve figure is made at this peak strip over noise
const scurve_snr = 7.0;
# bins of the hit position within the strip for the S-curve
const nbins = 20;


### scan ###

function scan()
    results = Dict();
    println("\n| width | peak S/N | ", join([rpad(n, 22) for n in option_names], " | "), " |")
    println("|---|---|", "---|"^3)
    for w in widths, snr in peak_snrs
        A = snr*amp_noise/peak_fraction(w);
        sim = simulate(number_events=events_per_cell, amp_signal=A, cluster_width=w);
        results[(w, snr)] = sim;
        cells = String[];
        for (name, key) in zip(option_names, option_keys)
            csig, cmn, out = summary("w=$(w) S/N=$(snr) $(name)", sim[key]);
            push!(cells, rpad("σ $(csig)  out $(out) %", 22));
        end
        println("| $(w) | $(snr) | ", join(cells, " | "), " |");
    end
    return results
end

# core sigma and outlier fraction against cluster width, one column per S/N
function plot_scan(results)
    pl = plot(layout=grid(2, length(peak_snrs)), size=(500*length(peak_snrs), 800),
              plot_title="threshold $(n_threshold) sigma, $(events_per_cell) events per point, noise rms $(amp_noise) ± $(round(amp_noise/3, digits=1))",
              legend=:topright, xscale=:log2, xticks=(widths, string.(widths)), xlab="cluster sigma / strips",
              left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    for (j, snr) in enumerate(peak_snrs)
        for (i, (name, key)) in enumerate(zip(option_names, option_keys))
            csig = Float64[]; outl = Float64[];
            for w in widths
                arr = results[(w, snr)][key];
                core = arr[abs.(arr) .< 2.0];
                push!(csig, std(core)); push!(outl, 100*(1-length(core)/length(arr)));
            end
            plot!(pl[j], widths, csig, label=name, color=option_colors[i], marker=option_markers[i], markersize=5, lw=2,
                  title="peak strip / noise = $(snr)", ylab="core sigma (|res|<2) / strips", ylim=(0, 1.2))
            plot!(pl[length(peak_snrs)+j], widths, outl, label=name, color=option_colors[i], marker=option_markers[i], markersize=5, lw=2,
                  ylab="|res| > 2 strips / %")
        end
    end
    savefig(pl, "clusternoisr_scan_n$(Int(n_threshold*100)).png")
    pl
end

# S-curves: mean residual against the hit position within a strip, one panel
# per cluster width at one S/N
function plot_scurve(results)
    edges = LinRange(0, strip_size, nbins+1);
    centres = 0.5*(edges[1:end-1] .+ edges[2:end]);
    ps = plot(layout=grid(1, length(widths)), size=(500*length(widths), 450),
              plot_title="mean residual against hit position within the strip, peak strip / noise = $(scurve_snr), threshold $(n_threshold) sigma",
              legend=:topright, xlab="hit position within strip / strip size", ylab="mean residual / strips",
              left_margin=8Plots.mm, bottom_margin=5Plots.mm)
    for (j, w) in enumerate(widths)
        sim = results[(w, scurve_snr)];
        frac = mod.(sim.hits, strip_size);
        for (i, (name, key)) in enumerate(zip(option_names, option_keys))
            arr = sim[key];
            m = [mean(arr[(edges[b] .<= frac .< edges[b+1]) .& (abs.(arr) .< 2.0)]) for b in 1:nbins];
            plot!(ps[j], centres, m, label=name, color=option_colors[i], marker=option_markers[i], markersize=4, lw=2,
                  title="cluster sigma $(w) strips")
        end
        hline!(ps[j], [0.0], color=:gray, linestyle=:dash, label="")
    end
    savefig(ps, "clusternoisr_scurve_SNR$(Int(scurve_snr))_n$(Int(n_threshold*100)).png")
    ps
end


### main ###

t1 = time_ns();
results = scan();
println("scan took $((time_ns()-t1)/1.0e9) s")
plot_scan(results);
plot_scurve(results);
println("### total time $((time_ns()-t0)/1.0e9) s ###")
