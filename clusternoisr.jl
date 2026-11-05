#!/usr/bin/env julia
################################################################################
# A tool to simulate the behaviour of a strip plane (e.g. in a gaseous         #
# detector), where noise and signal are treated in different ways: Strips      #
# below the noise threshold (a multiple 3-5 of the noise σ of the strip) are   #
# set to a value of zero (zero suppression).                                   #
# Then, the hit position (centre of gravity) is calculated from the surviving  #
# strips in three ways: with their entire (pedestal-corrected) amplitude, with #
# 1σ subtracted, or with the full threshold nσ subtracted. I wrote this to     #
# investigate the effects either has on the resolution (here, the residual     #
# distribution width) and the bias of the detector.                            #
#                                                                              #
# The noise and threshold model follows what GEM/TPC readouts do (references   #
# at the functions): the pedestal (baseline) of each channel is subtracted,    #
# what remains is Gaussian noise with mean zero and a channel-specific rms,    #
# and the zero suppression keeps samples above n times that rms with their     #
# pedestal-corrected amplitude.                                                #
#                                                                              #
# usage:  julia clusternoisr.jl                                                #
# runs one simulation with the parameters below and plots the residual        #
# histograms of the three options. The parameters can be overridden from the  #
# environment, e.g. AMP_SIGNAL=1200 AMP_NOISE=20 NUMBER_EVENTS=4000 julia ...  #
# scan.jl includes this file and runs the simulation over a parameter grid.    #
#                                                                              #
# (C) 2026, Philipp Bielefeldt, Uni Bonn                                       #
# released under the terms of Mozilla Public License (see file)                #
################################################################################
t0=time_ns();


### parameters ###

# each can be overridden by the environment variable of the same name in
# capitals, e.g. AMP_SIGNAL=1200
param(name, default) = parse(typeof(default), get(ENV, uppercase(name), string(default)));

# number of events to loop over
const number_events = param("number_events", 1000);
# can be set to 1 without loss of generality
const strip_size = 1;
# width (Gauss sigma) of the signal cluster, in units of the strip size
const cluster_width = param("cluster_width", 3.5);
# number of strips in the detector (size of arrays)
const number_strips = 255;
# in a.u., total charge of the hit (the integral of its Gauss over all strips);
# the strip under the hit centre gets peak_fraction(cluster_width) of it
const amp_signal = param("amp_signal", 600.0);
# rms of the noise (after pedestal subtraction), mean over the strips
const amp_noise = param("amp_noise", 10.0);
# zero-suppression threshold in multiples of the strip's own noise rms:
#  - 3.0:
#    - COMPASS GEM trackers, Altunbas NIM A 490 (2002) 177, sec. 7
#    - ALICE TPC, Alme et al. NIM A 622 (2010) 316, sec. 10.3.1
#  - 3.5: Bonn GEM-TPC cosmics, Berger PhD 2015, sec. VIII.1
#  - 4.5:
#    - T2K TPCs, Abgrall et al. NIM A 637 (2011) 25, sec. 7.1.1
#    - FOPI GEM-TPC, Boehmer et al. NIM A 737 (2014) 214, tab. 1)
const n_threshold = param("n_threshold", 3.0);
# the hit position is distributed uniformly over the plane, excluding a margin
# at both edges (a fixed position makes the residual distribution a set of
# discrete values, and makes the result depend on the noise levels of the few
# strips around that position)
const hit_margin = 20*strip_size;
# seed of the random number generator (noise rms of the strips, hits, noise)
const seed = param("seed", 43);


### includes ###

t1 = time_ns();
using Plots;
using SpecialFunctions;
using Statistics;
using Random;
t2 = time_ns();
println("includes took $((t2-t1)/1.0e6) ms")


### functions ###

# generate a random gauss/normal distribution
function rand_norm(mu, sigma, N=1)
    mu .+ sigma .* randn(N)
end


# the fraction of the signal on strip i: the signal is a Gauss centred at mu
# with standard deviation sigma and unit integral, and the strip receives its
# integral from the beginning of the strip (x1) to the end of the strip (x2),
# where the first strip is from 0. to strip_size and so on. The integral of
# the normal distribution is Φ(x) = (1 + erf(x/√2))/2.
function rand_erf(i, mu, sigma)
    x1 = (i-1)*strip_size;
    x2 = i*strip_size;
    0.5*(erf((x2-mu)/(sqrt(2)*sigma))-erf((x1-mu)/(sqrt(2)*sigma)))
end

# fraction of amp_signal on the strip under the hit centre (hit at the centre
# of strip 1): the peak strip amplitude is peak_fraction(sigma)*amp_signal.
# This, not amp_signal/amp_noise, decides whether the threshold cuts into
# the cluster core.
peak_fraction(sigma) = rand_erf(1, 0.5*strip_size, sigma);

# fill all strips with Gaussian noise of mean zero and the strip's rms.
# This is the baseline-subtracted signal of a channel without a hit: the
# pedestal (the mean of the baseline) is subtracted per channel in the
# front-end/ADC firmware, the remaining fluctuation is Gaussian and its rms
# is the noise. Negative values occur and are kept here; the zero
# suppression below removes them.
function set_noise!(strips_arr, rms_arr)
    strips_arr .= [rand_norm(0.0, ai)[1] for ai in rms_arr];
end

# signal on plane
function set_signal!(arr, mu, A, sigma)
    arr .= A.*rand_erf.(1:number_strips, mu, sigma);
end

# truncate an array
# this works mostly like a cluster finder: the strip with the largest amplitude
# seeds the cluster, which is the group of adjacent strips above threshold
# (COMPASS, Altunbas et al. 2002, sec. 7: "a cluster is defined as the group of
# adjacent strips with a content exceeding threshold", seeded at the channel
# with the highest content). Works up and down from the seed until it hits an
# entry equal zero (<1e-18), returns an array of only these strips.
# this is the input for the centre of gravity calculation.
function trunc_arr(arr)
    N = argmax(arr); # number of bin with the largest amplitude
    z = 1e-16;  # used as a replacement for z (set to machine precision?)
    ret = zeros(length(arr)); # returned array
    # loop to the left
    for i in N:-1:1
        if arr[i] < z
            break;
        end
        ret[i] = arr[i]; # copy values until reaching zero
    end
    # loop to the right
    #NOTE: will probably use arr[N] twice in the loop (but that's ok)
    for i in N:+1:length(arr)
        if arr[i] < z
            break;
        end
        ret[i] = arr[i]; # copy values until reaching zero
    end
    return ret
end

# calculate the centre of gravity
# the CoG is th mean of all bin-heights * bin-positions
# first, it has some cluster finding: starts at the maximum bin, goes left and
# right until it finds a bin with value 0.0 and stops there
function get_cog(arr)
    tarr = trunc_arr(arr);
    if sum(tarr)==0
        return 0.0
    end
    sum((0.5:(length(tarr)-0.5)) .* tarr)/sum(tarr)
end


### the simulation ###

# one simulation: a strip plane with random noise rms per strip, number_events
# hits. Returns the residuals (Monte Carlo truth hit position minus the
# reconstructed CoG) of the three options, and the hit positions:
#  res_0s: CoG from the pedestal-corrected amplitudes of the strips above
#          threshold (the treatment of all readouts referenced above)
#  res_1s: from these amplitudes minus 1σ of the strip
#  res_nc: minus the full threshold nσ ("noise corrected")
function simulate(; number_events=number_events, amp_signal=amp_signal,
                    amp_noise=amp_noise, cluster_width=cluster_width,
                    n_threshold=n_threshold, seed=seed)
    Random.seed!(seed);
    # "strip plane" containing the noise rms of each strip (sigma_strip): the
    # rms of the pedestal-subtracted fluctuations of a channel without signal.
    # the channels of one front-end differ in their noise (the spread is set
    # to 1/3 of the mean here), the rms is positive by construction
    noise_rms = abs.(rand_norm(amp_noise, amp_noise/3.0, number_strips));
    # zero-suppression threshold of each strip, n times its noise rms. The
    # pedestal is assumed to be subtracted already, so it is not added here.
    threshold = n_threshold.*noise_rms;

    res_0s = zeros(number_events);
    res_1s = zeros(number_events);
    res_nc = zeros(number_events);
    hits   = zeros(number_events);
    # keeps Gaussian noise for each pad (re-calculated per event)
    enoise_arr = zeros(number_strips);
    # a gaus smeared over several pads
    signal_arr = zeros(number_strips);
    # only zero suppression / 1σ subtracted / threshold subtracted
    cutted_0s_arr = zeros(number_strips);
    cutted_1s_arr = zeros(number_strips);
    cutted_nc_arr = zeros(number_strips);

    for c in 1:number_events
        # randomly positioned hit (MC truth position), uniform over the plane
        hit_mu = hit_margin + (number_strips*strip_size - 2*hit_margin)*rand();
        set_noise!(enoise_arr, noise_rms);
        set_signal!(signal_arr, hit_mu, amp_signal, cluster_width);
        # combine Gaussian noise and signal to measured (pedestal-corrected) data
        data_arr = enoise_arr.+signal_arr
        # zero suppression with the threshold of each strip, n times its own noise rms:
        # all entries in data that are smaller than the threshold are suppressed/set
        # to zero. The surviving strips keep their pedestal-corrected amplitude
        # (the readouts above), or are reduced by 1σ or by the threshold
        cutted_0s_arr .= [d < t ? 0.0 : d   for (d,t,s) in zip(data_arr, threshold, noise_rms)]; # only suppress noisy strips
        cutted_1s_arr .= [d < t ? 0.0 : d-s for (d,t,s) in zip(data_arr, threshold, noise_rms)]; # cut 1 sigma
        cutted_nc_arr .= [d < t ? 0.0 : d-t for (d,t,s) in zip(data_arr, threshold, noise_rms)]; # cut n sigma
        # caclulate the residual
        # not the pull: there is not good meaure for the measurement uncertainty,
        # and it was only a constant value anyway ...
        res_0s[c] = hit_mu - get_cog(cutted_0s_arr);
        res_1s[c] = hit_mu - get_cog(cutted_1s_arr);
        res_nc[c] = hit_mu - get_cog(cutted_nc_arr);
        hits[c]   = hit_mu;
    end
    return (res_0s=res_0s, res_1s=res_1s, res_nc=res_nc, hits=hits, noise_rms=noise_rms)
end


### output ###

# the three options, in fixed order, with a fixed colour and marker each
const option_names   = ["only zeros suppressed", "1 sigma subtracted", "threshold subtracted"];
const option_keys    = [:res_0s, :res_1s, :res_nc];
const option_colors  = ["#2a78d6", "#eb6834", "#1baf7a"];
const option_markers = [:circle, :utriangle, :square];

# width and mean of the residuals; clusters that were truncated at a strip
# below threshold or lost altogether give residuals of several strips, so the
# fraction of residuals outside +/-2 strips and the rms of the rest ("core") are
# given as well
function summary(name, arr)
    sig, mn = round(sqrt(var(arr)),digits=3), round(mean(arr),digits=3);
    core = arr[abs.(arr) .< 2.0];
    out = round(100*(1-length(core)/length(arr)), digits=1);
    csig, cmn = round(sqrt(var(core)),digits=3), round(mean(core),digits=3);
    println(rpad(name, 32), " sigma: ", sig, "  mean: ", mn, "  |  core (|res|<2): sigma ", csig, "  mean ", cmn, "  outliers ", out, " %");
    return csig, cmn, out
end

# residual histograms of the three options for one simulation
function make_plot(sim; xlim=5.0, amp_signal=amp_signal, amp_noise=amp_noise,
                   cluster_width=cluster_width, n_threshold=n_threshold)
    # signal of the strip under the hit centre (the peak strip) in units of the
    # mean noise: this is the quantity that decides whether the threshold cuts
    # into the cluster core (amp_signal/amp_noise is the *total* signal over noise)
    peak = amp_signal*peak_fraction(cluster_width);
    println("total signal / noise: ", amp_signal/amp_noise, "   peak strip / noise: ", round(peak/amp_noise, digits=1), "   threshold: ", n_threshold, " sigma")
    pl = plot(layout=grid(1,3), size=(1500,500), legend=false,
              plot_title="peak strip / noise $(round(peak/amp_noise, digits=1)) (total charge / noise $(round(amp_signal/amp_noise, digits=1))), threshold $(n_threshold) sigma, cluster sigma $(cluster_width) strips, $(length(sim.hits)) events")
    for (i, (name, key)) in enumerate(zip(option_names, option_keys))
        arr = sim[key];
        histogram!(pl[i], arr, bins=LinRange(-xlim,xlim,100), title=name, xlab="residual / #strips", color=option_colors[i], linecolor=option_colors[i]);
        hmax = max(pl[i][1][:y][1:6:end]...);
        sig, mn, out = summary(name, arr);
        annotate!([(xlim,0.8*hmax,text("sigma: $(sig) \nmean: $(mn)",:right))], subplot=i)
    end
    savefig(pl, "clusternoisr_peakSNR$(Int(round(peak/amp_noise)))_n$(Int(n_threshold*100)).png")
    pl
end


### main (not run when included from scan.jl) ###

if abspath(PROGRAM_FILE) == @__FILE__
    t1 = time_ns();
    sim = simulate();
    t2 = time_ns();
    println("main loop took $((t2-t1)/1.0e9) s")
    make_plot(sim);
    println("### total time $((time_ns()-t0)/1.0e9) s ###")
end
