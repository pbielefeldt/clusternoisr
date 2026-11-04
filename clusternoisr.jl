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
# (C) 2026, Philipp Bielefeldt, Uni Bonn                                       #
# released under the terms of Mozilla Public License (see file)                #
################################################################################
t0=time_ns();


### parameters ###

# number of events to loop over
const number_events = 2000;
# can be set to 1 without loss of generality
const strip_size = 1;
# width of the signal cluster, in units of the strip size
const cluster_width = 5;
# number of strips in the detector (size of arrays)
const number_strips = 255;
# in a.u., height of the signal
const amp_signal = 200;
# rms of the noise (after pedestal subtraction), mean over the strips
const amp_noise = 10;
# zero-suppression threshold in multiples of the strip's own noise rms:
#  - 3.0:
#    - COMPASS GEM trackers, Altunbas NIM A 490 (2002) 177, sec. 7
#    - ALICE TPC, Alme et al. NIM A 622 (2010) 316, sec. 10.3.1
#  - 3.5: Bonn GEM-TPC cosmics, Berger PhD 2015, sec. VIII.1
#  - 4.5:
#    - T2K TPCs, Abgrall et al. NIM A 637 (2011) 25, sec. 7.1.1
#    - FOPI GEM-TPC, Boehmer et al. NIM A 737 (2014) 214, tab. 1)
const n_threshold = 3.0;
# the hit position is distributed uniformly over the plane, excluding a margin
# at both edges (a fixed position makes the residual distribution a set of
# discrete values, and makes the result depend on the noise levels of the few
# strips around that position)
hit_margin = 20*strip_size;


### includes ###

t1 = time_ns();
using Plots;
using SpecialFunctions;
using Statistics;
using Random;
Random.seed!(43);
t2 = time_ns();
println("includes took $((t2-t1)/1.0e6) ms")


### functions ###

# generate a random gauss/normal distribution
function rand_norm(mu, sigma, N=1)
    mu .+ sigma .* randn(N)
end


# the signal, centred at mu, with a Gauss sigma
# it is received as the integral over the Gauss distribution - which is the
# error function from the beginning of the strip (x1) to the end of the strip
# (x2), where the first strip is from 0. to strip_size and so on.
function rand_erf(i, mu, sigma=cluster_width)
    x1 = (i-1)*strip_size;
    x2 = i*strip_size;
    erf((x2-mu)/sigma)-erf((x1-mu)/sigma)
end

# "strip plane" containing the rms noise rms information:
# the channels of one front-end differ in their noise (the spread is set 
# to 1/3 of the mean here), the rms is positive by construction
noiselevel = abs.(rand_norm(amp_noise, amp_noise/3.0, number_strips));

# fill all strips with Gaussian noise of mean zero and the strip's rms. 
# This is the baseline-subtracted signal of a channel without a hit: the 
# pedestal (the mean of the baseline) is subtracted per channel in the 
# front-end/ADC firmware, the remaining fluctuation is Gaussian and its rms 
# is the noise. Negative values occur and are kept here; the zero 
# suppression below removes them.
function set_noise!(strips_arr, nlevel_arr)
    strips_arr .= [rand_norm(0.0, ai)[1] for ai in nlevel_arr];
end

# signal on plane
function set_signal!(arr, mu, A=amp_signal)
    arr .= A.*rand_erf.(1:number_strips, mu);
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


### storage arrays ###

# the "final output" array of data
# residuals here are the difference between Monte Carlo Truth centre of gravity
# and the reconstructed CoGs
# we have three: the CoG from the pedestal-corrected amplitudes of the strips
# above threshold (residuals_0s_arr, the treatment of all readouts referenced
# above), from these amplitudes minus 1σ of the strip (residuals_1s_arr), and
# minus the full threshold nσ (residuals_nc_arr, "noise corrected")
residuals_nc_arr = zeros(number_events);
residuals_1s_arr = zeros(number_events);
residuals_0s_arr = zeros(number_events);
hit_positions    = zeros(number_events);
# keeps Gaussian noise for each pad (re-calculated per event)
enoise_arr = zeros(number_strips);
# a gaus smeared over several pads
signal_arr = zeros(number_strips);
# threshold subtracted
cutted_nc_arr = zeros(number_strips);
# 1σ subtracted
cutted_1s_arr = zeros(number_strips);
# only zero suppression
cutted_0s_arr = zeros(number_strips);


### main loop over events ###

t1 = time_ns();
for c in 1:number_events
    # randomly positioned hit (MC truth position)
    # uniform over the plane
    hit_mu = hit_margin + (number_strips*strip_size - 2*hit_margin)*rand();
    set_noise!(enoise_arr, noiselevel);
    set_signal!(signal_arr, hit_mu);
    # combine Gaussian noise and signal to measured (pedestal-corrected) data
    data_arr = enoise_arr.+signal_arr
    # zero suppression with the threshold of each strip, n times its own noise rms
    noise_cut = n_threshold.*noiselevel;
    # all entries in data that are smaller than the threshold are suppressed/set
    # to zero. The surviving strips keep their pedestal-corrected amplitude
    # (the readouts above), or are reduced by 1σ or by the threshold
    cutted_0s_arr .= [d < t ? 0.0 : d   for (d,t,s) in zip(data_arr, noise_cut, noiselevel)]; # only suppress noisy strips
    cutted_1s_arr .= [d < t ? 0.0 : d-s for (d,t,s) in zip(data_arr, noise_cut, noiselevel)]; # cut 1 sigma
    cutted_nc_arr .= [d < t ? 0.0 : d-t for (d,t,s) in zip(data_arr, noise_cut, noiselevel)]; # cut n sigma
    # caclulate the residual
    # not the pull: there is not good meaure for the measurement uncertainty,
    # and it was only a constant value anyway ...
    residual_0s = (hit_mu - get_cog(cutted_0s_arr));
    residual_1s = (hit_mu - get_cog(cutted_1s_arr));
    residual_nc = (hit_mu - get_cog(cutted_nc_arr));
    # for every event, write out the residual to histo
    residuals_0s_arr[c] = residual_0s;
    residuals_1s_arr[c] = residual_1s;
    residuals_nc_arr[c] = residual_nc;
    hit_positions[c]    = hit_mu;

    # some timing
    tX = (time_ns()-t1)/1.0e9;
    if 10*c % number_events == 0
        println("$c events with $(round(c/tX, digits=2)) Hz")
    end
end
t2 = time_ns();
println("main loop took $((t2-t1)/1.0e9) s")
t1 = time_ns();

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
    return csig, cmn
end

function make_plot(;xlim=2.0)
    # signal of the strip under the hit centre (the peak strip) in units of the
    # mean noise: this is the quantity that decides whether the threshold cuts
    # into the cluster core (amp_signal/amp_noise is the *total* signal over noise)
    peak = amp_signal*rand_erf(1, 0.5*strip_size); # hit at the centre of strip 1
    println("total signal / noise: ", amp_signal/amp_noise, "   peak strip / noise: ", round(peak/amp_noise, digits=1), "   threshold: ", n_threshold, " sigma")
    pl = plot(layout=grid(1,3), size=(1500,500), legend=false)
    names = ["only zeros suppressed", "1 sigma subtracted", "threshold ($(n_threshold) sigma) subtracted"];
    arrs = [residuals_0s_arr, residuals_1s_arr, residuals_nc_arr];
    for (i, (name, arr)) in enumerate(zip(names, arrs))
        histogram!(pl[i], arr, bins=LinRange(-xlim,xlim,100), xlab="residual ($(name)) / #strips");
        hmax = max(pl[i][1][:y][1:6:end]...);
        sig, mn = summary(name, arr);
        annotate!([(xlim,0.8*hmax,text("sigma: $(sig) \nmean: $(mn)",:right))], subplot=i)
    end
    savefig(pl, "clusternoisr_SNR$(Int(round(amp_signal/amp_noise)))_n$(n_threshold).png")
    pl
end

t2 = time_ns();
println("### total time $((t2-t0)/1.0e9) s ###")

make_plot()
