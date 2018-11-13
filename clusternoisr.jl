#!/usr/bin/env julia

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
const amp_noise = 10;

# the MC truth of hit width
hit_sigma = strip_size*0.667;


### includes ###
t1 = time_ns();
using Plots;
using SpecialFunctions;
using Statistics;
using Random;
#using Compat, Distributions;

Random.seed!(43);
t2 = time_ns();
println("includes took $((t2-t1)/1.0e6) ms")

### functions ###

# generate a random gauss/normal distribution
function rand_norm(mu, sigma, N=1)
    mu .+ sigma .* randn(N)
end

# get a random distribution following exponential decay for every strip (pink
# noise)
rand_exp(N=1) = -1.0*log.(rand(N));

# the signal, centred at mu, with a Gauss sigma
# it is received as the integral over the Gauss distribution - which is the
# error function from the beginning of the strip (x1) to the end of the strip
# (x2), where the first strip is from 0. to strip_size and so on.
function rand_erf(i, mu, sigma=cluster_width)
    x1 = (i-1)*strip_size;
    x2 = i*strip_size;
    erf((x2-mu)/sigma)-erf((x1-mu)/sigma)
end

# "strip plane" containing the mean noise information
noiselevel = rand_norm(amp_noise, amp_noise/3.0, number_strips);

# fill all strips with exponential noise around their noise levels
function set_noise!(strips_arr, nlevel_arr)
    # note: the abs() is a bit of an unphysical hack ...
    strips_arr .= [rand_exp(1)[1]*abs(ai) for ai in nlevel_arr];
end

# signal on plane
function set_signal!(arr, mu, A=amp_signal)
    arr .= A.*rand_erf.(1:number_strips, mu);
end

# truncate an array
# this works mostly like a cluster finder: it is assumed the central strip is
# part of the cluster (!!!11!). works up and down from it until it hits an
# entry equal zero (<1e-18), returns an array of only the central strips.
# this is the input for the centre of gravity calculation.
function trunc_arr(arr)
    N = Int(trunc(length(arr)/2)); # number of central bin
    z = 1e-16;  # used as a replacement for z (set to machine precision?)
    ret = zeros(2*N); # returned array

    # loop to the left
    for i in N:-1:1
        if arr[i] < z
            break;
        end
        ret[i] = arr[i]; # copy values until reaching zero
    end

    # loop to the right
    #NOTE: will probably use arr[N] twice in the loop (but that's ok)
    for i in N:+1:2*N
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
    sum((0.5:(length(tarr)-0.5)) .* tarr)/sum(tarr)
end
# function get_cog(arr)
#
#     # contains all values of arr that are neighbours to the cluster
#     # v[i][1] is the position (i.e. strip number)
#     # v[i][2] is the weight (i.e. arr entry)
#     v = [];
#     max_bin_n = findmax(arr)[2];
#
#     i::Int = max_bin_n; # counter
#     max_i::Int = length(arr);
#
#     # left of maximum
#     while i > 0
#         # fill data to v (until minimum is reached)
#         if arr[i] > 0.0
#             push!(v, (i, arr[i]));
#             i = i-1;
#         else
#             break
#         end
#     end
#
#     # reset i (don't count arr[i] twice, though)
#     i = max_bin_n+1;
#
#     # right of maximum
#     while i < max_i
#         # fill data to v (until minimum is reached)
#         if arr[i] > 0.0
#             push!(v, (i, arr[i]));
#             i = i+1;
#         else
#             break
#         end
#     end
#     # cog in units of strip_size
#     # sum((0.5:(length(v)-0.5)) .* v)/sum(v) #TODO: Only works if v is sorted by i
#     cog = 0.0;
#     sumx = 0.0;
#     for x in v
#         cog = cog + ((x[1]-0.5)*x[2]);
#         sumx = sumx + x[2];
#     end
#     cog/sumx
# end


### storage arrays ###
# the "final output" array of data
# residuals here are the difference between Monte Carlo Truth centre of gravity
# and the reconstructed CoGs
# we have two: one under the assumption that the CoG is calculated with cutted
# data (residuals_nc_arr), one where all signal (above the pedestal, which we
# do not account for) is used (residuals_0s_arr)
residuals_nc_arr = zeros(number_events);
residuals_0s_arr = zeros(number_events);

# keeps exponential noise for each pad (re-calculated per event)
enoise_arr = zeros(number_strips);
# a gaus smeared over several pads
signal_arr = zeros(number_strips);
# noise cut
cutted_nc_arr = zeros(number_strips);
# only zero suppression
cutted_0s_arr = zeros(number_strips);


### main loop over events ###
t1 = time_ns();
for c in 1:number_events
    # randomly positioned hit (MC truth position)
    # hit_mu = rand_norm((number_strips+1)*strip_size/2.0, hit_sigma);
    hit_mu = (number_strips+1)*strip_size/2.0;

    set_noise!(enoise_arr, noiselevel);
    set_signal!(signal_arr, hit_mu);

    # combine exponential noise and signal to measured data
    data_arr = enoise_arr.+signal_arr

    # bar(data_arr)
    # it is usual to have 3×noise as cutoff
    noise_cut = 3.0*amp_noise;

    # data array without the noise
    # all entries in data that are smaller than the noise cut are suppressed/set
    # to zero, those above the cutoff amplitude are reduced by cutoff
    cutted_nc_arr = [d < noise_cut ? 0 : d-noise_cut for d in data_arr]; # cut 3 sigma
    cutted_0s_arr = [d < noise_cut ? 0 : d for d in data_arr]; # only suppress noisy strips

    # caclulate the residual
    # not the pull: there is not good meaure for the measurement uncertainty,
    # and it was only a constant value anyway ...
    #println("cog ", get_cog(cutted_nc_arr))
    residual_nc = (hit_mu[1] - get_cog(cutted_nc_arr));
    residual_0s = (hit_mu[1] - get_cog(cutted_0s_arr));

    # for every event, write out the pull to histo
    # push!(residuals_nc_arr, residual_nc); #TODO: This seems to have many 0's -- why?
    # push!(residuals_0s_arr, residual_0s);
    residuals_nc_arr[c] = residual_nc;
    residuals_0s_arr[c] = residual_0s;

    # some timing
    tX = (time_ns()-t1)/1.0e9;
    if 10*c % number_events == 0
        println("$c events with $(round(c/tX, digits=2)) Hz")
    end
end
t2 = time_ns();
println("main loop took $((t2-t1)/1.0e9) s")

#t1 = time_ns();
function make_plot(;xlim=5.0)
    pl = plot(layout=grid(1,2), size=(1000,500), legend=false)
    h1 = histogram!(pl[1], residuals_nc_arr, bins=LinRange(-xlim,xlim,100), xlab="residual (noise corrected, 3 sigma) / #strips");
    h2 = histogram!(pl[2], residuals_0s_arr, bins=LinRange(-xlim,xlim,100), xlab="residual (only zeros suppressed) / #strips");
    h1max = max(pl[1][1][:y][1:6:end]...);
    h2max = max(pl[2][1][:y][1:6:end]...);
    # @show h1max, h2max
    sig1, sig2 = round(sqrt(var(residuals_nc_arr)),digits=3), round(sqrt(var(residuals_0s_arr)),digits=3);
    mean1, mean2 = round(mean(residuals_nc_arr),digits=3), round(mean(residuals_0s_arr),digits=3);
    annotate!([(xlim,0.8*h1max,text("sigma: $(sig1) \nmean: $(mean1)",:right))], subplot=1)
    annotate!([(xlim,0.8*h2max,text("sigma: $(sig2) \nmean: $(mean2)",:right))], subplot=2)
end

# plot(
#     layout=grid(2,3), legend=false,
#     bar([0:255], noiselevel, xlabel="noise level"),
#     bar([0:255], enoise_arr, xlabel="exp. noise"),
#     bar([0:255], signal_arr, xlabel="signal"),
#     bar([0:255], data_arr, xlabel="signal+noise"),
#     bar([0:255], cutted_0s_arr, xlabel="0 supressed"),
#     bar([0:255], cutted_nc_arr, xlabel="amp cutted")
#     )

t2 = time_ns();
println("### total time $((t2-t0)/1.0e9) s ###")

make_plot()
