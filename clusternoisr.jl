#!/usr/bin/env julia

### parameters ###

# number of events to loop over
const number_events = 10000;

# can be set to 1 without loss of generality
const strip_size = 1;

# width of the signal cluster, in units of the strip size
const cluster_width = 3;

# number of strips in the detector (size of arrays)
const number_strips = 30;

# in a.u., height of the signal
const amp_signal = 50;
const amp_noise = 1;


### includes ###
using Plots;
using SpecialFunctions;
using Statistics;
#using Compat, Random, Distributions;

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
    strips_arr .= [rand_exp(1)[1]*amp_noise*abs(ai) for ai in nlevel_arr];
end

function set_signal!(arr, mu, A=amp_signal)
    # signal on plane
    arr .= A.*rand_erf.(1:number_strips, mu);
end

# calculate the centre of gravity
# the CoG is th mean of all bin-heights * bin-positions
# caution! expects an array from 1 to number_strips
#TODO: it would be better to have a cluster finder here!
function get_cog(arr)
    sum((0.5:(length(arr)-0.5)) .* arr)/sum(arr)
end

# using an array as a histogram
function hfill!(arr, x, xstart, xend, nbins, amp=1)
    len = (xend-xstart);
    pos = x/len;
    bin::Int = floor(Int, 1.5+pos*nbins) ;
    #println("filling bin ", bin);
    arr[bin] += amp;
end

enoise_arr = zeros(number_strips);
signal_arr = zeros(number_strips);
cutted_nc_arr = zeros(number_strips); # noise cut
cutted_0s_arr = zeros(number_strips); # only zero suppression

# the MC truth of hit width
hit_sigma = strip_size*2.0;

# the "final output" array of data
# residuals here are the difference between Monte Carlo Truth centre of gravity
# and the reconstructed CoGs
# we have two: one under the assumption that the CoG is calculated with cutted
# data (residuals_nc_arr), one where all signal (above the pedestal, which we
# do not account for) is used (residuals_0s_arr)
const xstart = -20;
const xend = 30;
const nbins = 500;
residuals_nc_arr = zeros(number_events);
residuals_0s_arr = zeros(number_events);

xarr = collect(xstart : (xend-xstart)/nbins : xend-((xend-xstart)/nbins)); # x-axis for the bar chart

### main loop over events ###
for c in 1:number_events
    # randomly positioned hit (MC truth position)
    hit_mu = rand_norm((number_strips+1)*strip_size/2.0, hit_sigma);
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
    residual_nc = (hit_mu[1] - get_cog(cutted_nc_arr));
    residual_0s = (hit_mu[1] - get_cog(cutted_0s_arr));

    # for every event, write out the pull to histo
    # push!(residuals_nc_arr, residual_nc); #TODO: This seems to have many 0's -- why?
    # push!(residuals_0s_arr, residual_0s);
    residuals_nc_arr[c] = residual_nc;
    residuals_0s_arr[c] = residual_0s;
end

function make_plot(;xlim=2.0)
    pl = plot(layout=grid(1,2), size=(1000,500), legend=false)
    h1 = histogram!(pl[1], residuals_nc_arr, bins=LinRange(-xlim,xlim,200), xlab="pull (noise corrected, 3 sigma)");
    h2 = histogram!(pl[2], residuals_0s_arr, bins=LinRange(-xlim,xlim,200), xlab="pull (only zeros suppressed)");
    h1max = max(pl[1][1][:y][1:6:end]...);
    h2max = max(pl[2][1][:y][1:6:end]...);
    # @show h1max, h2max
    sig1, sig2 = round(sqrt(var(residuals_nc_arr)),digits=3), round(sqrt(var(residuals_0s_arr)),digits=3);
    mean1, mean2 = round(mean(residuals_nc_arr),digits=3), round(mean(residuals_0s_arr),digits=3);
    annotate!([(xlim,0.8*h1max,text("sigma: $(sig1) \nmean: $(mean1)",:right))], subplot=1)
    annotate!([(xlim,0.8*h2max,text("sigma: $(sig2) \nmean: $(mean2)",:right))], subplot=2)
end

make_plot()
