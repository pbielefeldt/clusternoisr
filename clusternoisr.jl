#!/usr/bin/env julia

### parameters ###

# number of events to loop over
const number_events = 2000;

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
noiselevel = rand_norm(amp_noise, amp_noise, number_strips);

# fill all strips with exponential noise around their noise levels
function set_noise!(strips_arr, nlevel_arr)
    # noisy plane
    strips_arr = [rand_exp(1)*amp_noise*ai for ai in nlevel_arr];
end

function set_signal!(arr, mu, A=amp_signal)
    # signal on plane
    arr .= A.*rand_erf.(1:number_strips, mu);
end

# calculate the centre of gravity
# the CoG is th mean of all bin-heights * bin-positions
# caution! expects an array from 1 to number_strips
function get_cog(arr)
    ret = [];
    tot_amp = 0;
    for x in 1:length(arr)
        amp = arr[x];
        push!(ret, (x*strip_size)*amp);
        tot_amp += amp;
    end
    mean(ret)/tot_amp
end

# using an array as a histogram
function hfill!(arr, x, amp=1)
    bin::Int = 1+(x÷(strip_size*number_strips)) ;
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
const xstart = -16;
const xend = 48;
residuals_nc_arr = zeros(1+xend-xstart);
residuals_0s_arr = zeros(1+xend-xstart);

xarr = [xstart:xend]; # I still hope there is a sane way for histos in Julia …?

### main loop over events ###
for c in 1:number_events
    # randomly positioned hit (MC truth position)
    hit_mu = rand_norm((number_strips+1)*strip_size/2.0, hit_sigma);

    set_noise!(enoise_arr, noiselevel);
    set_signal!(signal_arr, hit_mu);

    # combine exponential noise and signal to measured data
    data_arr = enoise_arr.+signal_arr

    # it is usual to have 3×noise as cutoff
    noise_cut = 3.0*amp_noise;

    # data array without the noise
    # all entries in data that are smaller than the noise cut are suppressed/set
    # to zero, those above the cutoff amplitude are reduced by cutoff
    cutted_nc_arr = [d < noise_cut ? 0 : d-noise_cut for d in data_arr]; # cut 3 sigma
    cutted_0s_arr = [d < noise_cut ? 0 : d for d in data_arr]; # only suppress noisy strips

    # calculate centre of gravity for arr_0 and arr_sigma, write pulls
    pull_nc = (hit_mu[1] - get_cog(cutted_nc_arr))/hit_sigma; #TODO: mu is seen as array of length 1 >.<
    pull_0s = (hit_mu[1] - get_cog(cutted_0s_arr))/hit_sigma;

    # for every event, write out the pull to histo
    #push!(residuals_nc_arr, pull_nc);
    #push!(residuals_0s_arr, pull_0s);
    hfill!(residuals_nc_arr, pull_nc);
    hfill!(residuals_0s_arr, pull_0s);
end

pull_nc_hist = bar(residuals_nc_arr, xarr, xlabel="pull (noise cut applied)");
pull_0s_hist = bar(residuals_0s_arr, xarr, xlabel="pull (only zero suppression)");

# draw
plot(pull_nc_hist, pull_0s_hist, layout=(1,2), legend=false)
