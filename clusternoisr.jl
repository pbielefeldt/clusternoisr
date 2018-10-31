#!/usr/bin/env julia

### parameters ###

# number of events to loop over
const number_events = 2;

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
using Statistics

### functions ###

# get a random distribution following exponential decay for every strip (pink
# noise)
rand_exp(N=1) = -amp_noise*log.(rand(N))

# the signal, centred at mu, with a Gauss sigma
# it is received as the integral over the Gauss distribution - which is the
# error function from the beginning of the strip (x1) to the end of the strip
# (x2), where the first strip is from 0. to strip_size and so on.
function rand_erf(i, mu, sigma=cluster_width)
    x1 = (i-1)*strip_size;
    x2 = i*strip_size;
    erf((x2-mu)/sigma)-erf((x1-mu)/sigma)
end

function set_noise!(arr)
    # noisy plane
    arr .= rand_exp(number_strips);
end

function set_signal!(arr, mu, A=amp_signal)
    # signal on plane
    arr .= A.*rand_erf.(1:number_strips, mu);
end


enoise_arr = fill(0.0,number_strips);
signal_arr = fill(0.0,number_strips);
#cutted_arr = fill(0.0,number_strips);

residuals_arr = [];

### main loop over events ###
for c in 1:number_events
    # randomly positioned hit
    mu = number_strips*strip_size/2;

    set_noise!(enoise_arr);
    set_signal!(signal_arr, mu);

    # combine exponential noise and signal to measured data
    data_arr = enoise_arr.+signal_arr;

    noise_cut = 3.0*amp_noise;
    cutted_arr = [d < noise_cut ? 0 : d-noise_cut for d in data_arr]

    push!(residuals_arr, (mean(cutted_arr),mean(data_arr)))
end

residuals_arr
