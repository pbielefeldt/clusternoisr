#!/usr/bin/env julia

### parameters ###

# can be set to 1 without loss of generality
const strip_size = 1;

# width of the signal cluster, in units of the strip size
const cluster_width = 3;

# number of strips in the detector (size of arrays)
const number_strips = 30;

# in a.u., height of the signal
const amp_signal = 50;


### includes ###
using Plots;
using SpecialFunctions;


### functions ###

# get a random distribution following exponential decay for every strip (pink
# noise)
rand_exp(N=1) = -log.(rand(N))

# the signal, centred at mu, with a Gauss sigma
# it is received as the integral over the Gauss distribution - which is the
# error function from the beginning of the strip (x1) to the end of the strip
# (x2), where the first strip is from 0. to strip_size and so on.
function rand_erf(i,mu=number_strips*strip_size/2,sigma=cluster_width)
    x1 = (i-1)*strip_size;
    x2 = i*strip_size;
    erf((x2-mu)/sigma)-erf((x1-mu)/sigma)
end

enoise_arr = rand_exp(30);
signal_arr = amp_signal*rand_erf.(1:30);
data_arr = enoise_arr.+signal_arr;
bar(data_arr)
