clusternoisr
============

A small script in Julia to investigate the impact of noise on cluster size and -resolution.
Since it is unclear (at least to me …), how "noise" has to be cut off correctly, I started this toy Monte Carlo simulation.


Problem
-------

It is disputed how, in a particle detector, noise should be treated in cluster reconstruction.
The terms used here, per electronic channel (pad/strip/…) `i`:

* __pedestal__ `p_i`: the baseline of the channel, the mean reading without signal. It is measured per channel and subtracted in the front-end electronics; all amplitudes below are pedestal-subtracted, so the baseline is at zero.
* __noise__ `e_i`: the fluctuation of one reading around the baseline. It is Gaussian with mean zero and can be positive or negative; one measurement reads `d_i = s_i + e_i` for a signal `s_i`.
* __noise rms__ `σ_i`: the rms of `e_i`, measured from many readings without signal. It differs from channel to channel.
* __threshold__ `t_i = n σ_i`: the zero-suppression cut (`p_i + n σ_i` in raw ADC counts). `n` is typically between 3 and 5.

It is usual to remove all channels from analysis that have a reading (i.e. physics signal plus noise) below their threshold.

Now, for the remaining channels, it is unclear whether or not the nσ cut applied to the signal should also be applied to the analysed data. 
If done so, the remaining signal is smoother, else, the edges will be sharper.
A "smoother" (gaussian) distribution of the signal over neighboring strips/pads might make it easier to fit to them for a precise definition of the centre of the hit.
However, in reality, gauss-fits to determine the centre of a hit are rarely done.
Mostly, an averaging algorithm (Steiner's theorem or similar) is usually applied to calculate a centre of gravity of the hit.
For such calculations, the "sharper" shape of the cut ditribution might lead to smaller uncertainties.

In brief: Fitting to a broad, smooth signal can yield the best results, while for averaging a smaller number of strips/pads might be more beneficial.

Note: It might seem counter-intuitive, but yes, a broaders signal (over more than one strip), generally yields in a higher spatial resolution -- up to a certain point at least.

Since the nσ cut is somewhat arbritrary, many argue it is more accurate to only cut 1σ (the assumed noise) from the signal. 
It can also be argued that the entire (pedestal-corrected) signal should be taken for finding the centre of gravity; that is what the readouts referenced above do.

In all cases, only the strips with `d_i ≥ t_i` are used; the options differ in the amplitude these strips enter the position calculation (fit or centre of gravity) with:
1. The **signal over pedestal** `d_i`: nothing is subtracted from the surviving strips. This is what the readouts referenced below store.
2. The **signal minus 1σ** `d_i - σ_i`.
3. The **signal over threshold** `d_i - t_i = d_i - n σ_i`.


Input
-----

A 2D __strip plane__ is assumed as geometry (default: 255 strips wide).
Each strip is associated with a randomly chosen __noise rms__ `σ_i` (`noise_rms`; per simulation, not per event), that is gaussianly shaped around `amp_noise`, and the threshold `t_i = n_threshold σ_i` (`threshold`).
Take `number_events` runs in the loop. 
In each, do the following:
 * For each strip, calculate an event-specific noise. That is *Gaussian* with mean zero and the noise rms of this strip, i.e. the pedestal-corrected reading of a channel without signal. 
 * Place a hit at a position drawn uniformly over the strip plane (excluding a margin at the edges). It has a *gauss shape* with standard deviation `cluster_width` (in strips) and total charge `amp_signal` (the integral over all strips); each strip receives the integral of the Gauss over its width. The strip under the hit centre gets the fraction `erf(1/(2√2 cluster_width))` of the total charge, e.g. 11 % for `cluster_width = 3.5`, 38 % for 1 strip, 95 % for 0.25 strips. This *peak strip* amplitude over the noise rms, not `amp_signal/amp_noise`, decides whether the threshold cuts into the cluster.
 * Add up signal and Gaussian noise (`data_arr = enoise_arr.+signal_arr`)


Noise Suppression
-----------------

All strips with an amplitude below their threshold `t_i = n_threshold σ_i` are set to zero (zero suppression). The surviving strips are then treated in three ways:

__Option 1__ ("only zeros suppressed"): keep the pedestal-corrected amplitude. This is what the referenced readouts store.

__Option 2__ ("1σ subtracted"): subtract the noise rms `σ_i` of the strip.

__Option 3__ ("noise corrected"): subtract the full threshold `t_i`.


Cluster Finding
---------------

Strips far away from the hit can pass the threshold by a noise fluctuation alone: the probability per strip is `1 - Φ(n)`, about 0.13 % for `n = 3`, i.e. roughly one noise strip every three events on a 255-strip plane.
Such things of course always happen in real data.
In this very comprehensive Monte Carlo study, a simplistic approach is used to handle the situation:
The strip with the largest amplitude seeds the cluster, and the cluster is the group of adjacent strips above threshold.
From the seed, the function `trunc_arr(arr)` iterates to the left and right, until it hits a bin with amplitude zero. 
All bins beyont this bin are discarded.


Relevant Read
-------------

- COMPASS GEM trackers: Altunbas NIM A 490 (2002) 177, sec. 7
- ALICE TPC: Alme et al. NIM A 622 (2010) 316, sec. 10.3.1
- FOPI GEM-TPC: Boehmer et al. NIM A 737 (2014) 214, sec. 2.2
- Bonn GEM-TPC: Berger PhD 2015, sec. VIII.1
- T2K TPCs: Abgrall et al. NIM A 637 (2011) 25, sec. 7.1.1


Example Outcome
---------------

`julia clusternoisr.jl` prints the width (sigma) and mean of the residual distribution for the three options and saves the histograms as `clusternoisr_peakSNR<peak strip/amp_noise>_n<100 n_threshold>.png`.
The parameters at the top of the script can be overridden from the environment, e.g. `NUMBER_EVENTS=4000 AMP_SIGNAL=1200 AMP_NOISE=20 julia clusternoisr.jl`.

`julia scan.jl` (which includes `clusternoisr.jl`) scans the cluster width and the peak strip signal over noise instead, prints a table of core sigma and outlier fraction per cell, and saves `clusternoisr_scan_n<100 n_threshold>.png` (core sigma and outlier fraction against cluster width) and `clusternoisr_scurve_SNR<peak strip/amp_noise>_n<100 n_threshold>.png` (mean residual against the hit position within a strip).

The figures of the earlier version (exponential noise, global 3×noise cut, fixed hit position on a strip boundary, two options) are kept for reference:
![large amp - snr 30](https://user-images.githubusercontent.com/13000622/48422474-0402c580-e75f-11e8-8738-3a04c55b6b8d.png)
![small amp - snr 15](https://user-images.githubusercontent.com/13000622/48422905-dd915a00-e75f-11e8-9d16-8b36298fd535.png)

The spiky right-hand distributions of these figures are the discrete values the residual takes when the hit sits exactly on a strip boundary and one edge strip survives the cut while the other does not; with the hit position distributed over a strip they disappear.
Some more plots of that version are shown in issue #17.
