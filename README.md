clusternoisr
============

A small script in Julia to investigate the impact of noise on cluster size and -resolution.
Since it is unclear (at least to me …), how "noise" has to be cut off correctly, I started this toy Monte Carlo simulation.


Problem
-------

It is disputed how, in a particle detector, noise should be treated in cluster reconstruction.
Usually, all electronic channels (pads/strips/…) of the detector have __noise__ of a certain degree. 
Furthermore, they have __pedestals__ (baselines), which are measured per channel and subtracted in the front-end electronics; what remains is a Gaussian fluctuation with mean zero whose rms is the noise of the channel.
It is usual to remove all channels from analysis that have signal (i.e. physics signal plus noise) smaller than a multiple `n` of the noise rms of the channel; with `n` anywhere between 3 and 5 typically.

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

There are three options:
1. The **signal over threshold** (`pedestal + nσ`) is used to calculate the hit's mean position (either by fitting a function to the histogram, or by averaging for a centre of gravity)
2. The **signal over pedestal** is used to calculate the position: all strips with signal > threshold are used, but nσ is not subtracted from (or re-added to) their values
3. The **actual signal** is used, i.e. all strips with signal < threshold are removed, and then the mean is calculated from the signal value (without any subtraction, or with pedestal + nσ re-added) 


Input
-----

A 2D __strip plane__ is assumed as geometry (default: 255 strips wide).
Each strip is associated with a randomly chosen __noise level__ (its noise rms; per simulation, not per event), that is gaussianly shaped around `amp_noise`.
Take `number_events` runs in the loop. 
In each, do the following:
 * For each strip, calculate an event-specific noise. That is *Gaussian* with mean zero and the noise level of this strip as rms, i.e. the pedestal-corrected reading of a channel without signal. 
 * Place a hit at a position drawn uniformly over the strip plane (excluding a margin at the edges). It has a *gauss shape*, with an amplitude of `amp_signal` and a width of `cluster_width`.
 * Add up signal and Gaussian noise (`data_arr = enoise_arr.+signal_arr`)


Noise Suppression
-----------------

All strips with an amplitude below `n_threshold` times their own noise level are set to zero (zero suppression). The surviving strips are then treated in three ways:

__Option 1__ ("only zeros suppressed"): keep the pedestal-corrected amplitude. This is what the referenced readouts store.

__Option 2__ ("1σ subtracted"): subtract one noise level of the strip.

__Option 3__ ("noise corrected"): subtract the full threshold, `n_threshold` noise levels of the strip.


Cluster Finding
---------------

As can be seen from the figure above, there are some bins with entries that apparently do not belong to the signal.
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

The script prints the width (sigma) and mean of the residual distribution for the three options and saves the histograms as `clusternoisr_SNR<amp_signal/amp_noise>_n<n_threshold>.png`.

The figures of the earlier version (exponential noise, global 3×noise cut, fixed hit position on a strip boundary, two options) are kept for reference:
![large amp - snr 30](https://user-images.githubusercontent.com/13000622/48422474-0402c580-e75f-11e8-8738-3a04c55b6b8d.png)
![small amp - snr 15](https://user-images.githubusercontent.com/13000622/48422905-dd915a00-e75f-11e8-9d16-8b36298fd535.png)

The spiky right-hand distributions of these figures are the discrete values the residual takes when the hit sits exactly on a strip boundary and one edge strip survives the cut while the other does not; with the hit position distributed over a strip they disappear.
Some more plots of that version are shown in issue #17.
