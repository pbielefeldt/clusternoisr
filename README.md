clusternoisr
============

A small script in Julia to investigate the impact of noise on cluster size and -resolution.
Since it is unclear (at least to me …), how "noise" has to be cut off correctly, I started this toy Monte Carlo simulation.


Problem
-------

It is disputed how, in a particle detector, noise should be treated in cluster reconstruction.
Usually, all electronic channels (pads/strips/…) of the detector have __noise__ of a certain degree. 
Furthermore, they have __pedestals__ (baselines), which are measured per channel and subtracted in the front-end electronics; what remains is a Gaussian fluctuation with mean zero whose rms is the noise of the channel (Vandenbroucke, PhD 2012, sec. 4.5.1; Berger et al., NIM A 869 (2017) 180, sec. 12.2; Alme et al., NIM A 622 (2010) 316, sec. 10.3.1).
It is usual to remove all channels from analysis that have signal (i.e. physics signal plus noise) smaller than a multiple `n` of the noise rms of the channel; `n` is 3 for the COMPASS GEM trackers and the ALICE TPC, 3.5 for the Bonn GEM-TPC cosmics (Berger, PhD 2015, sec. VIII.1) and 4.5 for the T2K TPCs and the FOPI GEM-TPC (Abgrall et al., NIM A 637 (2011) 25, sec. 7.1.1; Böhmer et al., NIM A 737 (2014) 214, tab. 1).
Now, for the remaining channels, it is unclear whether or not the nσ cut applied to the signal should also be applied to the analysed data. 
If done so, the remaining signal is smoother. 
Else, the edges will be sharper.
However, since the nσ cut is somewhat arbritrary, many argue it is more accurate to only cut 1σ (the assumed noise) from the signal. 
It can also be argued that the entire (pedestal-corrected) signal should be taken for finding the centre of gravity; that is what the readouts referenced above do.


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

(The figure below is from the earlier version of the script with exponential noise, a global cut and two options; it shows the arrays for a high SNR of 100 so that the peak is overly visible.)
![all arrays for large peak](https://user-images.githubusercontent.com/13000622/48423934-e2570d80-e761-11e8-95ed-b877430366f2.png)


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
