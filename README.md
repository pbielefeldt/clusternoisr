clusternoisr
============

A small script in Julia to investigate the impact of noise on cluster size and -resolution.
Since it is unclear (at least to me …), how "noise" has to be cut off correctly, I started this toy Monte Carlo simulation.


Problem
-------

It is disputed how, in a particle detector, noise should be treated in cluster reconstruction.
Usually, all electronic channels (pads/strips/…) of the detector have __noise__ of a certain degree. 
Furthermore, they have __pedestals__, which are not analysed here.
It is usual to remove all channels from analysis that have signal (i.e. physics signal plus noise) smaller than three times the assumed noise on each channel. 
Now, for the remaining channels, it is unclear whether or not the 3σ cut applied to the signal should also be applied to the analysed data. 
If done so, the remaining signal is smoother. 
Else, the edges will be sharper.
However, since the 3σ cut is somewhat arbritrary, many argue it is more accurate to only cut 1σ (the assumed noise) from the signal. 
It can also be argued that the entire signal should be taken for finding the centre of gravity.


Input
-----

A 2D __strip plane__ is assumed as geometry (default: 255 strips wide).
Each strip is associated with a randomly chosen __noise level__ (per simulation, not per event), that is gaussianly shaped around `amp_noise`.
Take `number_events` runs in the loop. 
In each, do the following:
 * For each strip, calculate an event-specific noise. That is *exponential* around the noise level of this strip. 
 * Roughly in the centre of the strip plane, place a hit. It has a *gauss shape*, with an amplitude of `amp_signal` and a width of `hit_sigma`.
 * Add up signal and exponential noise (`data_arr = enoise_arr.+signal_arr`)


Noise Suppression
-----------------

The `data_arr` is treated in two ways:

__Option 1__: Subtract three times the noise level for each strip. Strips with an amp<0 are set to amp=0. This is called "noise corrected".

__Option 2__: Only set channels to zero that have an amplitude less than three times the noise level. Hence, *only* do zero suppression. 

The outcome is shown in this figure, taken with a high SNR of 100 (so that the peak is overly visible):
![all arrays for large peak](https://user-images.githubusercontent.com/13000622/48423934-e2570d80-e761-11e8-95ed-b877430366f2.png)


Cluster Finding
---------------

As can be seen from the figure above, there are some bins with entries that apparently do not belong to the signal.
Such things of course always happen in real data.
In this very comprehensive Monte Carlo study, a simplistic approach is used to handle the situation:
It is assumed that the central bin `N = Int(trunc(length(arr)/2))` always is within the cluster.
From this bin, the function `trunc_arr(arr)` iterates to the left and right, until it hits a bin with amplitude zero. 
All bins beyont this bin are discarded.


Example Outcome
---------------

For a well-behaved situation with `amp_signal = 300` and `amp_noise = 10`, it looks like this:
![large amp - snr 30](https://user-images.githubusercontent.com/13000622/48422474-0402c580-e75f-11e8-8738-3a04c55b6b8d.png)

The case for the smallest amp I could simulate (`amp_signal = 150`) yields this image:
![small amp - snr 15](https://user-images.githubusercontent.com/13000622/48422905-dd915a00-e75f-11e8-9d16-8b36298fd535.png)

Some more plots are shown in issue #17.
