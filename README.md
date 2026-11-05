clusternoisr
============

A small script in Julia to investigate the impact of noise on cluster size and -resolution.
Since it is unclear (at least to me …), how "noise" has to be cut off correctly, I started this toy Monte Carlo simulation.

This README describes the problem, what the tool does and how to use it, and summarises the outcome.


Problem
-------

It is disputed how, in a particle detector, noise should be treated in cluster reconstruction.
The terms used here, per electronic channel (pad/strip/…) `i`:

* __pedestal__ `p_i`: the baseline of the channel, the mean reading without signal. It is measured per channel and subtracted in the front-end electronics; all amplitudes below are pedestal-subtracted, so the baseline is at zero.
* __noise__ `e_i`: the fluctuation of one reading around the baseline. It is Gaussian with mean zero and can be positive or negative; one measurement reads `d_i = s_i + e_i` for a signal `s_i`.
* __noise rms__ `σ_i`: the rms of `e_i`, measured from many readings without signal. It differs from channel to channel.
* __threshold__ `t_i = n σ_i`: the zero-suppression cut (`p_i + n σ_i` in raw ADC counts). `n` is typically between 3 and 5.

It is usual to remove all channels from analysis that have a reading (i.e. physics signal plus noise) below their threshold.
For the remaining channels, it is unclear whether the amplitude should enter the position calculation as is, or whether the noise, or the threshold, should be subtracted first.
Only the strips with `d_i ≥ t_i` are used in all cases; the options differ in the amplitude these strips enter the centre of gravity with:

1. The **signal over pedestal** `d_i`: nothing is subtracted from the surviving strips. This is what the readouts referenced below store.
2. The **signal minus 1σ** `d_i - σ_i`.
3. The **signal over threshold** `d_i - t_i = d_i - n σ_i`.

Intuitively, subtracting the threshold makes the weight of a strip go to zero continuously as it drops out of the cluster, while without subtraction a strip just above threshold enters with its full amplitude `t_i`; on the other hand, subtracting a constant from every strip distorts the weights of a narrow cluster.
Which effect wins is the question.


What the tool does
------------------

`clusternoisr.jl` simulates a plane of 255 strips:

* Each strip gets a __noise rms__ `σ_i`, drawn once per simulation from a Gauss around `amp_noise` with a spread of one third of that.
* Per event, a __hit__ is placed uniformly over the plane (excluding a margin at the edges). Its charge is a Gauss with standard deviation `cluster_width` (in strips) and total charge `amp_signal`; each strip receives the integral of the Gauss over its width. The strip under the hit centre gets the fraction `erf(1/(2√2 cluster_width))`: 11 % for `cluster_width = 3.5`, 38 % for 1, 95 % for 0.25. This *peak strip* amplitude over the noise rms, not `amp_signal/amp_noise`, decides whether the threshold cuts into the cluster.
* Each strip reads signal plus Gaussian noise with its `σ_i`; strips below `n_threshold σ_i` are set to zero.
* The surviving strips are weighted in the three ways above, and the __centre of gravity__ is calculated from the group of adjacent non-zero strips around the strip with the largest amplitude (COMPASS-style cluster finder: the cluster ends at the first strip below threshold).
* The __residual__ is the Monte Carlo truth position minus the centre of gravity, in strips.

Clusters that are split at a strip below threshold, or seeded on a noisy strip far away, give residuals of several to a hundred strips, which would dominate a plain rms.
Therefore the __resolution__ is reported as the rms of the residuals within ±2 strips, and the fraction outside as __lost clusters__.


Usage
-----

```
julia clusternoisr.jl
```
runs one simulation with the parameters at the top of the script and saves the residual histograms of the three options as `clusternoisr_peakSNR<peak strip/noise>_n<100 n_threshold>.png`.
The parameters can be overridden from the environment:
```
NUMBER_EVENTS=4000 AMP_SIGNAL=1200 AMP_NOISE=20 julia clusternoisr.jl
```
(`CLUSTER_WIDTH`, `N_THRESHOLD` and `SEED` likewise).

```
julia scan.jl
```
includes `clusternoisr.jl` and runs it over a grid of cluster widths and peak strip over noise (set at the top of `scan.jl`), prints a table of resolution and lost clusters per cell, and saves `clusternoisr_scan_n<100 n_threshold>.png` (resolution and lost clusters against cluster width) and `clusternoisr_scurve_SNR<peak strip/noise>_n<100 n_threshold>.png` (mean residual against the hit position within a strip).

Requires the Julia packages `Plots`, `SpecialFunctions`, `Statistics`, `Random`.


Results
-------

Residual histograms of the three options for a wide hit (σ = 3.5 strips) with the peak strip at 6.8 times the noise (left) and at 3.8 times the noise (right), 4000 events each:

![peak strip / noise = 6.8](https://github.com/user-attachments/assets/92a63ab1-077b-4d33-b346-2ccee182a9f4)
![peak strip / noise = 3.8](https://github.com/user-attachments/assets/b288abcc-bbfa-4e9b-a279-15850f04ce8c)

With the signal well above threshold, subtracting the threshold gives the best resolution (0.61 vs 0.67 strips) and the fewest lost clusters.
With the threshold cutting into the cluster, the resolution is the same for all options (1 strip), but one in five clusters is lost, fewest again with the threshold subtracted.

Resolution and lost clusters against the cluster width, for three values of the peak strip over noise (`scan.jl`, 10000 events per point):

![resolution and lost clusters](https://github.com/user-attachments/assets/bd3cdb33-ee0a-4ba9-aa92-b94be4fc15d9)

Mean residual against the hit position within a strip (0 = strip edge, 0.5 = strip centre), the bias of the centre of gravity:

![mean residual against hit position](https://github.com/user-attachments/assets/c57ffad3-8c0d-4c4a-8612-41f27a78b908)

In brief:

* Subtracting the __threshold__ loses the fewest clusters in every setting, and gives the best resolution for clusters wider than about one strip. It is the robust choice.
* Subtracting __nothing__ is fragile: a noisy channel just above its own (high) threshold enters with its full amplitude and can seed the cluster. For clusters narrower than a strip it can give the better resolution.
* Subtracting __1σ__ is never the best and never the worst option.
* For sub-strip clusters, the best weighting depends on the cluster width at the ±20 % level: subtracting the threshold takes most of the weight off a neighbour strip that sits just above it.

The old version of the tool (exponential noise, global 3×noise cut, fixed hit position on a strip boundary, two options) is discussed in issue #17.


Relevant Read
-------------

- COMPASS GEM trackers: Altunbas NIM A 490 (2002) 177, sec. 7
- ALICE TPC: Alme et al. NIM A 622 (2010) 316, sec. 10.3.1
- FOPI GEM-TPC: Boehmer et al. NIM A 737 (2014) 214, sec. 2.2
- Bonn GEM-TPC: Berger PhD 2015, sec. VIII.1
- T2K TPCs: Abgrall et al. NIM A 637 (2011) 25, sec. 7.1.1
