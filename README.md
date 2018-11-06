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
