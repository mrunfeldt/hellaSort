# hellaSort
Semi-automatic neuronal spike sorting GUI written in MATLAB. 

This interactive GUI permits the user to explore their waveforms and either use a gaussian mixture model to cluster spikes or manually cluster in a 2D PCA space. 

The user can vary the number of clusters and assess the impact of the chosen clusters on various metrics, including: the evolution of several waveform parameters over time, ISI distribution,  pairwise correlations, crosscorrelations, quality metrics (Kleinfeld et al.,2011), average waveform +/- std, and binned spike rate over time.

The main "hellaSort.m" function calls a number of modules/functions that can be independently modified. 

To use the software, download and place into your MATLAB directory: 
(1) subroutines (i.e. functions) from the "subroutines" branch <https://github.com/mrunfeldt/hellaSort/tree/master/subroutines> and place them within your MATLAB directory. 
(2) the main functions "hellaSort.m" and "raw2nev.m"

"RUNME_load_sort_exampleData.m" is a script that you can use to run hellaSort on some example neural data ("exampleData.mat") that is provided.

Have fun and post modifications!

Written by Melissa Runfeldt, initiated March 2015 in the lab of Brian Malone at UCSanFransisco
