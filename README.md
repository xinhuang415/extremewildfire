# extremewildfire
NCL and matlab code for the analysis of the extreme wildfires.

This code uses satellite retrievals on daily burned areas to identify the synoptic-scale variation of wildfires in one specific fire-prone region and then analyze the temporal variation of key meteorological parameters derived from ERA5 reanalysis data. If you have any comments or suggestions regarding the code, please share them with us.

Empirical mode decomposition (EMD) was first applied to daily burned areas over one fire-prone region. After decomposing the long-term time series into a set of intrinsic mode functions (IMF) of burned area, the power spectrum of each IMF by the fast Fourier transform (FFT) method is calculated. For the purpose of extracting the periodicity at the synoptic scale, the power spectrum was obtained by averaging the IMFs with periods less than 20 days. 

To investigate the temporal evolution of key meteorological parameters around extreme wildfires, which is identified by daily burned area greater than 90th percentile, the fifth generation ECMWF atmospheric reanalysis (ERA5) reanalysis data is collected. Then, two-week variations of main meteorological parameters in one fire-prone region, including relative humidity, wind speed, and precipitation around the extreme-fire day are calculated.

Software requirements:
NCL code tested with NCL version 6.4 and matlab code tested with R2020.
