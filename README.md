# OMC package: An Occurrence Markov Chain model for precipitation
The OMC package consists of two Matlab function: OMC_fit and OMC_sim to be used in 
simulating daily precipitation occurrence and intensity. Occurrence follows a 
variable-order Markov Chain, and intensity follows a gamma-gamma mixture model. 

The method follows the description in:
Daniel J. Short Gianotti, Bruce T. Anderson, and Guido D. Salvucci, 2014: 
The Potential Predictability of Precipitation Occurrence, Intensity, and 
Seasonal Totals over the Continental United States. J. Climate, 27, 6904–6918.

The code is hosted at:
http://github.com/dgianotti/OMC-precip

You are welcome to use this code in any non-commercial context, however
I ask that you please cite the code using the following DOI:

An example citation style is:


If for some reason you are not able to cite source code, please cite the following 
paper:
Daniel J. Short Gianotti, Bruce T. Anderson, and Guido D. Salvucci, 2014: 
The Potential Predictability of Precipitation Occurrence, Intensity, and 
Seasonal Totals over the Continental United States. J. Climate, 27, 6904–6918.


This software is licensed under the GNU General Public License v3.0.
See the accompanying License.txt file for details.


## Documentation for OMC_fit
function [selected_dims, occ_params, int_params] = OMC_fit(precip_obs,varargin)
%OMC_FIT
%
% OMC_fit(precip_obs) fits an Occurrence-based Markov Chain model to the
% precipitation data precip_obs. precip_obs should be an N_years x 365 days
% matrix of data, with the first year in the first row, and the first day
% of the year in the first column (usually Jan 1, but specific dates aren't
% important as long as the days are sequential -- water years work fine,
% for example).
%
% OMC_fit returns three variables for use in analysis and/or simulation of
% precipitation data using OMC_sim:
%
% selected_dims [1 x 365 vector]. Each entry of selected dims corresponds
% to a given column of precip_obs, and the values represent the optimal
% memory to maintain the auto-correlation in precipitation occurrence for
% that day of the year. For example, if selected_dims(4) = 3, then the
% dimensionality of the 4th day of the year (presumably Jan 4) is 3 days,
% or equivalently 2 lags; so occurrence data from the previous two days
% (Jan 2-3) is needed to appropriately determine the probability of
% precipitation occurrence on Jan 4.  The conditional probabilities of
% occurrence for each day of the year are calculated separately (although
% using pooled data for better estimation and more smoothly-varying
% parameter values), and together the 365 models constitute a
% variable-order Markov Chain in occurrence.
%
% occ_params [365 x 1 cell array]. Each entry in the occ_params cell array
% gives the conditional probability of precipitation occurrence given a
% specific occurrence history. For example, if selected_dims(4) = 3, then
% there are two lagged values of occurrence used to determine the
% probability of occurrence for this day, and the possible histories of
% occurrence are 00, 01, 10, 11, where 00 means no rain either of the last
% two days, 01 means no rain two days ago and rain yesterday, etc. The
% values saved in occ_params{4} then are 
% [p001|00x; p011|01x; p101|10x;
% p111|11x], where p001|00x means then probability of occurrence given two
% dry days, etc. Conditional probabilities are determined from the pooled
% data in a window around the given day of year.
%
% int_params [5 x 365 matrix]. Each column of int_params gives the
% parameters for a mixture model of two gamma distributions. This type of
% model was generally selected as the most appropriate by AICc for modeling
% precipitation from the USHCN daily precipitation dataset for 774 stations
% from [cite: Daniel J. Short Gianotti, Bruce T. Anderson, and Guido D.
% Salvucci, 2014: The Potential Predictability of Precipitation Occurrence,
% Intensity, and Seasonal Totals over the Continental United States. J.
% Climate, 27, 6904–6918.] Each column is in the following order: [k1;
% theta1; k2; theta2; weight1] where weight1 is the weight for the FIRST
% distribution (and weight2 = 1-weight1).
%
%
% OMC_fit(precip_obs, max_dim) overrides the maximum dimensionality tested
% (selected using the corrected AICc) with the integer value max_dim.
% max_dim must be greater than 0 [default value 5 days, equivalent to 4
% lags].
%
%
% OMC_fit(precip_obs, max_dim, pooling_size) overrides the pooling window
% size used for parameter estimation [default vaule 31 days]. pooling_size
% must be positive and odd. A pooling size of 1 uses a single day's
% observed values in estimating occurrence and intensity parameters for the
% model for that day.  Pooling sizes ranging from 15-45 days seem to be
% appropriate typically for 80-120 years of daily data [cite: Daniel J.
% Short Gianotti, Bruce T. Anderson, and Guido D. Salvucci, 2014: The
% Potential Predictability of Precipitation Occurrence, Intensity, and
% Seasonal Totals over the Continental United States. J. Climate, 27,
% 6904–6918.]. Too large of a window will limit the amplitude of
% seasonality, and too small of a window will lead to poor calibration of
% model parameters.
%
%
% OMC_fit requires the statistical toolbox and the optimization toolbox.



## Documentation for OMC_sim
function precip_sim = OMC_sim(selected_dims, occ_params, int_params, N_years, N_sims)
%OMC_sim
%
% OMC_sim(selected_dims, occ_params, int_params, N_years, N_sims) simulates
% an ensemble of N_sims simulations, each N_years long, using the model
% parameters returned by OMC_fit. Typical usage is:

% > [selected_dims, occ_params, int_params] = OMC_fit(precip_obs);
%
% > precip_sim = OMC_sim(selected_dims, occ_params, int_params,
% size(precip_obs,1), 1000);
%
% precip_sim is an N_years * N_sims by 365 matrix of data, with each
% N_years rows representing a separately simulated ensemble member.
% Different ensemble members should not have any daily correlation from the
% end of one run to the begining of the next, so
% OMC_sim(p_obs,dims,op,ip,100,10) will be slightly different statistically
% than OMC_sim(p_obs,dims,op,ip,1000,1), but only in terms of the
% auto-correlation at the break between simulations, not in terms of mean
% or interannual variance.
% 
% OMC_sim uses one year as burn-in for the random-number generator and
% occurrence auto-correlation chain, but this is removed prior to returning
% to the user.
%
