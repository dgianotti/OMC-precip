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

% The method follows the description in: 
% Daniel J. Short Gianotti, Bruce T. Anderson, and Guido D. Salvucci, 2014:
% The Potential Predictability of Precipitation Occurrence, Intensity, and
% Seasonal Totals over the Continental United States. J. Climate, 27,
% 6904–6918.
% 
% The code is hosted at: http://github.com/dgianotti/OMC-precip
% 
% You are welcome to use this code in any non-commercial context, however I
% ask that you please cite the code using the following DOI:
% 
%
% An example citation style is:
% 
% 
% If for some reason you are not able to cite source code, please cite the
% following paper: Daniel J. Short Gianotti, Bruce T. Anderson, and Guido
% D. Salvucci, 2014: The Potential Predictability of Precipitation
% Occurrence, Intensity, and Seasonal Totals over the Continental United
% States. J. Climate, 27, 6904–6918.

pooling = 31;
max_dims = 5;

if nargin > 1
    max_dims = varargin{1};
    assert(max_dims>=1 && round(max_dims)==max_dims)
    fprintf('Using max dims = %d.\n',max_dims);
end
if nargin > 2
    pooling = varargin{2};
    if pooling < 1 || mod(pooling,2)~=1
        error('OMC_fit expects the input pooling_size to be a positive, odd integer. Aborting!');
    end
    fprintf('Using pooling_size = %d.\n',pooling);
end


occ_obs = precip_obs>0;

selected_dims = zeros(1,365);
int_params = zeros(5,365);
occ_params = cell(365,1);

for i = 1:365
    fprintf('Fitting model for DOY %d...\n',i);
    AICs = inf(max_dims,1);
    
    for dim = 1:max_dims
        pooled_data = pool_data(occ_obs, pooling, dim-1, i);
        if dim == 1
            p_occ = 1-sum(pooled_data==0)/numel(pooled_data);
            n=1;
        else
            patterns = dec2bin(0:(2^(dim-1)-1))-'0'; % Each row is a pattern, like [0,0,0; 0,0,1; 0,1,0;... 1,1,1]
            n = size(patterns,1);
            p_occ = zeros(n,1);
            % Convert the pooled data to patterns:
            pooled_patterns = bin2dec(int2str(pooled_data(:,1:(end-1))));
            
            for k = 1:n
                p_occ(k) = sum(pooled_patterns==(k-1) & pooled_data(:,end)>0)...
                    / sum(pooled_patterns==(k-1));
            end
        end
        
        un_pooled = pool_data(occ_obs,1,dim-1,i);
        
        pdfs = zeros(size(un_pooled,1),1);
        
        
        if dim == 1
            pdfs(un_pooled(:,end)==1) = p_occ;
            pdfs(un_pooled(:,end)==0) = 1-p_occ;
        else
            for k = 1:n
                pdfs( all(bsxfun(@eq,patterns(k,:),un_pooled(:,1:(end-1))),2) & un_pooled(:,end)==1) = p_occ(k);
                pdfs( all(bsxfun(@eq,patterns(k,:),un_pooled(:,1:(end-1))),2) & un_pooled(:,end)==0) = 1-p_occ(k);
            end
        end
        LL = sum(log(pdfs));
        AICs(dim) = 2*n/2 - LL + n*(n/2+1)/(length(pdfs)-n/2-1);
        
        if AICs(dim) == min(AICs);
            selected_dims(i) = dim;
            occ_params{i} = p_occ;
        end
    end
    
    
    % Then re-pool the data and pick the double-gamma parameters:
    pooled_data = pool_data(precip_obs,pooling,0,i);
    int_data = pooled_data(pooled_data>0);
    int_params(:,i) = fit_mixed_gamma(int_data);
    
end

end % function

% % % % % % % % % % % % % % % % % % % %
function Y = pool_data(X,pl,lg,dy)
%POOL_DATA
% The function pool_data takes four inputs (X,pl,lg,dy), and returns a
% matrix of pooled data.
%
%   X is an N x D matrix of data, where N is the number of years and D is
%   the number of days per year;
%
%   pl is an odd positive integer, specifying the width of the pooling
%   window (i.e. pl = 7 means a pooling window of 7 days with the given day
%   at the center);
%
%   lg is a non-negative integer specifying the number of lags. If lg = 0,
%   then the returned Y will be a vector of days with no lags. If lg = 2,
%   then Y will have three columns (the first being lag 2, then lag 1, then
%   the present value);
%
%   dy is the day of the calendar year, and is an integer between 1 and D
%   (typically between 1 and 365). This is the day that is the center of
%   the pooling window.

% To deal with wrap-around problems (having lags from 100 years in the
% future...), we'll put a row of NaNs at one end of the data matrix and
% then throw out any pooled vectors with NaNs in them:

% Shift the data and extract the values we want:
shift_size = (pl+1)/2+lg-dy;
% Normally here I would use the awesome
%
% paren = @(x, varargin) x(varargin{:});
%
% to make use of pointer indexing and not copy the matrix so often,
% but I don't want to redefine it a billion times (since this is going
% to be called over and over again, so I'm going to make use of some
% arcane inner-Matlab voodoo:
X_subset = paren( ShiftXdays([X;nan(1,size(X,2))],shift_size),...
    ':', 1:(pl+lg) );
% X_subset is now an (N+1) x (pooling+lag) matrix.

% We need to turn it into something with lg+1 columns, this is a
% nightmare to read, but is vectorized at least:
cols = bsxfun(@plus, repmat((1:(lg+1))',...
    [1, size(X_subset,2)-lg]), ((1:(size(X_subset,2)-lg))-1)); % These are the column indices we want
Y = reshape(X_subset(:,cols)',[(lg+1),size(cols,2)*size(X_subset,1)])';

% Throw out the ones with NaNs (which would othwise wrap around from
% the first year to the last):
Y(any(isnan(Y),2),:) = [];
end

% % % % % % % % % % % % % % % % % % % %
function [params,neg_log_like] = fit_mixed_gamma(data)
% fit mixed gamma fits a mixture gamma distribution (2 weighted dists) to
% the data, requiring that the distribution's mean and variance fit the 
% data's mean and variance. The returned vector, params, is 5 x 1 in the
% following order: [k1; theta1; k2; theta2; weight1] where weight1 is the
% weight for the FIRST distribution (and weight2 = 1-weight1).

if isempty(data)
    params = nan(5,1);
    neg_log_like = inf;
    return
end

options = optimset('Algorithm','active-set','Display','notify','MaxFunEvals',2000);

params1 = gamfit(data);
var_obs = var(data);
mean_obs = mean(data);

guess_vector=[params1(1), params1(2), 0.5]; % our initial guess
lower_limit = [0.01, 0.01, 0.01];
upper_limit = [100, 100, 0.99];

[optimal_vec,~,exit_flag]=fmincon(@(x) gam_mix_likeC(x,data,mean_obs,var_obs),guess_vector,[],[],[],[],lower_limit,upper_limit,[],options);

if (exit_flag > 0) % Good news!
    [neg_log_like,k2,theta2]=gam_mix_likeC(optimal_vec,data,mean_obs,var_obs);
    params = [optimal_vec(1); optimal_vec(2); k2; theta2; optimal_vec(3)];
    return;
else % didn't converge!!!
    % try patternsearch! 
    
    % removing the patternsearch part because it uses the global
    % optimization toolbox. Can un-comment to use if toolbox is available
    
%     [optimal_vec,~,exit_flag]=patternsearch(@(x) gam_mix_likeC(x,data,mean_obs,var_obs),guess_vector,[],[],[],[],lower_limit,upper_limit);
%     fprintf('Patternsearch exit flag: %i.\n',exit_flag);
%     [neg_log_like,k2,theta2]=gam_mix_likeC(optimal_vec,data,mean_obs,var_obs);
%     params = [optimal_vec(1); optimal_vec(2); k2; theta2; optimal_vec(3)];
end

% If patternsearch worked, great, if not, try a little for loop for fmincon
% over weighting values:
if exit_flag <= 0 % Ugghh, this is slow!
    params = zeros(5,1);
    best_neg_log_like = inf;
    for i = 1:19
        weight = 0.05*i;
        [optimal_vec,~,exit_flag]=fmincon(@(x) gam_mix_likeC_fixed_weight(x,weight,data,mean_obs,var_obs),guess_vector,[],[],[],[],lower_limit,upper_limit,[],options);
        [neg_log_like,k2,theta2]=gam_mix_likeC([optimal_vec,weight],data,mean_obs,var_obs);
        if (neg_log_like < best_neg_log_like) && (exit_flag > 1)
            params = [optimal_vec(1); optimal_vec(2); k2; theta2; weight];
            best_neg_log_like = neg_log_like;
        end
    end
    fprintf('For loop neg_log_like gives %d.\n',best_neg_log_like);
    if ~isinf(best_neg_log_like)
        neg_log_like = best_neg_log_like;
        return;
    end
end

% If we're still going... then just fit a single gamma
fprintf('Using a simple gamma!\n');
params = [params1, 1, 1, 1]; % Set it to be weighted entirely to the first gamma dist
likelihood = gampdf(data,params1(1),params1(2));
neg_log_like = -sum(log(max(likelihood,1e-16)));

end % function

% % % % % % % % % % % % % % % % % % % %
function [neg_log_like,k2,theta2]=gam_mix_likeC(params1,data,mean_obs,var_obs)
% This determines the log likelihood from the data, given params1, a vector
% of three parameters: [k1, theta1, weight].

% Determine the likelihood of the data given the parameters for the first
% gamma distribution:
likelihood1 = gampdf(data,params1(1),params1(2));

mu1 = params1(1).*params1(2);
v1 = params1(1).*params1(2).^2;

wt = params1(3);

mu2 = (mean_obs-wt.*mu1)./(1-wt);

v2 = (var_obs-wt.*(mu1-mean_obs).^2-wt.*v1-(1-wt).*(mu2-mean_obs).^2)./(1-wt);

k2 = mu2.^2./v2; % k2
theta2 = v2./mu2; % theta2

likelihood2 = gampdf(data,k2,theta2);

likelihood = wt.*likelihood1+(1-wt).*likelihood2;

neg_log_like = -sum(log(max(likelihood,1e-16)));
end % function






