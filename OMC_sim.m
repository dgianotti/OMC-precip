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


max_lags = max(selected_dims)-1;
hist = rand(N_sims,max_lags) > mean(occ_params{1}(:));
precip_sim = zeros(N_sims,N_years*365 +365); % with the extra +365 as a first year to throw out
noise = rand(size(precip_sim));

for i = 1:size(precip_sim,2)
    doy = mod(i-1,365)+1;
    
    if doy == 1
        fprintf('Simulating data for year %d...\n',floor(i/365)+1);
    end
    
    if selected_dims(doy) == 1
        precip_sim( noise(:,i) < occ_params{doy}  , i) = 1;
    else % Need to look at history:
        patterns = dec2bin(0:(2^(selected_dims(doy)-1)-1))-'0'; % Each row is a pattern, like [0,0,0; 0,0,1; 0,1,0;... 1,1,1]
        n = size(patterns,1);
        for k = 1:n
            pattern_matches = all( bsxfun(@eq,hist(:, (end-selected_dims(doy)+2):end),patterns(k,:)),2);
            precip_sim( pattern_matches & noise(:,i) < occ_params{doy}(k), i) = 1;
        end
        hist = [hist(:,2:end),precip_sim(:,i)];
    end
end

%% Add in intensity!

% First we want to throw out the first year and reshape precip_sim to be
% (N_years*1000) x 365:
precip_sim = reshape(precip_sim(:,366:end)',[365,N_sims*N_years])';
n1 = size(precip_sim,1);

for i = 1:365
    % Assign which gamma distribution gets selected:
    noise = rand(n1,2);
    % For those in gamma1:
    precip_sim( (noise(:,1) < int_params(5,i) ) & (precip_sim(:,i)==1) , i)...
        = gaminv( noise( (noise(:,1) < int_params(5,i) ) & (precip_sim(:,i)==1) , 2), int_params(1,i), int_params(2,i) );
    % For those in gamma2:
    precip_sim( (noise(:,1) >= int_params(5,i) ) & (precip_sim(:,i)==1) , i)...
        = gaminv( noise( (noise(:,1) >= int_params(5,i) ) & (precip_sim(:,i)==1) , 2), int_params(3,i), int_params(4,i) );   
end

end % function