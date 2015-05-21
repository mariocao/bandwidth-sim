function [ results, sorted_data ] = analyzeAggrBw( data, timesteps )
%ANALYZEPROBABILITIES Analyses the aggregated bandwidth
%   Analyses the aggregated bandwidth taking into account the period for each data vlue

if (nargin ~=2)
    error('usage:  analyzeAggrBw( capacities, timesteps )');
end

results = zeros(1,length(data));

[sorted_data, indices] = sort(data,'ascend');
sorted_timesteps = timesteps(indices);

cum_times = cumsum(sorted_timesteps);

cum_percent = cum_times/sum(sorted_timesteps);

results = ((sum(sorted_timesteps)-(cum_times))/sum(sorted_timesteps));

% for i=1:length(probabilities)
%    indices = find(cum_percent>=probabilities(i));
%    if ( ~ isempty(indices) )
%        results(i) = sorted_caps(indices(1));
%    else
%        results(i) = 0;
%    end
% end

% fprintf('-> Capacity 100%%: %.2f Mbps \n',r_mean(100)/1e6);
% fprintf('-> Capacity 98%%: %.2f Mbps \n',r_mean(98)/1e6);
% fprintf('-> Capacity 95%%: %.2f Mbps \n',r_mean(95)/1e6);


end

