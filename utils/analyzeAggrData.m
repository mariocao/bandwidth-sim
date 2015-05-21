function [ results, sorted_data ] = analyzeAggrData( data, timesteps )
%ANALYZEPROBABILITIES Analyses the aggregated data
%   Analyses the aggregated data taking into account the period for each data value

if (nargin ~=2)
    error('usage:  analyzeAggrBw( capacities, timesteps )');
end

[sorted_data, indices] = sort(data,'ascend');
sorted_timesteps = timesteps(indices);

results = sorted_timesteps/sum(sorted_timesteps);

end

