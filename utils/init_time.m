function [ init_time ] = init_time( transientTime, sampleSize )
%WEB_OFF Init time of traffic source
%   Initial time of traffic source following a uniform distr.

%Output in BYTES
if nargin == 1
    sampleSize=1;
end

init_time = unifrnd(0,transientTime,sampleSize,1);

end

