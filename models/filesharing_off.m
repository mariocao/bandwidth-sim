function [ off_time ] = filesharing_off( sampleSize )
%FILESHARING_OFF Calculates the time off of a FILE SHARING model
%   *** Calculates the time off (packet interarrival time) of a 
%   FILE SHARING model.
%
% USAGE:
%       [ off_time ] = filesharing_on( sampleSize )
%
% INPUT:
%       - sampleSize:       size of output off times
%
% OUTPUT:
%       - off_time:         the inactivity time
%
% Author: Mario Cao Cueto
% DIT-UPM, Universidad Politécnica de Madrid
% eMail: mcao@dit.upm.es
% Copyright 2014 Mario Cao Cueto.


%% FILE SHARING: Toff
if nargin == 0
    sampleSize=1;
end

off_time = zeros(sampleSize,1);

end

