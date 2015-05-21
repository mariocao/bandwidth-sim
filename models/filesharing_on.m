function [ size ] = filesharing_on( sampleSize )
%FILESHARING_OFF Calculates the bytes during the on time of a FILE SHARING model
%   *** Calculates the bytes to be consume during the activity time (Ton)
%   of a FILE SHARING model.
%
% USAGE:
%       [ size ] = filesharing_on( sampleSize )
%
% INPUT:
%       - sampleSize:   size of output off times
%
% OUTPUT:
%       - size:         size of the resources to be consumed
%
% Author: Mario Cao Cueto
% DIT-UPM, Universidad Politécnica de Madrid
% eMail: mcao@dit.upm.es
% Copyright 2014 Mario Cao Cueto.


%% FILE SHARING: BYTES during Ton
%Output in BYTES
if nargin == 0
    sampleSize=1;
end

% Infinite file sizes
size = Inf(sampleSize,1);

end

