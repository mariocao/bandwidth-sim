function [ total_size ] = gaming_on(  )
%GAMING_ON Calculates the bytes to consume given by an ONLINE GAMING model
%   *** Calculates the bytes to be consume during the activity time (Ton)
%   of an ONLINE GAMING model.
%
% USAGE:
%       [ total_size ] = gaming_on( )
%
% INPUT:
%       - sampleSize:   the size of the sample to be 
%
% OUTPUT:
%       - total_size:   the packet size (in bytes) during activity period
%
% Author: Mario Cao Cueto
% DIT-UPM, Universidad Politécnica de Madrid
% eMail: mcao@dit.upm.es
% Copyright 2014 Mario Cao Cueto.

% ONLINE MODEL:
% Srinivasan, R., Zhuang, J., Jalloul, L., Novak, R., & Park, J. (2008). 
% IEEE 802.16 m evaluation methodology document (EMD).
% http://ieee802.org/16/tgm/docs/80216m-08_004r2.pdf


%% ONLINE GAMING: BYTES during Ton

server_packet = 0;
while(server_packet<=0)
    server_packet = evrnd(330,82);
end

total_size=server_packet;

end

