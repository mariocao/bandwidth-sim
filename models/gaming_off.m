function [ off_time ] = gaming_off( ton_duration )
%GAMING_OFF Calculates the time off of an ONLINE GAMING model
%   *** Calculates the time off (packet interarrival time) of an 
%   ONLINE GAMING model given a certain Ton (previous activity duration).
%
% USAGE:
%       [ off_time ] = gaming_off( ton_duration )
%
% INPUT:
%       - ton_duration:     the time of the previous acitivity period
%
% OUTPUT:
%       - off_time:         the interarrival time in seconds
%
% Author: Mario Cao Cueto
% DIT-UPM, Universidad Politécnica de Madrid
% eMail: mcao@dit.upm.es
% Copyright 2014 Mario Cao Cueto.

% ONLINE MODEL:
% Srinivasan, R., Zhuang, J., Jalloul, L., Novak, R., & Park, J. (2008). 
% IEEE 802.16 m evaluation methodology document (EMD).
% http://ieee802.org/16/tgm/docs/80216m-08_004r2.pdf


%% ONLINE GAMING: Toff (Packet IAT)

off_time_aux = evrnd(50,4.5)*1e-3;

% Avoid 'strange' negative values.
while(off_time_aux<0)
    off_time_aux = evrnd(50,4.5)*1e-3;
end

% If Toff-Ton is smaller than 0, an online gaming update is required.
off_time = off_time_aux - ton_duration;
if (off_time < 0)
    off_time = 0;
end

end

