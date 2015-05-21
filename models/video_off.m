function [ off_time ] = video_off( ton_duration )
%VIDEO_OFF Calculates the time off of an VIDEO model
%   *** Calculates the time off (packet interarrival time) of an 
%   VIDEO model given a certain Ton (previous activity duration).
%
% USAGE:
%       [ off_time ] = video_off( ton_duration )
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

seconds_per_frames = 1;

off_time = seconds_per_frames - ton_duration;

if (off_time < 0)
    off_time = 0;
end

end

