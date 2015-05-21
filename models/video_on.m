function [ total_size ] = video_on(  )
%VIDEO_ON Summary of this function goes here
%   Detailed explanation goes here

fps=25;
a=1.2;
L=4800;
H=26100;

uniform_samples = rand(fps,1);

frame_samples = (-( (uniform_samples*(H^a)-uniform_samples*(L^a)-(H^a)) / ((H^a)*(L^a)) )).^(-1/a);

total_size = sum(frame_samples);

end

