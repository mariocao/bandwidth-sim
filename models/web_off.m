function [ reading_time ] = web_off( N, M )
%WEB_OFF Summary of this function goes here
%   Detailed explanation goes here

if nargin == 0
    N=1;
    M=1;
elseif nargin <= 1
    M=1;
end
    
%READING TIME
reading_time_max=10000; 
reading_time=lognrnd(-0.495204,2.7731,N,M);
reading_time=min(reading_time,reading_time_max);

end

