function [ array ] = distributeValues( values, percentages, N )
%DISTRIBUTEVALUES Assign specific values to an array given its percentages.
%   ***Calculates the values of an array of size N, which appear following
%   the percentages given.
% 
% USAGE:
%   [ array ] = distributeValues( values, percentages, N )
%
% INPUT:
%   values      - values of the array
%   percentages - appearance porcentages of the values
%   N           - size of the result array
%
% OUTPUT:
%   array       - resulting array


if (nargin ~= 3)
    error('usage:  distributeValues( values, percentages, N )');
end

if (length(values)~=length(percentages))
    error('Values and percentages dimensions should be equal');
end

array = zeros(N,1);
counter=0;
for i=1:(length(values)-1),
    elementsNumber=round(N*percentages(i));
    array(counter+1:counter+elementsNumber)=values(i);
    counter=counter+elementsNumber;
end
array(counter+1:end)=values(i+1);
array = array(randperm(size(array,1)),:);
    
end

