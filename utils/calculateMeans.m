function [ mean_per_row, mean_per_row_active ] = calculateMeans( data, timesteps )
%CALCULATEMEANS Calculates the mean of discrete data
%   DCalculates the mean of discrete data taking into account the period (time step) for each value


if (size(data,2)~=length(timesteps))
   if ( size(data,2)==1 && size(data,1)==length(timesteps) )
       %Hypothesis: only one user
       data = data';
   else
       fprintf('\nError output for traceability: \n');
       data
       timesteps
       error('Data columns should agree with the length of timesteps.');
   end
end

timesteps_normalized = timesteps/sum(timesteps);

%Mean
mean_per_row = sum(bsxfun(@times,data,timesteps_normalized),2);

%Mean during activity period
mean_per_row_active=zeros(size(data,1),1);
timesteps_filtered = bsxfun(@times,timesteps,data~=0);
%Avoid NaNs for users without consumption
ind_not_nan = sum(timesteps_filtered,2)~=0;
mean_per_row_active(ind_not_nan) = sum(data(ind_not_nan,:).*bsxfun(@rdivide,timesteps_filtered(ind_not_nan,:),sum(timesteps_filtered(ind_not_nan,:),2)),2);

% %Alternative to previous lines:(if mean_per_row==0, then mean_per_row_active=0)
% aux_nan_ind = mean_per_row==0;
% mean_per_row_active(aux_nan_ind)=0;

end

