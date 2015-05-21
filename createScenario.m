function [simScenario, bandwidths, data] = createScenario(N_households, concurrency, modelConfigFile, DEBUG_LEVEL, RANDOMIZE)
%CREATESCENARIO Creates a pseudo-random scenario from input parameters
%   The scenario is created from the input parameters of the number of
%   households and some models specified in the given file.


%% INIT

if (nargin < 3 || nargin > 5 )
    error('usage:  createScenario(N_households, concurrency, modelConfigFile, DEBUG_LEVEL, RANDOMIZE)')
elseif nargin ==3
    DEBUG_LEVEL = 0;
    RANDOMIZE = 0;
elseif nargin == 4
    RANDOMIZE = 0;
end

FG_SERVICE_ID = 1;
BG_SERVICE_ID = 2;

%Load model parameters
fid = fopen(modelConfigFile);
vars = textscan(fid, '%s', 'Delimiter', '', 'CommentStyle', '%');
cellfun(@evalc, vars{1}, 'UniformOutput', false);
fclose(fid);

%Add utils path
addpath('utils');


%% APPLICATION OF SCENARIO MODELs

% Random columns of a matrix:
% test(randperm(size(test,1)),:)
% Random rows of a matrix:
% test(:,randperm(size(test,2)))

% ------------------------------------------------------------------------
%  Users per household
% ------------------------------------------------------------------------
% Hypothesis: Capacity in households (last mile) does not limit user's bandwidth
if (RANDOMIZE)
    %OLD version: random number of users per household
    users_in_households = randsample((1:length(users_per_households_ratios))', N_households, true, users_per_households_ratios);
    total_users = sum(users_in_households);
else
    %total_users = ceil(N_households*users_per_household_ratio);
    users_in_households = distributeValues(1:length(users_per_households_ratios), users_per_households_ratios, N_households);
    total_users = sum(users_in_households);
end

% Household IDs for each user
user_household_id = zeros(total_users,1);
index=0;
for i=1:length(users_in_households),
    user_household_id(index+1:index+users_in_households(i)) = i;
    index=index+users_in_households(i);
end


% ------------------------------------------------------------------------
%  Bandwidth per household
% ------------------------------------------------------------------------
if (RANDOMIZE)
    household_bandwidth = randsample( (subscriber_bandwidth)', N_households, true, subscriber_ratio);
    users_bandwidth =zeros(total_users,1);
    j=1;
    for i=1:length(household_bandwidth),
        users_bandwidth(j:j+users_in_households(i)-1)=household_bandwidth(i);
        j=j+users_in_households(i);
    end
    clear i j;
else
    
%OLD
%     users_bandwidth = zeros(total_users,1);
%     counter=0;
%     for i=1:(length(subscriber_bandwidth)-1),
%         elementsNumber=round(total_users*subscriber_ratio(i));
%         users_bandwidth(counter+1:counter+elementsNumber)=subscriber_bandwidth(i);
%         counter=counter+elementsNumber;
%     end
%     users_bandwidth(counter+1:end)=subscriber_bandwidth(i+1);
%     users_bandwidth(randperm(size(users_bandwidth,1)),:);

%CURRENT
%    users_bandwidth = distributeValues(subscriber_bandwidth,subscriber_ratio,total_users); %MARIO CAMBIAR!!! hacerlo por casas y no usuarios!
    
%NEW
    household_bandwidth = distributeValues(subscriber_bandwidth,subscriber_ratio,N_households);
    users_bandwidth = zeros(total_users,1);
    j=1;
    for i=1:length(household_bandwidth),
        users_bandwidth(j:j+users_in_households(i)-1)=household_bandwidth(i);
        j=j+users_in_households(i);
    end
    clear i j;
    
end


% ------------------------------------------------------------------------
%  User Typology (0/1/2/3/4)
% ------------------------------------------------------------------------

if (RANDOMIZE)
    %Discard Non-Users with Internet Access --> 0:Non-Users, 1:Internet Users
    internet_users = randsample( [0;1], total_users, true, [nonusers_ratio 1-nonusers_ratio]);
    internet_users_indices = internet_users(:,1)== 1;

    % User Profiles -->  1/2/3/4; 0:Non-User
    profiles_users = zeros(total_users,1);
    profiles_users(internet_users_indices) = randsample(4, sum(internet_users_indices), true, profiles_ratios);
else
    %Discard Non-Users with Internet Access --> 0:Non-Users, 1:Internet Users
    internet_users = distributeValues([0:1],[nonusers_ratio 1-nonusers_ratio],total_users);
    internet_users_indices = internet_users(:,1)== 1;
    
    % User Profiles -->  1/2/3/4; 0:Non-User
    profiles_users = zeros(total_users,1);
    profiles_users(internet_users_indices) = distributeValues([1:4],profiles_ratios,sum(internet_users_indices));
end


% ------------------------------------------------------------------------
%  User Activity Models for each user type and service
% ------------------------------------------------------------------------

% Concurrency Rates Calculation
% alpha = concurrency*(24*60)/((1-nonusers_ratio)*sum(profiles_connection_ratios.*profiles_connection_duration.*profiles_ratios));
% profiles_concurrency_rates = alpha/(24*60)*(profiles_connection_ratios.*profiles_connection_duration);
alphas = profiles_connection_ratios.*profiles_connection_duration.*profiles_ratios;
alphas_norm = alphas/sum(alphas);
profiles_concurrency_rates = alphas_norm*concurrency./(profiles_ratios*(1-nonusers_ratio));

% User Activity --> 1:ON, 0:OFF
if (RANDOMIZE)
    active_users = zeros(total_users,1);
    for i=1:length(profiles_concurrency_rates),
        aux_profile_indices = profiles_users(:,1)== i;
        active_users(aux_profile_indices) = randsample( [1;0], sum(aux_profile_indices==1), true, [profiles_concurrency_rates(i) 1-profiles_concurrency_rates(i)] );    
    end
    clear i aux_profile_indices;
else
    active_users = zeros(total_users,1);
    for i=1:length(profiles_concurrency_rates),
        aux_profile_indices = profiles_users(:,1)== i;
        if (profiles_concurrency_rates(i)<=1)
            active_users(aux_profile_indices) = distributeValues( [1 0], [profiles_concurrency_rates(i) 1-profiles_concurrency_rates(i)] ,sum(aux_profile_indices==1) );
        else
            active_users(aux_profile_indices) = 1;
        end
    end
    clear i aux_profile_indices;
end



%User Services Activity 1:ON, 0:OFF
aux_active_indices = active_users(:,1)== 1;
active_user_services = zeros(total_users,length(profiles_ratios));

if (RANDOMIZE)
    for i=1:length(profiles_ratios),
        %Indices
        %aux_profile_indices = profiles_users(:,1)== i;
        aux_active_profile_indices = (profiles_users(:,1)==i) & (aux_active_indices);

        %sum(aux_active_profile_indices)

        aux_bg_serv_indices = find(services_concurrency_types==BG_SERVICE_ID);
        aux_fg_serv_indices = find(services_concurrency_types==FG_SERVICE_ID);

        if (sum(aux_active_profile_indices)>0)
            %Background services
            for j=1:size(aux_bg_serv_indices),
                active_user_services(aux_active_profile_indices,j) = randsample( [1;0], sum(aux_active_profile_indices), true, [profiles_services_matrix(i,j) 1-profiles_services_matrix(i,j)] );
            end
            clear j;

            %Foreground services
            aux_services_being_used_indices= randsample(aux_fg_serv_indices', sum(aux_active_profile_indices), true, profiles_services_matrix(i,aux_fg_serv_indices));
            for j=1:length(aux_fg_serv_indices),
                active_user_services(aux_active_profile_indices,aux_fg_serv_indices(j))=aux_services_being_used_indices==aux_fg_serv_indices(j);        
            end;
            clear j;
        end;

    end;
    clear i;
else
    for i=1:length(profiles_ratios),
        %Indices
        %aux_profile_indices = profiles_users(:,1)== i;
        aux_active_profile_indices = (profiles_users(:,1)==i) & (aux_active_indices);

        %sum(aux_active_profile_indices)

        aux_bg_serv_indices = find(services_concurrency_types==BG_SERVICE_ID);
        aux_fg_serv_indices = find(services_concurrency_types==FG_SERVICE_ID);

        if (sum(aux_active_profile_indices)>0)
            %Background services
            for j=1:size(aux_bg_serv_indices),
                active_user_services(aux_active_profile_indices,j) = distributeValues( [1 0], [profiles_services_matrix(i,j) 1-profiles_services_matrix(i,j)] , sum(aux_active_profile_indices) );
            end
            clear j;

            %Foreground services
            aux_services_being_used_indices = distributeValues( aux_fg_serv_indices, profiles_services_matrix(i,aux_fg_serv_indices) , sum(aux_active_profile_indices));
            %aux_services_being_used_indices = randsample(aux_fg_serv_indices', sum(aux_active_profile_indices), true, profiles_services_matrix(i,aux_fg_serv_indices))
            for j=1:length(aux_fg_serv_indices),
                active_user_services(aux_active_profile_indices,aux_fg_serv_indices(j))=aux_services_being_used_indices==aux_fg_serv_indices(j);        
            end;
            clear j;
        end;

    end;
    clear i;
end


%% RESULTs

data = [user_household_id internet_users profiles_users active_users active_user_services];
simScenario = active_user_services(aux_active_indices,:);
bandwidths = [user_household_id(aux_active_indices) users_bandwidth(aux_active_indices)];


%% TRACES

% Results traces (2)
if (DEBUG_LEVEL>=1)
   fprintf('\tActive Users: %d/%d (%.2f%%)\n',sum(active_users),total_users,sum(active_users)/total_users*100); 
end

% Low-level debugging traces (2)
if (DEBUG_LEVEL>=2)
    fprintf('\t\tNon-Internet Users (Type=0): %d \n', sum(profiles_users==0));
    for i=1:length(profiles_ratios),
        fprintf('\t\tActive Type %d Users (%s): %d of %d (%.1f%% of actives) \n',i, profiles_names{i}, sum(profiles_users(aux_active_indices)==i), sum(profiles_users==i),sum(profiles_users(aux_active_indices)==i)/sum(active_users)*100);
    end
    clear i;
    
    for i=1:length(services_concurrency_types),
        fprintf('\t\tUsers using service %d (%s): %d (%.1f%% of actives) \n', i, services_names{i}, sum(active_user_services(:,i)), sum(active_user_services(:,i))/sum(active_users)*100);
    end
    
end

% High-level debugging traces (3)
if (DEBUG_LEVEL>=3)
    columns={'household_id','internet_user' 'user_profile' 'active_users' 'File Sharing' 'Video' 'Web and others' 'Online Gaming'};
    s1 = sprintf('%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t\n', columns{:});
    s2 = sprintf('%12d\t%12d\t%12d\t%12d\t%12d\t%12d\t%12d\t%12d\n', data');
    fprintf('\n\t\tScenario Details:\n\n')
    disp(horzcat(s1,s2))
    clear debug_mat columns s1 s2 debug_results;
end

%Remove path
rmpath('utils');

end
