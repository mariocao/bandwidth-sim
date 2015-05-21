function [ map_aggr_bw, map_simult_users , map_simult_services, activity_factor ...
    ] = bandwidthSim( sim_scenario, household_bandwidth, sim_until, transient_period, stationary_start, debug_level, gui_on, tolerance )
%BANDWITHSIM Simulation of bandwidth allocation for a network scenario.
%   ***Simulates the bandwidth allocation for a network scenario defined by
%   a set of users and a set of user services.
%
% USAGE:
%   [timesteps, res_Ceq_users, res_Ceq_per_service]=bandwidthSim(C,data,period,[DEBUG,tolerance])
%
% INPUT:
%   C           - shared capacity
%   data        - m-by-n matrix of scenario data (m-users by n-services)
%   period      - simulation period
%   DEBUG       - debug level for traceability
%   transient   - transient period to be considered
%   tolerance   - mathematical tolerance for capacity operations
%
% OUTPUT:
%   timesteps      - simulation time steps
%   sample_Ceq_user1  - sampled capacity available for each user (ON+OFF)
%   sample_Ceq_services1 - sampled capacity available for each active user service (ON+OFF)
%   sample_Ceq_user2  - capacity available for each users for each step (ON)
%   sample_Ceq_services2 - capacity available for each active user service (ON)
%
% EXAMPLES:
%   tbc
%
% HISTORY:
% version 0.1: Release 17-Jul-2014
% version 0.2: Release 14-Aug-2014
%
% See also: runSim

% Author: Mario Cao Cueto
% DIT-UPM, Universidad Politécnica de Madrid
% eMail: mcao@dit.upm.es
% Copyright 2014 Mario Cao Cueto.
% $First created: 07-Jul-2014
% $Revision: 0.2 $ $Date: 17-Aug-2014 $


%% INIT

DEBUG_LEVEL=0;
TOLERANCE = 1e-8;
GUI_ON=0;

if nargin < 5
    error('usage:  bandwidthSim( sim_scenario, household_bandwidth, sim_until, transient_period, stationary_start, [debug_level, gui_on, tolerance] )')
elseif nargin == 6
    DEBUG_LEVEL = debug_level;
elseif nargin == 7
    DEBUG_LEVEL = debug_level;
    GUI_ON = gui_on;
elseif nargin == 8
    DEBUG_LEVEL = debug_level;
    TOLERANCE = tolerance;
end

SERVICE_ON=1;
SERVICE_OFF=2;

addpath('models');
addpath('utils');

N_users = size(sim_scenario,1);
N_services = size(sim_scenario,2);
users_services_status = zeros(N_users,N_services);
users_remaining_bytes = zeros(N_users,N_services);
users_remaining_off_time = Inf(N_users,N_services);

% Counter for conditional Toff times
users_ton_start = Inf(N_users,N_services);


%% PRINT INITIAL DATA

printline = '------------------------------------------------------------------------\n';

if (DEBUG_LEVEL >0)
    fprintf(printline);
    fprintf('Bandwidth Allocation Simulation\n');
    fprintf('\tConcurrent users: %d\n', N_users);
    fprintf('\tConcurrent services: %d\n', sum(any(sim_scenario)));
    fprintf('\tTotal of concurrent user services: %d\n', sum(sum(sim_scenario)));
    fprintf(printline);
end


%% INITIALIZATION OF SERVICE USERS (t=0): all users starts using their active services

%TextProgressBar
if (DEBUG_LEVEL>=1 && DEBUG_LEVEL<3)
    textprogressbar('Initializing services:    ');
end

if (GUI_ON) 
    h = waitbar(0,'Initializing services...','Name','Bandwitdh Allocation Simulation');
end

%For each user
for user_ind=1:N_users,
    %For each user service
    for service_ind=1:N_services,
        if (sim_scenario(user_ind,service_ind))
            users_services_status(user_ind,service_ind)=SERVICE_OFF;
            switch service_ind
                case 1 %P2P
                    users_remaining_off_time(user_ind,service_ind)=init_time(transient_period);
                case 2 %VIDEO
                    users_remaining_off_time(user_ind,service_ind)=init_time(transient_period);
                case 3 %WEB
                    users_remaining_off_time(user_ind,service_ind)=init_time(transient_period);
                case 4 %GAMING
                    users_remaining_off_time(user_ind,service_ind)=init_time(transient_period);
                otherwise
                    error('Error: unknown service');
            end
            
            if (DEBUG_LEVEL>2)
                fprintf('(INIT) user %d init time of %.4f msec for service %d \n',user_ind,users_remaining_off_time(user_ind,service_ind)*1e3,service_ind);
            end
            
        end
    end
    if (DEBUG_LEVEL>=1 && DEBUG_LEVEL<3)
        textprogressbar(round(user_ind/N_users*100));
    end
    if (GUI_ON) 
        waitbar(user_ind/N_users,h);
    end
end
clear user_ind service_ind;

if (DEBUG_LEVEL>=1 && DEBUG_LEVEL<3)
    textprogressbar(' done');
end
if (GUI_ON) 
    waitbar(1,h);
end

%% FIRST ITERATION of BANDWIDTH ALLOCATION (valid for both ON and OFF init)

% %Data for given t
% current_cap_users = zeros(N_users,1);
% users_ON_mask = any(users_remaining_bytes,2); %sum(users_remaining_bytes')'~=0
% current_cap_users(users_ON_mask) = household_bandwidth(users_ON_mask,2);
% clear users_ON_mask;

%Data for given t: NEW version with household capacity limitation
users_ON_mask = any(users_remaining_bytes,2);
current_cap_users = zeros(N_users,1);
user_sharefactor = arrayfun(@(t)nnz(household_bandwidth(users_ON_mask,1)==t), household_bandwidth(users_ON_mask,1));
current_cap_users(users_ON_mask) = household_bandwidth(users_ON_mask,2)./user_sharefactor;
clear users_ON_mask;

%Compute list of events (timestamps) and identify next event
current_cap_services = zeros(N_users, N_services);
user_on_stop_time = Inf(N_users,N_services);
for user_ind=1:N_users,
    aux_active_services=users_remaining_bytes(user_ind,:)~=0;
    aux_cap_per_service=current_cap_users(user_ind)/sum(aux_active_services);
    if (isnan(aux_cap_per_service))
        aux_cap_per_service=0;
    end
    current_cap_services(user_ind,:)=sim_scenario(user_ind,:)*aux_cap_per_service;
    user_on_stop_time(user_ind,aux_active_services)=users_remaining_bytes(user_ind,aux_active_services)*8/aux_cap_per_service;
end
clear user_ind;

user_event_times = min(user_on_stop_time,users_remaining_off_time);


%% SIMULATION LOOP until SIMULATION TIME IS EXCEEDED

%TextProgressBar
if (DEBUG_LEVEL>=1 && DEBUG_LEVEL<3)
    textprogressbar('Simulating ON/OFF events: ');
end
if (GUI_ON) 
    waitbar(0,h,'Simulating ON/OFF events...');
end

elapsed_time=0;
steps=1;

% %SAMPLED RESULTS: auxiliary variables for loop iterations
% aux_timesteps=[];
% aux_res_cap_users=[];
% aux_res_cap_users(1:size(current_cap_users,1),1)=current_cap_users;
% aux_res_cap_services=[];
% aux_res_cap_services(1:size(current_cap_services,1),1:size(current_cap_services,2))=current_cap_services;
% 
% %SAMPLED RESULTS: final results with memory pre-allocation
% sampleStep=0;
% sample_timestamps=zeros(sample_size,1);
% sample_cap_user1=zeros(N_users,sample_size);
% sample_cap_user2=zeros(N_users,sample_size);
% sample_cap_services1=zeros(N_users,N_services,sample_size);
% sample_cap_services2=zeros(N_users,N_services,sample_size);
% sample_user_alphas=zeros(N_users,sample_size);
% 
% %SAMPLE EXTRA RESULTs
% sample_simult_users=zeros(sample_size,1);
% sample_simult_services=zeros(sample_size,1);
% 
% %MAIN RESULT
% sample_total_bandwidth=zeros(sample_size,1);

%INFO. EXTRACTION MODULE params
map_aggr_bw = containers.Map('KeyType','double','ValueType','double');
map_simult_users = containers.Map('KeyType','double','ValueType','double');
map_simult_services = containers.Map('KeyType','double','ValueType','double');

aggr_activity_time = 0;
aggr_inactivity_time = 0;
first_already_saved = 0;


%MAIN LOOP with sampling capabilities
while ( true )
        
    %Identify next event
    clear rows columns;
    [rows, columns]=find(user_event_times==min(user_event_times(:)));
    %Avoid multiple equal next timestamps
    next_timestep = user_event_times(rows(1),columns(1));
    
%     % -------------------------
%     %  RESULTS SAMPLING MODULE
%     % -------------------------
%     if (DEBUG_LEVEL>=3)
%         fprintf('\tSampling monitor: Next simulation step %d at %.2f msec (+%.4f msec)\n',steps,(elapsed_time+next_timestep)*1e3,next_timestep*1e3);
%     end
%     
%     %Condition: next timestep is greater than next sample step
%     if (elapsed_time+next_timestep>=(sim_until/sample_size)*(sampleStep+1))
%         
%         aux_loops=0;
% 
%         %Loop until sample steps within the next timestep are considered
%         while (elapsed_time+next_timestep>=(sim_until/sample_size)*(sampleStep+1) && ...
%             (sim_until/sample_size)*(sampleStep+1)<=sim_until )
%         
%             sampleStep=sampleStep+1;
%            
%             %Detect simple or multiple sample steps
%             if ( (length(aux_timesteps)==1) && ( (aux_timesteps>=sim_until/sample_size) || (abs(aux_timesteps-sim_until/sample_size)<=TOLERANCE)  ) )
%                 aux_timesteps = sim_until/sample_size;
%             else
%                 aux_timesteps(end+1)=(sim_until/sample_size)*(sampleStep)-elapsed_time-aux_loops*(sim_until/sample_size);
%             end
%             
%             %Print debugging traces
%             if (DEBUG_LEVEL>=3)
%                 fprintf('\t--> Sample step %d reached (%.2f)\n',sampleStep,(sim_until/sample_size)*(sampleStep)*1e3);
%                 fprintf('\t\t Sum of timesteps within this SAMPLE STEP: %.4f\n',sum(aux_timesteps*1e3));
%                 fprintf('\t\t Timesteps: %.4f', aux_timesteps(1)*1e3);
%                 fprintf(' - %.4f',aux_timesteps(2:end)*1e3);
%                 fprintf('\n');
%             end
%            
%             %Users Capacities
%             [r11,r12]=calculateMeans(aux_res_cap_users,aux_timesteps);
%             sample_cap_user1(:,sampleStep)=r11;
%             sample_cap_user2(:,sampleStep)=r12;
%             
%             %User Services Capacities 
%             for i=1:N_services,
%                 service_caps_active=squeeze(aux_res_cap_services(sim_scenario(:,i)==1,i,:));
%                 %service_caps = squeeze(aux_res_cap_services(:,i,:));
%                 if (not(isempty(service_caps_active)))
%                     [r21,r22]=calculateMeans(service_caps_active,aux_timesteps);
%                 else
%                     r21=0;
%                     r22=0;
%                 end
%                 sample_cap_services1(sim_scenario(:,i)==1,i,sampleStep)=r21;
%                 sample_cap_services2(sim_scenario(:,i)==1,i,sampleStep)=r22;
%             end
%                       
%             %Alphas: Activity factors
%             sample_user_alphas(:,sampleStep) = sum(bsxfun(@times,aux_res_cap_users~=0,aux_timesteps/sum(aux_timesteps)),2);
%             
%             %Timestamps
%             sample_timestamps(sampleStep)=(sim_until/sample_size)*(sampleStep);
%             
%             %Extra results
%             %timesteps_normalized = timesteps/sum(timesteps);
%             sample_simult_users(sampleStep) = sum(sum(aux_res_cap_users~=0,1).*aux_timesteps/sum(aux_timesteps));
%             sample_simult_services(sampleStep) = sum(squeeze(sum(sum(aux_res_cap_services~=0,1),2))'.*aux_timesteps/sum(aux_timesteps));
%             
%             %MAIN RESULT: Total Demanded Capacity/Bandwitdth
%             sample_total_bandwidth(sampleStep) = sum(sum(aux_res_cap_users,1).*aux_timesteps/sum(aux_timesteps));
%             
%             %MARIO - CAMBIAR
% %             aggr_band = sum(aux_res_cap_users);
% %             if ( test_aggr_bandwidth.isKey(aggr_band) ) 
% %                 
% %             else
% %                 test_aggr_bandwidth(aggr_band) = 
% %             end
%             
% %             aux_res_cap_users
% %             aux_timesteps
%             
% %             test_aggr_bandwidth = [ test_aggr_bandwidth sum(aux_res_cap_users) ];
% %             test_aux_timesteps = [ test_aux_timesteps aux_timesteps ];
%             
%             %UPDATE variables for next loop run
%             aux_timesteps=elapsed_time+next_timestep-(sim_until/sample_size*sampleStep);
%             
%             aux_res_cap_users=aux_res_cap_users(:,end);
%             aux_res_cap_services=aux_res_cap_services(:,:,end);
%             aux_loops=aux_loops+1;
%             
%             if (DEBUG_LEVEL>=3)
%                 fprintf('\t\t Portion of Timestep(s) left: %.4f\n',aux_timesteps*1e3);
%             end
%         end
%         
%     %No sample step reached... saving next_timestep
%     else 
%         aux_timesteps(end+1) = next_timestep;
%     end
%     
%     % -------------------------
%     %  END SAMPLING MODULE
%     % -------------------------
    
    % --------------------------------
    %  START INFO. EXTRACTION MODULE
    % --------------------------------
    
    %Aggregated Bandwidth
    aggr_bw = sum(current_cap_users);
    %Simult. Users (Active & Inactives)
    aggr_actives = sum(current_cap_users~=0);
    aggr_inactives = sum(current_cap_users==0);
    %Simult. Services
    simult_services = sum(sum(current_cap_services~=0,1),2);
    
    aggr_ts = 0;
    
    if (elapsed_time+next_timestep>=sim_until)
        if (first_already_saved)
            %Last timestep value: sim_until - elapsed_time
            %disp('last');
            aggr_ts = sim_until - elapsed_time;
            
            if ( map_aggr_bw.isKey(aggr_bw) )
                map_aggr_bw(aggr_bw) = map_aggr_bw(aggr_bw) + aggr_ts;
            else
                map_aggr_bw(aggr_bw) = aggr_ts;
            end
            
            if ( map_simult_users.isKey(aggr_actives) )
                map_simult_users(aggr_actives) = map_simult_users(aggr_actives) + aggr_ts;
            else
                map_simult_users(aggr_actives) = aggr_ts;
            end
            
            if ( map_simult_services.isKey(simult_services) )
                map_simult_services(simult_services) = map_simult_services(simult_services) + aggr_ts;
            else
                map_simult_services(simult_services) = aggr_ts;
            end
            
        else
            %Last timestep and FIRST value (Mostly infrequent): sim_until - elapsed_time
            %disp('last and first');
            aggr_ts = sim_until - stationary_start;
            
            map_aggr_bw(aggr_bw) = sim_until - stationary_start;
            map_simult_users(aggr_actives) = sim_until - stationary_start;
            map_simult_services(simult_services) = sim_until - stationary_start;
            
        end
        
        aggr_activity_time = aggr_activity_time + aggr_ts*aggr_actives;
        aggr_inactivity_time = aggr_inactivity_time + aggr_ts*aggr_inactives;
        
        break;
        
    elseif ( elapsed_time+next_timestep>=stationary_start )
        if  (first_already_saved)
            %Standard timestep value: next_timestep
            %disp('normal');
            aggr_ts = next_timestep;
            
            if ( map_aggr_bw.isKey(aggr_bw) )
                map_aggr_bw(aggr_bw)=map_aggr_bw(aggr_bw)+aggr_ts;
            else
                map_aggr_bw(aggr_bw) = aggr_ts;
            end
            
            if ( map_simult_users.isKey(aggr_actives) )
                map_simult_users(aggr_actives) = map_simult_users(aggr_actives) + aggr_ts;
            else
                map_simult_users(aggr_actives) = aggr_ts;
            end
            
            if ( map_simult_services.isKey(simult_services) )
                map_simult_services(simult_services) = map_simult_services(simult_services) + aggr_ts;
            else
                map_simult_services(simult_services) = aggr_ts;
            end
            
        else 
            %First timestep value: elapsed + next_timestep - stationary_start
            %disp('first');
            aggr_ts = elapsed_time + next_timestep - stationary_start;
            
            map_aggr_bw(aggr_bw) = aggr_ts;
            map_simult_users(aggr_actives) = aggr_ts;
            map_simult_services(simult_services) = aggr_ts;
            first_already_saved=1;
            
        end
        
        aggr_activity_time = aggr_activity_time + aggr_ts*aggr_actives;
        aggr_inactivity_time = aggr_inactivity_time + aggr_ts*aggr_inactives;
        
    end
    
    % --------------------------------
    %  END INFO. EXTRACTION MODULE
    % --------------------------------
    
%     %Finish until Simulation Time Input  ---->>> IS BEING MODIFIED!!!!!!
%     if (elapsed_time+next_timestep>=sim_until)
%         %Next line commented preserve data for printing debugging traces
%         %elapsed_time=sim_until;
%         
%         %AQUI MARIO - añadir 
% %         current_cap_users
% %         next_timestep = sim_until - elapsed_time
%           
%         aggr_bw = sum(current_cap_users);
%         if ( aggr_bw_map.isKey(aggr_bw) )
%         	aggr_bw_map(aggr_bw)=aggr_bw_map(aggr_bw)+next_timestep;
%         else
%             aggr_bw_map(aggr_bw) = next_timestep;
%         end
%         
% %         aggr_bw_map.keys()
% %         sum(cell2mat(aggr_bw_map.values()))
%         
%         break;
%     else
%         %AQUI MARIO - añadir guardar capacidad agregada!!!
% %         current_cap_users
% %         next_timestep
%         
%         if(elapsed_time+next_timestep>=transient_period*2)
%             
%             
%         end
%         elapsed_time+next_timestep
%         next_timestep
%         
%         aggr_bw = sum(current_cap_users);
%         if ( aggr_bw_map.isKey(aggr_bw) )
%            aggr_bw_map(aggr_bw)=aggr_bw_map(aggr_bw)+next_timestep;
%         else
%             aggr_bw_map(aggr_bw) = next_timestep;
%         end
%         
%     end
    
    %Update variables of simulation
    elapsed_time = elapsed_time+next_timestep;
        
    %TextProgressBar
    if (DEBUG_LEVEL>=1 && DEBUG_LEVEL<3)
        textprogressbar(round(elapsed_time/sim_until*100));
    end
    if (GUI_ON) 
        waitbar(elapsed_time/sim_until,h);
    end
        
    %Update bytes consumed
    users_remaining_bytes=users_remaining_bytes-current_cap_services*next_timestep/8;
    %Update remaining off time
    users_remaining_off_time=users_remaining_off_time-next_timestep;
    
    %Detect triggered event(s) (it considers multiple event triggers!)
    for i=1:length(rows),
            user_ind = rows(i);
            service_ind =columns(i);
            
            %DETECTION OF NEW OFF EVENT(s)
            if (users_services_status(user_ind,service_ind)==SERVICE_ON)
                
                if (abs(users_remaining_bytes(user_ind,service_ind))<=TOLERANCE)
                    users_remaining_bytes(user_ind,service_ind)=0;
                elseif (abs(users_remaining_bytes(user_ind,service_ind))>=TOLERANCE)
                    error('ERROR detecting new OFF event (remaining_bytes:%d) \nPlease try to change the tolerance value.\n)',users_remaining_bytes(user_ind,service_ind));
                end
                    
                if ( users_remaining_off_time(user_ind,service_ind)~=Inf(1))
                   error('ERROR detecting new OFF event (remaining_off_time:%d) \n)',users_remaining_off_time(user_ind,service_ind)); 
                end

                if (DEBUG_LEVEL>=3)
                    fprintf('(OFF) user %d stops service %d after %.4f msec \t\t(step %d, elapsed: %.4f sec)\n',user_ind,service_ind,next_timestep*1e3,steps,elapsed_time);
                end
                
                %New status: OFF
                users_services_status(user_ind,service_ind)=SERVICE_OFF;
                
                %v25
                
                switch service_ind
                    case 1
                        users_remaining_off_time(user_ind,service_ind)=filesharing_off();
                    case 2
                        users_remaining_off_time(user_ind,service_ind)=video_off(elapsed_time-users_ton_start(user_ind,service_ind));
                    case 3
                        users_remaining_off_time(user_ind,service_ind)=web_off();
                    case 4
                        users_remaining_off_time(user_ind,service_ind)=gaming_off(elapsed_time-users_ton_start(user_ind,service_ind));
                    otherwise
                        error('Error: unknown service');
                end
                
                %v25
                users_ton_start(user_ind,service_ind)=Inf(1);                
                
                if (DEBUG_LEVEL>=4)
                    fprintf('\t -> OFF Time: %.2f msec\n', users_remaining_off_time(user_ind,service_ind)*1e3);
                end
                
            %DETECTION OF NEW ON EVENT
            elseif (users_services_status(user_ind,service_ind)==SERVICE_OFF)
           
                if (users_remaining_bytes(user_ind,service_ind)~=0 || users_remaining_off_time(user_ind,service_ind)~=0)
                   error('ERROR detecting new ON event (remaining_bytes:%d, remaining_off_time:%.2f\n)',users_remaining_bytes(user_ind,service_ind),users_remaining_off_time(user_ind,service_ind)); 
                end
                
                %Now OFF time is infinite because is ON again
                users_remaining_off_time(user_ind,service_ind)=Inf(1);
                                
                if (DEBUG_LEVEL>=3)
                    fprintf('(ON)  user %d starts service %d after %.4f msec \t(step %d, elapsed: %.4f sec)\n',user_ind,service_ind,next_timestep*1e3,steps,elapsed_time);
                end
                
                %New status: ON
                users_services_status(user_ind,service_ind)=SERVICE_ON;
                
                %v25
                users_ton_start(user_ind,service_ind)=elapsed_time;
                
                switch service_ind
                    case 1
                        users_remaining_bytes(user_ind,service_ind)=filesharing_on();
                    case 2
                        users_remaining_bytes(user_ind,service_ind)=video_on();
                    case 3
                        users_remaining_bytes(user_ind,service_ind)=web_on();
                    case 4
                        users_remaining_bytes(user_ind,service_ind)=gaming_on();
                    otherwise
                        error('Error: unknown service');
                end
                
                if (DEBUG_LEVEL>=4)
                    fprintf('\t -> Remaining_bytes: %.2f KB\n', users_remaining_bytes(user_ind,service_ind)/1e3);
                end
                
            %DETECTION OF INACTIVE TRIGGERED EVENT
            else 
                error('ERROR: event triggered for inactive service');
            end
    end
    clear i user_ind service_ind;
    
%     %Data for given t
%     current_cap_users = zeros(N_users,1);
%     users_ON_mask = any(users_remaining_bytes,2);
%     current_cap_users(users_ON_mask) = household_bandwidth(users_ON_mask,2);
%     clear users_ON_mask;
    
    %Data for given t: NEW version with household capacity limitation
    users_ON_mask = any(users_remaining_bytes,2);
    current_cap_users = zeros(N_users,1);
    user_sharefactor = arrayfun(@(t)nnz(household_bandwidth(users_ON_mask,1)==t), household_bandwidth(users_ON_mask,1));
    current_cap_users(users_ON_mask) = household_bandwidth(users_ON_mask,2)./user_sharefactor;
    clear users_ON_mask;
    
    %Compute list of events (timestamps) and identify next event
    current_cap_services = zeros(N_users, N_services);
    user_on_stop_time = Inf(N_users,N_services);

    %First: ON services
    for user_ind=1:N_users,
        if (any(users_remaining_bytes(user_ind,:)))
            aux_active_services=users_remaining_bytes(user_ind,:)~=0;
            aux_cap_per_service=current_cap_users(user_ind)/sum(aux_active_services);
            current_cap_services(user_ind,:)=aux_cap_per_service*aux_active_services;
            user_on_stop_time(user_ind,aux_active_services)=users_remaining_bytes(user_ind,aux_active_services)*8/aux_cap_per_service;
        end
    end
    clear user_ind;
    
    user_event_times = min(user_on_stop_time,users_remaining_off_time);

    %Results recollection
%     res_Ceq_users(1:size(Ceq_users,1),end+1) = Ceq_users;
%     res_Ceq_per_service(1:size(Ceq_users_service,1),1:size(Ceq_users_service,2),end+1) = Ceq_users_service;
    %mean for one user and service: mean(results_Ceq_per_service(i,j,:))
    
%     %NEW SAMPLED RESULTS
%     aux_res_cap_users(1:size(current_cap_users,1),end+1)=current_cap_users;
%     aux_res_cap_services(1:size(current_cap_services,1),1:size(current_cap_services,2),end+1) = current_cap_services;
    
    steps=steps+1;
    
end


%% TRACES

%TextProgressBar (1)
if (DEBUG_LEVEL>=1 && DEBUG_LEVEL<3)
    textprogressbar(100);
    textprogressbar(' done');
end
if (GUI_ON)
    waitbar(1,h);
    delete(h); 
end

%High-Level debugging traces (3)
if (DEBUG_LEVEL >=4)
    fprintf('(---) next event trigger of %.2f msec (at %.2f sec) exceeds sim time (%.2f sec)\n',next_timestep*1e3,elapsed_time+next_timestep,sim_until);
    fprintf('(END) simulation time exceeded after %.2f msec \t(step %d, elapsed: %.2f sec)\n',(sim_until-elapsed_time)*1e3,steps,sim_until);
    fprintf(printline);
end

% if (DEBUG_LEVEL >=5)
%     sample_cap_user1
%     sample_cap_user2
% end
% 
% if (DEBUG_LEVEL >=5)
%     sample_cap_services1
%     sample_cap_services2
% end


% %Results (1)
% if (DEBUG_LEVEL >= 1)
%     fprintf('Simulation of %.2f sec finished after %d steps (in %d samples) \n',sim_until,steps,sample_size);
% end

%NEW Results (1)
if (DEBUG_LEVEL >= 1)
    fprintf('Simulation of %.2f sec finished after %d steps \n',sim_until,steps);
end


% r_capacities = cell2mat(map_aggr_bw.keys());
% r_timesteps = cell2mat(map_aggr_bw.values());


% [aggbw_results,aggbw_caps] = analyzeAggrBw(r_capacities,r_timesteps);
% save('test.mat','r_capacities','r_timesteps','aggbw_results');
% 
% %OUTPUT (cambiar mario)
% vertcat(aggbw_caps/1e6,aggbw_results*100)
% aggr_activity_time/(aggr_activity_time+aggr_inactivity_time)
% 
% calculateMeans(cell2mat(map_simult_users.keys()),cell2mat(map_simult_users.values()))


activity_factor = aggr_activity_time/(aggr_activity_time+aggr_inactivity_time);


rmpath('utils');
rmpath('models');

end






