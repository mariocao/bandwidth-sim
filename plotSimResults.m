function [ figures ] = plotSimResults( handles )
%PLOTSIMRESULTS Plot simulation results
%   Shows 3 graphs with the simulation results


%% FIGURE 1: Scenario Details
figure(1);
orient tall
% Data
userbwData = histc( handles.users_bandwidth(:,2), unique(handles.users_bandwidth(:,2)) );
usertypesData = histc(handles.scenarioData(handles.scenarioData(:,4)==1,3),1:max(handles.scenarioData(:,3)));
userServicesData = zeros(size(handles.simScenario,2),1);
for i=1:size(handles.simScenario,2),
    userServicesData(i)=sum(handles.simScenario(:,i));
end

%Bandwidths
subplot(3,1,1);
graph1=bar(1:length(userbwData),userbwData,'LineWidth',1);
hTitle1 = title('User capacities (Mbps)');
hYLabel1 = ylabel('Users (-)');
hXLabel1 = xlabel('Capacity (Mbps)');
set(gca,'XTickLabel', unique(handles.users_bandwidth(:,2))/1e6);
set(gca, ...
  'FontName'   , 'Helvetica', ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );


%User Types
subplot(3,1,2);
graph2=bar(1:max(handles.scenarioData(:,3)),usertypesData);
hTitle2 = title('Active user profiles');
hYLabel2 = ylabel('Users (-)');
hXLabel2 = xlabel('User profile (-)');
set(gca, ...
  'FontName'   , 'Helvetica', ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );


%User Types
subplot(3,1,3);
graph3=bar(1:size(handles.simScenario,2),userServicesData);
hTitle3 = title('Application types being used');
hYLabel3 = ylabel('Users (-)');
hXLabel3 = xlabel('Application type (-)');
set(gca, ...
  'FontName'   , 'Helvetica', ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );

set([ ...
    hXLabel1, hXLabel2, hXLabel3, ...
    hYLabel1, hYLabel2, hYLabel3], ...
    'FontName'   , 'Helvetica', ...
    'FontSize'   ,  9);
set( [hTitle1, hTitle2, hTitle3], ...
    'FontName'   , 'Helvetica', ...
    'FontSize'   , 11);
set([graph1, graph2, graph3], ...
  'LineWidth'   , 1    );
orient tall


%% FIGURE 2: SIMULTANEOUS USERs and SERVICEs
figure(2);

%Simultaneous Users
subplot(2,1,1);

graph1=bar(handles.simult_users_data,handles.simult_users_percent*100);
hTitle  = title('Percent of simultaneous active users');
hXLabel = xlabel('Users (-)');
hYLabel = ylabel('Percent (%)');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica', ...
    'FontSize'   ,  9);
set( hTitle                    , ...
    'FontSize'   , 11);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );


%Simultaneous Services
subplot(2,1,2);

graph2=bar(handles.simult_services_data,handles.simult_services_percent*100);
hTitle  = title('Percent of simultaneous applications');
hXLabel = xlabel('Applications (-)');
hYLabel = ylabel('Percent (%)');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica', ...
    'FontSize'   ,  9);
set( hTitle                    , ...
    'FontSize'   , 11);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'in'     , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );

set([graph1, graph2],  ...
  'LineWidth'   , 1 , ...
  'BarWidth'    , 1 );


%% FIGURE 3: Grade of Service (GoS)
figure(3);

%GoS

graph = stairs(handles.aggbw_caps/1e6,handles.aggbw_gos*100);

hTitle  = title('Grade of Service (GoS)');
hXLabel = xlabel('Bandwidth (Mbps)');
hYLabel = ylabel('GoS (%)');
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'Helvetica', ...
    'FontSize'   ,  9);
set( hTitle                    , ...
    'FontSize'   , 11);
set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'XGrid'       , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'LineWidth'   , 1         );
set(graph, ...
  'LineWidth'   , 1    );


end

