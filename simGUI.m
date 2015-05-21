function varargout = simGUI(varargin)
%SIMGUI M-file for simGUI.fig
%      SIMGUI, by itself, creates a new SIMGUI or raises the existing
%      singleton*.
%
%      H = SIMGUI returns the handle to a new SIMGUI or the handle to
%      the existing singleton*.
%
%      SIMGUI('Property','Value',...) creates a new SIMGUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to simGUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SIMGUI('CALLBACK') and SIMGUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SIMGUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help simGUI

% Last Modified by GUIDE v2.5 07-Jan-2015 16:37:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @simGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before simGUI is made visible.
function simGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

loadState(handles);

handles.modelFilename = get(handles.model_edit,'String');
handles.households=str2double(get(handles.households_edit,'String'));
handles.concurrency = str2double(get(handles.concurrency_edit,'String'))/100;
handles.sampleSize=str2double(get(handles.samplesize_edit,'String'));
handles.sim_until=str2double(get(handles.simtime_edit,'String'));
handles.transient_period = str2double(get(handles.initialtime_edit,'String'));
handles.stationary_start = str2double(get(handles.stationary_edit,'String'));

handles.debuglevel = 2;
handles.lastDir = [];

handles.probabilities = 1:100;

title(handles.userbw_axes,'User capacities (Mbps)');
title(handles.usertypes_axes,'Active user profiles');
title(handles.userservices_axes,'Service types being used');
title(handles.gos_axes,'Grade of Service (GoS)');

% Choose default command line output for simGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes simGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = simGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in run_pushbutton.
function run_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to run_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic
set(handles.action_text,'String','Simulating...');
set(hObject,'Enable','off');

if (handles.debuglevel>=1)
    fprintf('Simulation Scenario Details\n');
    %fprintf('\tTotal Shared Bandwidth: %d Mbps\n',handles.C/1e6);
end

%UNCOMMENT
[handles.simScenario,handles.users_bandwidth, handles.scenarioData] = createScenario(handles.households, handles.concurrency, strcat(handles.modelFilename), handles.debuglevel);

% %SAMPLE SCENARIO with ONLY WEB users
% web_user = [0 0 1 0];
% handles.simScenario = bsxfun(@times,ones(handles.households,4),web_user);
% handles.users_bandwidth = [(1:handles.households)' ones(handles.households,1)*10e6];
% scenario_user= [ 1 4 1 0 0 1 0];
% handles.scenarioData = [(1:handles.households)' bsxfun(@times,ones(handles.households,7),scenario_user)];

% %SAMPLE SCENARIO FOR BUGs TRACEABILITY (1 video user)
% handles.simScenario = [1 1 0 0; 1 1 0 0];
% handles.users_bandwidth = [10e6; 10e6];
% handles.scenarioData = [1 4 1 0 1 1 0; 1 4 1 0 1 1 0];

% %SAMPLE SCENARIO FOR BUGs TRACEABILITY (2 users: video + p2p)
% handles.simScenario = [1     1     0     0; 1     0     1     0];
% handles.users_bandwidth = [40e6; 40e6];
% handles.scenarioData = [1     2     1     1 1 0 0; 1     3     1     1 0 1 0];

if (isempty(handles.simScenario))
    disp('No Active Users');
    set(hObject,'Enable','on');
    set(handles.action_text,'String','Run Sim failed: No active users.');
    return;
end

if (handles.debuglevel>=1)
    % Print scenario preparation results    
    %fprintf('\tBandwidth per Active User: %.1f Kbps\n',handles.C/size(handles.simScenario,1)/1e3);
end

% Show scenario results
guidata(hObject, handles);
showScenarioResults(handles);

%NEW version of bandwidthSim.m
 [ handles.map_aggr_bw, handles.map_simult_users , handles.map_simult_services, handles.activity_factor ] = ...
    bandwidthSim(handles.simScenario, handles.users_bandwidth, handles.sim_until, handles.transient_period, handles.stationary_start, handles.debuglevel, handles.debuglevel);

%to delete
% vertcat(cell2mat(map_aggr_bw.keys())/1e6,cell2mat(map_aggr_bw.values()))


addpath('utils');

[handles.aggbw_gos, handles.aggbw_caps] = analyzeAggrBw(cell2mat(handles.map_aggr_bw.keys()),cell2mat(handles.map_aggr_bw.values()));

[handles.simult_users_percent, handles.simult_users_data] = analyzeAggrData(cell2mat(handles.map_simult_users.keys()),cell2mat(handles.map_simult_users.values()));
[handles.simult_services_percent, handles.simult_services_data] = analyzeAggrData(cell2mat(handles.map_simult_services.keys()),cell2mat(handles.map_simult_services.values()));

handles.simult_users = calculateMeans(cell2mat(handles.map_simult_users.keys()),cell2mat(handles.map_simult_users.values()));
handles.simult_services = calculateMeans(cell2mat(handles.map_simult_services.keys()),cell2mat(handles.map_simult_services.values()));

rmpath('utils');


% [handles.timestamps, handles.simult_users, handles.simult_services, handles.total_bandwidth, ...
%     handles.cap_users1, handles.cap_user_services1, handles.cap_users2, handles.cap_user_services2, handles.user_alphas ] = ...
%     bandwidthSim(handles.simScenario, handles.users_bandwidth, handles.sim_until, handles.transient_period, handles.stationary_start,  handles.sampleSize, handles.debuglevel, handles.debuglevel);

% [handles.results_global, handles.results_users, handles.results_services] = ... 
%     analyzeSimSamples(handles.simScenario, handles.timestamps, ... 
%     handles.simult_users, handles.simult_services, ...
%     handles.cap_users1, handles.cap_user_services1, handles.cap_users2, handles.cap_user_services2, handles.debuglevel);

%handles.r_prob_capacity  = ...
%    analyzeProbabilities(handles.timestamps, handles.transient_period, handles.probabilities,handles.total_bandwidth, 1);

%[handles.results_alpha_users, handles.results_alpha_user_types] = analyzeAlphas(handles.user_alphas, handles.scenarioData,1);

handles.elapsed_time = toc;

% showSimResults(handles);
showSimResults(handles);

guidata(hObject, handles);


set(handles.action_text,'String','Simulation successfully finished.');
% set(handles.elapsed_text,'String',sprintf('%.2f seconds',handles.elapsed_time));
set(handles.saveresults_pushbutton,'Enable','on');
set(handles.plot_pushbutton,'Enable','on');
set(hObject,'Enable','on');


% --------------------
%  CALLBACK FUNCTIONs
% --------------------
function initialtime_edit_Callback(hObject, eventdata, handles)
% hObject    handle to initialtime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of initialtime_edit as text
%        str2double(get(hObject,'String')) returns contents of initialtime_edit as a double
handles.transient_period = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function initialtime_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initialtime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function samplesize_edit_Callback(hObject, eventdata, handles)
% hObject    handle to samplesize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of samplesize_edit as text
%        str2double(get(hObject,'String')) returns contents of samplesize_edit as a double
handles.sampleSize = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function samplesize_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samplesize_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function simtime_edit_Callback(hObject, eventdata, handles)
% hObject    handle to simtime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of simtime_edit as text
%        str2double(get(hObject,'String')) returns contents of simtime_edit as a double
handles.sim_until = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function simtime_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to simtime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function households_edit_Callback(hObject, eventdata, handles)
% hObject    handle to households_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of households_edit as text
%        str2double(get(hObject,'String')) returns contents of households_edit as a double
handles.households = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function households_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to households_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function model_edit_Callback(hObject, eventdata, handles)
% hObject    handle to model_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of model_edit as text
%        str2double(get(hObject,'String')) returns contents of model_edit as a double

% --- Executes during object creation, after setting all properties.
function model_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in browse_pushbutton.
function browse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to browse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile( ...
{  '*.cfg',  'Model files (*.cfg)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Choose a Model Config File');
handles.modelFilename = strcat(pathname,filename);
set(handles.model_edit,'String',handles.modelFilename);
guidata(hObject, handles);


% --- Executes on slider movement.
function concurrency_slider_Callback(hObject, eventdata, handles)
% hObject    handle to concurrency_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = round(get(hObject,'Value'));
set(hObject,'Value', sliderValue);
set(handles.concurrency_edit,'String', num2str(sliderValue));

handles.concurrency = sliderValue/100;
guidata(hObject, handles);


function concurrency_edit_Callback(hObject, eventdata, handles)
% hObject    handle to concurrency_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of concurrency_edit as text
%        str2double(get(hObject,'String')) returns contents of concurrency_edit as a double
textValue = str2double(get(hObject,'String'));
if (textValue >100)
    textValue = 100;
elseif (textValue <0)
    textValue = 0;
end
set(hObject,'String',num2str(textValue));
set(handles.concurrency_slider,'Value',textValue);

handles.concurrency = textValue/100;
guidata(hObject, handles);


function stationary_edit_Callback(hObject, eventdata, handles)
% hObject    handle to stationary_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.stationary_start = str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
guidata(hObject,handles);
saveState(handles);

addpath('utils');
figs=getFigsNoGUI();
rmpath('utils');

for i=1:length(figs)
    if ishandle(figs(i))
        close(figs(i))
    end
end

delete(hObject);


% ----------------------
%  SAVE AND LOAD STATE
% ----------------------
function saveState(handles)

state.modelFilenameStr = get(handles.model_edit, 'String');
state.householdsVal = get(handles.households_edit, 'String');
state.concurrencyVal = get(handles.concurrency_edit, 'String');
state.samplesizeVal = get(handles.samplesize_edit, 'String');
state.simtimeVal = get(handles.simtime_edit, 'String');
state.initialtimeVal = get(handles.initialtime_edit, 'String');
state.stationaryStartVal = get(handles.stationary_edit, 'String');

save('state.mat','state');

function loadState(handles)
fileName = 'state.mat';
if exist([pwd filesep fileName],'file')
    load(fileName);
    set(handles.model_edit, 'String', state.modelFilenameStr);
    set(handles.households_edit, 'String', state.householdsVal);
    set(handles.concurrency_edit, 'String', state.concurrencyVal);
    set(handles.concurrency_slider,'Value', str2double(state.concurrencyVal));
    set(handles.samplesize_edit, 'String', state.samplesizeVal);
    set(handles.simtime_edit, 'String', state.simtimeVal);
    set(handles.initialtime_edit, 'String', state.initialtimeVal);
    set(handles.stationary_edit, 'String', state.stationaryStartVal);
    delete(fileName);
end


% -----------------------
%  SHOW in GUI FUNCTIONs
% -----------------------
function showSimParams(handles)
set(handles.model_edit,'String',handles.modelFilename);
set(handles.households_edit,'String',handles.households);
set(handles.concurrency_edit,'String',handles.concurrency*100);
set(handles.concurrency_slider,'Value',handles.concurrency*100);
set(handles.samplesize_edit,'String',handles.sampleSize);
set(handles.simtime_edit,'String',handles.sim_until);
set(handles.initialtime_edit,'String',handles.transient_period);
set(handles.stationary_edit,'String',handles.stationary_start);


function showScenarioResults(handles)
active_users = size(handles.simScenario,1);
mean_bandwidth_per_user = sum(handles.users_bandwidth(:,2))/active_users;
total_users = size(handles.scenarioData,1);

set(handles.activeusers_text,'String',sprintf('%d / %d (%.2f%%)',active_users,total_users,active_users/total_users*100));
set(handles.bandwidthperuser_text,'String',strcat(num2str(mean_bandwidth_per_user/1e6,'%.2f'),' Mbps'));

%hist(handles.usertypes_axes,handles.scenarioData(:,3),max(handles.scenarioData(:,3))+1);

userbwData = histc( handles.users_bandwidth(:,2), unique(handles.users_bandwidth(:,2)) );
usertypesData = histc(handles.scenarioData(handles.scenarioData(:,4)==1,3),1:max(handles.scenarioData(:,3)));
userServicesData = zeros(size(handles.simScenario,2),1);
for i=1:size(handles.simScenario,2),
    userServicesData(i)=sum(handles.simScenario(:,i));
end

bar(handles.userbw_axes,1:length(userbwData),userbwData);
axis(handles.userbw_axes,'tight');
set(handles.userbw_axes,'XTickLabel', unique(handles.users_bandwidth(:,2))/1e6);
bar(handles.usertypes_axes,1:max(handles.scenarioData(:,3)),usertypesData);
axis(handles.usertypes_axes,'tight');
bar(handles.userservices_axes,1:size(handles.simScenario,2),userServicesData);
axis(handles.userservices_axes,'tight');

title(handles.userbw_axes,'User capacities (Mbps)');
title(handles.usertypes_axes,'Active user profiles');
title(handles.userservices_axes,'Service types being used');

function showSimResults(handles)
set(handles.simultusers_text,'String',num2str(handles.simult_users,'%.2f'));
set(handles.simultservices_text,'String',num2str(handles.simult_services,'%.2f'));
set(handles.activityfactor_text,'String',num2str(handles.activity_factor,'%.2f'));

stairs(handles.gos_axes,handles.aggbw_caps/1e6,handles.aggbw_gos*100);
axis(handles.gos_axes,'tight');

title(handles.gos_axes,'Grade of Service (GoS)');
xlabel(handles.gos_axes,'BW (Mbps)');
ylabel(handles.gos_axes,'GoS (%)');
set(handles.gos_axes,'XGrid','on', 'YGrid','on');

set(handles.elapsed_text,'String',num2str(handles.elapsed_time));


% -----------------------
%  SAVE AND LOAD BUTTONS
% -----------------------

% --- Executes on button press in saveresults_pushbutton.
function saveresults_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveresults_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.action_text,'String','Saving simulation results...');
SAVE_PATH=strcat('saved-data','\',datestr(now,'yyyymmdd-HHMM'), ...
    ' -',32,get(handles.households_edit,'String'),'subs',32, ...
    get(handles.simtime_edit,'String'),'sec\');
handles.lastDir = SAVE_PATH;
guidata(hObject, handles);
if ~exist(SAVE_PATH, 'dir')
  mkdir(SAVE_PATH);
else
  disp('Warning: Simulation Directory already exists');
end

h = waitbar(0,'Saving simulation data..','Name','Bandwidth Allocation Simulation');
fprintf('-> Saving data in %s ...',SAVE_PATH);

modelFilename=handles.modelFilename;
households=handles.households;
concurrency=handles.concurrency;
sampleSize=handles.sampleSize;
sim_until=handles.sim_until;
transient_period=handles.transient_period;
stationary_start=handles.stationary_start;
save(strcat(SAVE_PATH,'sim-params'), 'modelFilename', 'households', 'concurrency', ...
    'sampleSize', 'sim_until', 'transient_period','stationary_start');
waitbar(0.2,h);

simScenario=handles.simScenario;
users_bandwidth=handles.users_bandwidth;
scenarioData=handles.scenarioData;
save(strcat(SAVE_PATH,'sim-scenario'),  'simScenario','users_bandwidth','scenarioData');
waitbar(0.4,h);

map_aggr_bw=handles.map_aggr_bw;
map_simult_users=handles.map_simult_users;
map_simult_services=handles.map_simult_services;
activity_factor=handles.activity_factor;
elapsed_time=handles.elapsed_time;
save(strcat(SAVE_PATH,'sim-data'), 'map_aggr_bw', 'map_simult_users', 'map_simult_services', ...
        'activity_factor', 'elapsed_time');
waitbar(0.6,h);

aggbw_gos=handles.aggbw_gos;
aggbw_caps=handles.aggbw_caps;
simult_users_percent=handles.simult_users_percent;
simult_users_data=handles.simult_users_data;
simult_services_percent=handles.simult_services_percent;
simult_services_data=handles.simult_services_data;
simult_users=handles.simult_users;
simult_services=handles.simult_services;
save(strcat(SAVE_PATH,'sim-results'), 'aggbw_gos', 'aggbw_caps', ...
    'simult_users_percent', 'simult_users_data', ...
    'simult_services_percent', 'simult_services_data', ...
    'simult_users', 'simult_services');
waitbar(0.8,h);

fprintf('done! \n');

waitbar(1,h);
delete(h);

set(handles.action_text,'String','Simulation results successfully saved.');

% --- Executes on button press in load_pushbutton.
function load_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.action_text,'String','Loading simulation data...');
dirName = uigetdir;
if (dirName~=0)
    fileNames = {'sim-params','sim-scenario','sim-data','sim-results'};
    error = 0;
    for i = 1:numel(fileNames)
        fullFileName=sprintf('%s\\%s.mat',dirName,fileNames{i});
        if exist(fullFileName, 'file')
            load(fullFileName);
        else
            error=1;
            break;
        end
    end

    if (error==0)
        %sim-params
        handles.modelFilename=modelFilename;
        handles.households = households;
        handles.concurrency = concurrency;
        handles.sampleSize = sampleSize;
        handles.sim_until = sim_until;
        handles.transient_period = transient_period;
        handles.stationary_start = stationary_start;

        %sim-scenario
        handles.simScenario = simScenario;
        handles.users_bandwidth = users_bandwidth;
        handles.scenarioData = scenarioData;

        %sim-data
        handles.map_aggr_bw = map_aggr_bw;
        handles.map_simult_users = map_simult_users;
        handles.map_simult_services = map_simult_services;
        handles.activity_factor = activity_factor;
        handles.elapsed_time=elapsed_time;

        %sim-results
        handles.aggbw_gos = aggbw_gos;
        handles.aggbw_caps = aggbw_caps;
        handles.simult_users_percent = simult_users_percent;
        handles.simult_users_data = simult_users_data;
        handles.simult_services_percent = simult_services_percent;
        handles.simult_services_data = simult_services_data;
        handles.simult_users = simult_users;
        handles.simult_services = simult_services;
                
        guidata(hObject, handles);
        
        showSimParams(handles);
        showScenarioResults(handles);
        showSimResults(handles);
        
        set(handles.plot_pushbutton,'Enable','on');
        set(handles.saveresults_pushbutton,'Enable','on');
%         set(handles.stationary_checkbox,'Enable','on');
        set(handles.action_text,'String','Simulation data successfully loaded.');
    else
        warndlg('Ups! There has been an error loading the simulation data.','Warning');
        set(handles.action_text,'String','Error while loading simulation data.');
    end
end


% ---------------------------
%  PLOT AND SAVE FIGS BUTTONs
% ----------------------------
% --- Executes on button press in plot_pushbutton.
function plot_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plot_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.action_text,'String','Plotting simulation results...');

% if (get(handles.stationary_checkbox,'Value')) %on
%     handles.stationary_start = str2double(get(handles.stationary_edit,'String'));
% else %off
%     handles.stationary_start = handles.transient_period;
%     set(handles.stationary_edit,'String',handles.transient_period);
% end
guidata(hObject, handles);

%Simultaneous Users and Services

plotSimResults(handles);

guidata(hObject, handles);

set(handles.action_text,'String','Simulation plotting finished.');
set(handles.savefigs_pushbutton,'Enable','on');


% --- Executes on button press in savefigs_pushbutton.
function savefigs_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to savefigs_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.action_text,'String','Saving figures...');
    SAVE_PATH = uigetdir(handles.lastDir);
    if (SAVE_PATH~=0)
        addpath('utils');
        figures = getFigsNoGUI();
        saveFigures(SAVE_PATH, figures, 1);
        rmpath('utils');
        set(handles.action_text,'String','Simulation figures successfully saved.');
    else
        set(handles.action_text,'String','Figures saving interrupted.');
    end
