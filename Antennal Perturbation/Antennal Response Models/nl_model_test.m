%% Get necessary data
inputDir_pref = '/home/dinesh/Dropbox/Ncbs/Moth Projects/Antennal Perturbation Experiments/Analysis/returnCurves/projAngles-filtered';
treatment = 'Control';
sam_freq = 1000;
dataset = load(fullfile(inputDir_pref,treatment));
% ws_fitdata = dataset.filtReturnCurves.moth3_30jun15.Trial_1.ws_0_mps;
% ws_data = dataset.returnCurves.moth3_30jun15.Trial_1.ws_0_mps;
ws_fitdata = dataset.filtReturnCurves.moth3_30jun15.Trial_2.ws_4_mps;
ws_data = dataset.returnCurves.moth3_30jun15.Trial_2.ws_4_mps;
step_number = 1;
startpoint = dataset.startPoints.moth3_30jun15.Trial_2.ws_4_mps.Mean_last(step_number);
output = ws_fitdata(:,step_number);
input = NaN(size(output,1),3);
input(:,1) = ws_data{step_number,1}(end,2);
input(:,2) = output;
input(:,3) = ones(size(output));
input(1:startpoint-1,3) = 0;
data = iddata(output,input,1/sam_freq);
clearvars -except input output data dataset sam_freq

% %% define the model
% FileName      = 'nl_integrator_model';
% Order = [1 3 1];
% Parameters = [1,1,1];
% InitialStates = output(1);
% Ts = 0;
% nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
%     'Name', 'nl_integrator');
% 
% % Set Input and Output names
% set(nlgr, 'InputName', {'SetPoint','Output','OnoFF'}, ...
%     'InputUnit', {'deg', 'deg', 'Logical' },            ...
%     'OutputName', 'Position', ...
%     'OutputUnit', 'deg',                         ...
%     'TimeUnit', 's');
% 
% % Set parameter values
% setpar(nlgr,'Name',{'A','B','C'});
% 
% nlgr = setinit(nlgr, 'Fixed',false);
% opt = nlgreyestOptions('Display', 'on','SearchMethod','lsqnonlin',...
%     'SearchOption',optimset('lsqnonlin'));
% 
% nlgr = nlgreyest(data, nlgr, opt);

%% Test models
% [nlgr,data] = nl_integrator(input, output, sam_freq);
% [nlgr,data] = nl_proportional(input, output, sam_freq);
% [nlgr,data] = nl_proportional_integral(input, output, sam_freq);
[nlgr,data] = nl_proportional_integral_differential(input, output, sam_freq);
compare(data,nlgr);