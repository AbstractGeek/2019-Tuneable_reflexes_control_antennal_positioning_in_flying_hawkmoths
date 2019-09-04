%% Get necessary data
inputDir_pref = '/home/dinesh/Dropbox/Ncbs/Moth Projects/Antennal Perturbation Experiments/Analysis/returnCurves/projAngles-filtered';
treatment = 'Control';
sam_freq = 1000;
% Get dataset
dataset = load(fullfile(inputDir_pref,treatment));
ws_data = dataset.returnCurves.moth3_30jun15.Trial_2.ws_4_mps;
step_number = 1;
% Extract data
inputPadding = 100;
returnPadding = 0;
curr_data = ws_data{step_number,1};
output = curr_data(inputPadding-returnPadding+1:end,1);
input = curr_data(inputPadding-returnPadding+1:end,2);  
data = iddata(output,input,1/sam_freq);
clearvars -except input output data dataset sam_freq
% 
% %% define the model
% FileName      = 'l_integrator_model';
% Order = [1 1 1];
% Parameters = {'Ki',1};
% Ts = 1/sam_freq;
% lgr = idgrey(FileName, Parameters, 'c',{},0,...
% 'InputName', 'SetPoint', ...
% 'InputUnit', 'deg',            ...
% 'OutputName', 'Position', ...
% 'OutputUnit', 'deg',                         ...
% 'TimeUnit', 's');
% 
% % lgr = setinit(lgr, 'Fixed',false);
% opt = greyestOptions('Display', 'on','SearchMethod','lsqnonlin',...
%     'SearchOption',optimset('lsqnonlin'),'InitialState','backcast');
% 
% lgr = greyest(data, lgr, opt);

%% Test models
% [lgr,data] = l_integrator(input, output, sam_freq);
% [lgr,data] = l_proportional(input, output, sam_freq);
% [lgr,data] = l_proportional_integral(input, output, sam_freq)
[lgr,data] = l_proportional_integral_differential(input, output, sam_freq)
compare(data,lgr);