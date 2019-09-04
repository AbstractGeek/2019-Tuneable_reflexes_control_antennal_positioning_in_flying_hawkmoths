function [lgr,data] = l_proportional_differential(input, output, sam_freq, InitialState, InitialParameters)
% function [nlgr,data] = l_proportional_integral(input, output, sam_freq)
% 
% 
% Dinesh Natesan, 18th May 2016

if nargin <= 4
   predict_flag = false; 
elseif nargin == 5
    predict_flag = true; 
end

%% Create iddata
data = iddata(output,input,1/sam_freq);
% Define names etc
data.InputName = 'Set-Point';
data.InputUnit =  'deg';
data.OutputName = 'Antennal-Angle';
data.OutputUnit = 'deg';
data.Tstart = 0;
data.TimeUnit = 's';


%% Model Definition
FileName = 'l_proportional_differential_model';

if (predict_flag)
    Parameters = {'Kp',InitialParameters(1);'Kd',InitialParameters(2)};
else
    Parameters = {'Kp',1;'Kd',1};   
end

lgr = idgrey(FileName, Parameters, 'c',{},0,...
'InputName', 'Set-Point', ...
'InputUnit', 'deg',            ...
'OutputName', 'Antennal-Angle', ...
'OutputUnit', 'deg',                         ...
'TimeUnit', 's');

if (predict_flag)
    lgr.Structure.Parameters(1).Free = false;
    lgr.Structure.Parameters(2).Free = false;
end

opt = greyestOptions('Display', 'off','SearchMethod','lsqnonlin',...
    'SearchOption',optimset('lsqnonlin'),'InitialState',InitialState);

lgr = greyest(data, lgr, opt);

end