function [lgr,data] = l_double_integrator(input, output, sam_freq, InitialState, Ki)
% function [nlgr,data] = nl_integrator(input, output, sam_freq)
% 
% 
% Dinesh Natesan, 17th May 2016

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
FileName = 'l_double_integrator_model';

if (predict_flag)
    Parameters = {'Ki', Ki};
else
    Parameters = {'Ki', 1};
end


lgr = idgrey(FileName, Parameters, 'c',{},0,...
'InputName', 'Set-Point', ...
'InputUnit', 'deg',            ...
'OutputName', 'Antennal-Angle', ...
'OutputUnit', 'deg',                         ...
'TimeUnit', 's');

if (predict_flag)
    lgr.Structure.Parameters(1).Free = false;
end

opt = greyestOptions('Display', 'off','SearchMethod','lsqnonlin',...
    'SearchOption',optimset('lsqnonlin'),'InitialState',InitialState);

lgr = greyest(data, lgr, opt);

end