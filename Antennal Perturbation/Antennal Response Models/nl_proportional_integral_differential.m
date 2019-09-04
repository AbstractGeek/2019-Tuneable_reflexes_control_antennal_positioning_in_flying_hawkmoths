function [nlgr,data] = nl_proportional_integral_differential(input, output, sam_freq)
% function [nlgr,data] = nl_proportional_integral_differential(input, output, sam_freq)
% 
% 
% Dinesh Natesan, 17th May 2016

%% Create iddata
data = iddata(output,input,1/sam_freq);
% Define names etc
data.InputName = {'Set-Point', 'Antennal Angle','ReleaseState'};
data.InputUnit =  {'deg', 'deg','boolean'};
data.OutputName = 'Antennal Angle';
data.OutputUnit = 'deg';
data.Tstart = 0;
data.TimeUnit = 's';


%% Model Definition
FileName = 'nl_proportional_integral_differential_model';
Order = [1 3 2];
Parameters = {ones(2,2),ones(2,1),ones(1,2),1}';
InitialStates = [output(1);output(2)-output(1)];
Ts = 0;
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'nl_proportional_integral_differential');

% Set Input and Output names
set(nlgr, 'InputName', {'Set-Point', 'Antennal Angle','ReleaseState'}, ...
    'InputUnit',  {'deg', 'deg','boolean'},            ...
    'OutputName', 'Antennal Angle', ...
    'OutputUnit', 'deg',                         ...
    'TimeUnit', 's');

% Set init and par names
nlgr = setinit(nlgr, 'Name', {'Antennal Angle','Change in Antennal Angle'});
nlgr = setinit(nlgr, 'Unit', {'deg','deg'});
nlgr = setpar(nlgr, 'Name', {'A','B','C','D'});
nlgr = setpar(nlgr, 'Unit', {'constant', 'constant', 'constant','constant'});


nlgr = setinit(nlgr, 'Fixed',false);
opt = nlgreyestOptions('Display', 'off','SearchMethod','lsqnonlin',...
    'SearchOption',optimset('lsqnonlin'));
opt.SearchOption.MaxIter = 100;

nlgr = nlgreyest(data, nlgr, opt);

end