function [nlgr,data] = nl_proportional(input, output, sam_freq)
% function [nlgr,data] = nl_proportional(input, output, sam_freq)
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
FileName = 'nl_proportional_model';
Order = [1 3 0];
Parameters = 1;
InitialStates = [];
Ts = 0;
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'nl_proportional');

% Set Input and Output names
set(nlgr, 'InputName', {'Set-Point', 'Antennal Angle','ReleaseState'}, ...
    'InputUnit',  {'deg', 'deg','boolean'},            ...
    'OutputName', 'Antennal Angle', ...
    'OutputUnit', 'deg',                         ...
    'TimeUnit', 's');

% Set init and par names
nlgr = setinit(nlgr, 'Name', 'Antennal Angle');
nlgr = setinit(nlgr, 'Unit', 'deg');
nlgr = setpar(nlgr, 'Name', 'D');
nlgr = setpar(nlgr, 'Unit', 'constant');


nlgr = setinit(nlgr, 'Fixed',false);
opt = nlgreyestOptions('Display', 'off','SearchMethod','lsqnonlin',...
    'SearchOption',optimset('lsqnonlin'));

nlgr = nlgreyest(data, nlgr, opt);

end