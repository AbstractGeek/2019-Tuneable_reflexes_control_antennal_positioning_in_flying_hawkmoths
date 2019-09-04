function [nlgr,data] = nl_proportional_integral(input, output, sam_freq)
% function [nlgr,data] = nl_proportional_integral(input, output, sam_freq)
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
FileName = 'nl_proportional_integral_model';
Order = [1 3 1];
Parameters = [1,1,1,1];
InitialStates = output(1);
Ts = 0;
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'nl_proportional_integral');

% Set Input and Output names
set(nlgr, 'InputName', {'Set-Point', 'Antennal Angle','ReleaseState'}, ...
    'InputUnit',  {'deg', 'deg','boolean'},            ...
    'OutputName', 'Antennal Angle', ...
    'OutputUnit', 'deg',                         ...
    'TimeUnit', 's');

% Set init and par names
nlgr = setinit(nlgr, 'Name', 'Antennal Angle');
nlgr = setinit(nlgr, 'Unit', 'deg');
nlgr = setpar(nlgr, 'Name', {'A','B','C','D'});
nlgr = setpar(nlgr, 'Unit', {'constant', 'constant', 'constant','constant'});


nlgr = setinit(nlgr, 'Fixed',false);
opt = nlgreyestOptions('Display', 'off','SearchMethod','lsqnonlin',...
    'SearchOption',optimset('lsqnonlin'));

nlgr = nlgreyest(data, nlgr, opt);

end