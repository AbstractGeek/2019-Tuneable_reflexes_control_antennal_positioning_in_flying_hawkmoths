%% Create a simple exponential return
time_constant = 200;    % in milliseconds
gain = 5;
set_point = 20;         % some arbitary input
sam_freq = 1000;

% Create input, output data
time_data = (0:1/sam_freq:1)';
error_data = gain * exp(-1 * time_data / (time_constant/sam_freq)); % a = 10, b = -5
input_data = set_point * ones( size ( time_data ) );
output_data = input_data - error_data;


%% Try fitting control models to it
modelpzs = [0,0; 1,0; 1,1; 2,0; 2,2]; % P, I, PI, PID
l_modelfns = {'l_proportional','l_integrator','l_proportional_integral','l_double_integrator','l_proportional_integral_differential'};

sys = cell(5,1);

% h = figure('Units','centimeters','Position',[1 21 29.7 21],'Visible','on');
% hold on;

for m=1:size(modelpzs,1)
    [sys{m},data] = feval(l_modelfns{m},input_data, output_data, sam_freq, 'estimate');       
%     bodeplot(sys{m});    
end

% legend(l_modelfns);
figure('Units','centimeters','Position',[1 21 29.7 21],'Visible','on');
compare(data,sys{1},sys{2},sys{3},sys{4},sys{5});

%% Print parameters
for m=1:size(modelpzs,1)
 fprintf('[%s %s] \n',...
        num2str(sys{m}.Report.Parameters.X0','%0.2f '),...
        num2str(sys{m}.Report.Parameters.ParVector','%0.2f '));
    
end