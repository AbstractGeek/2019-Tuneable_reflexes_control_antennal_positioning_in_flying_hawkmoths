function [returnCurves, perturbStatus] = extractAntennalReturnCurves(inputDir, outputDir,...
    treatment_name, angle_name, filter_vals)
%
% Extracts antennal return curves, subtracts the wing beat frequency from
% them and saves it for the next step.
%
% Dinesh Natesan
% 27th Feb 2016

% If no arguments, set default inputs
if (nargin < 1)
    inputDir = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Analysis/mat files';
    outputDir = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Analysis';
    % angle_name = {'HeadVect','HeadVect','HeadVect','HeadVect'};
    % treatment_name = {'Control','JO','Control','JO'};
    % filter_vals = [false, false, true, true];
    angle_name = {'HeadVect','HeadVect'};
    treatment_name = {'Control','JO'};
    filter_vals = [false, false];
    % angle_name = {'Perturbed','Perturbed','HeadVect','HeadVect'};
    % angle_name = {'HeadVect','HeadVect','Perturbed','Perturbed'};
    % treatment_name = {'Control','JO','Control','JO'};
end

% Switch off warning
warning('off','MATLAB:legend:IgnoringExtraEntries');  % Switch off warning temporarily.

% Plotting related defaults
figure_size = [1 21 29.7 21];       % A4 size
windspeed = [0,1.5,2.5,4];
wsNames = {'ws_0_mps'
    'ws_1_5mps'
    'ws_2_5mps'
    'ws_4_mps'};
stepNames = arrayfun(@(x) sprintf('Step-%g',x),[1,2,3,4,5],...
    'UniformOutput',false);
legend_names = {'Raw','Notch','Sine-Fit-1','Sine-Fit-2'};
plotData = 1;

% Set misc default parameters
folder_name =  'returnCurves (Filtered)';
%experiment_exceptions = {'moth3_16jul15','moth4_17jul15','moth7_23jun15','moth8_23jun15'};
experiment_exceptions = cell(5,5);
experiment_exceptions(:,1) = {'moth3_16jul15','moth4_17jul15','moth3_17jul15','moth7_23jun15','moth8_23jun15'}';
experiment_exceptions(:,2) = {true(5,1), true(5,1), false(5,1), false(5,1), false(5,1)}';
experiment_exceptions(:,3) = {true(5,1), true(5,1), false(5,1), [false,true,true,true,true]', false(5,1)}';
experiment_exceptions(:,4) = {true(5,1), true(5,1), false(5,1), [true, false, true, true, true]', false(5,1)}';
experiment_exceptions(:,5) = {true(5,1), true(5,1), true(5,1), [true, true, false, false, false], true(5,1)}';

experiment_exceptions = cell2table(experiment_exceptions,...
    'VariableNames',[{'Experiments'};wsNames]);
% experiment_selections = {'moth3_30jun15','moth3_31aug15','moth2_1jul15','moth3_2sep15'};

% Filtering constants
sam_freq = 1000;
cutoff_freq = 150;
notch_width = 5;
std_cutoff_freq = 3;
std_cutoff_start_point = 3;
mean_separation_cutoff = 0.90;
mean_cutoff = 0.25;
filter_order = 4;

% Obtain complete response return frames
complete_response_return_frames = cell2mat(arrayfun(@(x,y)...
    x+1:x+499,...
    [1500,2500,3500,4500,5500],'UniformOutput',false)')';

% Obtain parameters not to be cleared
initvars = who;
initvars = [initvars;'mf';'initvars'];

for mf=1:length(angle_name)
    % Clear parameters and load new dataset
    clearvars('-except',initvars{:})
    CalculatedAngles = load(fullfile(inputDir,...
        sprintf('CalculatedAngles_%s_%s_%s.mat',...
        iff(filter_vals(mf),'Filtered','Unfiltered'),...
        angle_name{mf},treatment_name{mf})));
    
    %% Get return curves
    perturb_end = [1500,2500,3500,4500,5500];
    % window = [100,300,160];
    % window = [100,499,160];
    % window = [100,400,250];
    % time = (0:sum(window(1:2))-1)';
    
    rawReturnCurves = struct;
    filtReturnCurves = struct;
    outDir = fullfile(outputDir,folder_name, sprintf('%s (%s-%s)',...
        treatment_name{mf},...
        iff(filter_vals(mf),'Filtered','Unfiltered'),angle_name{mf}));
    if ~isdir(outDir)
        mkdir(outDir);
    end
    h = figure('Units','centimeters','Position',figure_size,'Visible','off');
    
    % Start loop
    experiments = fieldnames(CalculatedAngles.projAngles);
    
    for i=1:length(experiments)
        
        % check if the experiment is an exception
        % if ~any(strcmp(experiment_selections,experiments{i}))
        if any(strcmp(experiment_exceptions.Experiments,experiments{i}))
            exception_flag = true;
            exception_index = find(strcmp(experiment_exceptions.Experiments,experiments{i}));
        else
            exception_flag = false;            
        end
        
        % Create a csv to save filtered sine waves
        fileID = fopen(fullfile(outDir,sprintf('%s.csv',experiments{i})),'w');
        
        % Get trials and create save directory
        trials = fieldnames(CalculatedAngles.projAngles.(experiments{i}));
        outdirExp = fullfile(outDir,(experiments{i}));
        if ~isdir(outdirExp)
            mkdir(outdirExp);
        end
        
        for j=1:length(trials)
            
            windspeeds = fieldnames(CalculatedAngles.projAngles.(experiments{i}).(trials{j}));
            % Initialize clean perturbation table
            perturbStatus.(experiments{i}).(trials{j}) = ...
                table(false(5,1),false(5,1),false(5,1),false(5,1),...
                'VariableNames',wsNames,'RowNames',stepNames);
            % Initialize alternate perturbation table (just in case)
            altPerturbStatus.(experiments{i}).(trials{j}) = ...
                table(false(5,1),false(5,1),false(5,1),false(5,1),...
                'VariableNames',wsNames,'RowNames',stepNames);
            
            % Check if shortened mode
            digitized_data_length = length(find(all(...
                ~isnan(CalculatedAngles.projAngles.(experiments{i}).(trials{j}).(windspeeds{1})),2)));
            if digitized_data_length < 2000
                fprintf('Found less than 2000 digitized frames in %s-%s. Skipping current trial.\n',...
                    experiments{i},trials{j});
                continue;
            elseif digitized_data_length > 2000 && digitized_data_length < 3000
                complete_response_mode = true;
            else
                complete_response_mode = false;
            end
            
            % Set window and define time
            window = [100,300,250];
            time = (0:sum(window(1:2))-1)';
            
            % Run through individual airflows
            for k=1:length(windspeeds)
                
                % Write into the csv file
                fprintf(fileID,'%s_%s\n',trials{j},windspeeds{k});
                
                % Initialize cell structures
                % Raw and filtered return responses
                rawReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k}) = ...
                    NaN(sum(window(1:2)),5);
                filtReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k}) = ...
                    NaN(sum(window(1:2)),5);
                returnCurveVars.(experiments{i}).(trials{j}).(windspeeds{k}).sinefit1 = ...
                    NaN(sum(window(1:2)),5);
                returnCurveVars.(experiments{i}).(trials{j}).(windspeeds{k}).sinefit2= ...
                    NaN(sum(window(1:2)),5);
                returnCurveVars.(experiments{i}).(trials{j}).(windspeeds{k}).selectedfit = ...
                    NaN(sum(window(1:2)),5);
                % Start points
                startPoints.(experiments{i}).(trials{j}).(windspeeds{k}) = ...
                    table(NaN(5,1),NaN(5,1),NaN(5,1),NaN(5,1),...
                    'VariableNames',{'Std_first','Std_last','Mean_first','Mean_last'});
                % Extracted return curves
                returnCurves.(experiments{i}).(trials{j}).(windspeeds{k}) = cell(5,1);
                altReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k}) = cell(5,1);
                
                % Filter if complete data
                if ~complete_response_mode
                    processedAngles = ButterFilt(interpolateMissingPoints(...
                        CalculatedAngles.projAngles.(experiments{i}).(trials{j}).(windspeeds{k})),...
                        sam_freq, cutoff_freq);
                end
                
                % Loop through individual perturbations
                for l=1:length(perturb_end)
                    
                    % Check if it falls in the exception category
                    if exception_flag
                       % Check if the step in this airflow falls in the category
                       if experiment_exceptions.(wsNames{k}){exception_index}(l)
                           % Skip step
                           continue;
                       end
                    end
                    
                    % Obtain current return trajectory
                    if complete_response_mode
                        current = CalculatedAngles.projAngles.(experiments{i}).(trials{j}).(windspeeds{k})...
                            (complete_response_return_frames(:,l),2);
                        if (any(isnan(current)) && ~all(isnan(current)))
                            % Try to interpolate and obtain points
                            current = CalculatedAngles.projAngles.(experiments{i}).(trials{j}).(windspeeds{k})...
                            (complete_response_return_frames(:,l),2);
                            % nans
                            current_nans = isnan(current);
                            % interpolate (simple spline)
                            current = spline(complete_response_return_frames(~current_nans,l),...
                                current(~current_nans), complete_response_return_frames(:,l));                            
                        end
                        current = ButterFilt(current,sam_freq,cutoff_freq);
                    else
                        current = processedAngles(perturb_end(l)-window(1)+1:perturb_end(l)+499,2);
                    end
                    
                    % Check if it is filled with NaN
                    if all(isnan(current))
                        fprintf('[%s-%s-%0.1fmps-%d] All NaN trajectory - skipping\n',...
                        experiments{i},trials{j},windspeed(k),l);
                        continue;
                    end
                    
                    
                    %% Fourier based wing beat frequency filtering
                    % Obtain change in current return trajectory and use it
                    % to get peak frequencies above 10 Hz, and above the
                    % mean + std_cutoff_freq*std(fsig)
                    change = diff(current);
                    [fsig,f] = FourierTransform(change,sam_freq,1);
                    fsig = abs(fsig);   % Ignore phase information
                    
                    % threshold everything below the cutoff
                    nfsig = abs(fsig);
                    upper = mean(fsig)+std_cutoff_freq*std(fsig);
                    lower = mean(fsig)-std_cutoff_freq*std(fsig);
                    nfsig(abs(fsig)<upper & abs(fsig)>lower) = 0;
                    
                    % Find the top two wing beat frequencies
                    freqs = nan(4,1);
                    % Find the top frequency in the wing beat range                    
                    freq_range = f((f>25)' & (f<50)')';
                    freq_signal = fsig((f>25)' & (f<50)');
                    [~,Ival] = max(freq_signal);
                    freqs(1) = freq_range(Ival);                    
                    % Find the next freq in the wing beat range
                    freq_signal = freq_signal(freq_range~=freqs(1));
                    freq_range = freq_range(freq_range~=freqs(1));                    
                    [~,Ival] = max(freq_signal);
                    freqs(2) = freq_range(Ival);                    
                    % Add two harmonics to the frequencies
                    freqs(3:4) = 2.*freqs(1:2);                    
                    
                    % Copy change into another variable for modification
                    fxn = current;
                    % define x,y for easy reference
                    x = (0:length(change)-1)'./sam_freq;
                    y = change;
                    
                    for m=1:length(freqs)
                        % filter using notch (notch width and filter order
                        % set above in defaults)
                        F0 = freqs(m)*2/sam_freq;   % normalized freq w.r.t nyquist
                        BW = notch_width/sam_freq;
                        d = design(fdesign.notch('N,F0,BW',filter_order,F0,BW));
                        sosg = coeffs(d);
                        
                        % Pad and filter (Padding helps in reducing edge
                        % effects)
                        padded = padarray(fxn,[size(current,1) 0],'symmetric');
                        pfxn = filtfilt(sosg.SOSMatrix,sosg.ScaleValues,padded);
                        
                        % Save filtered waveform as fxn
                        fxn = pfxn(size(current,1)+1:size(current,1)+size(current,1));
                    end
                    
                    %% Fit based wing beat frequency filtering
                    % Find sine waves of high powers in the wing beat
                    % frequency range and save them.
                    
                    % Find the first sine wave in the wing beat frequency range
                    freq_range1 = f((f>25)' & (f<50)')';
                    freq_signal1 = fsig((f>25)' & (f<50)');
                    [Mval1,Ival1] = max(freq_signal1);
                    fit_freq1 = freq_range1(Ival1);
                    fit_freq_sig1 = Mval1;
                    
                    % Fit the first sine wave to the saved frequency and its harmonic,
                    % and subtract them from the waveform.
                    myfit1 = fittype(sprintf('a1*sin(%f*x+c1)+a2*sin(%f*x+c2)',...
                        2*pi*fit_freq1,4*pi*fit_freq1));
                    ffit1 = fit(x,y,myfit1,'Startpoint',[0,0,0,0]);
                    fxf1 = change - ffit1(x);                    
                    
                    % Find the second sine wave in the wing beat frequency range
                    % Perform fourier again to better fit the changed
                    % waveform (instead of choosing two at the start).
                    [fsig_ssv,f_ssv] = FourierTransform(fxf1,sam_freq,1); %ssv = second sine wave
                    fsig_ssv = abs(fsig_ssv);   % Ignore phase information
                    freq_range2 = f_ssv((f_ssv>25)' & (f_ssv<50)')';
                    freq_signal2 = fsig_ssv((f_ssv>25)' & (f_ssv<50)');
                    [Mval2,Ival2] = max(freq_signal2);
                    fit_freq2 = freq_range2(Ival2);
                    fit_freq_sig2 = Mval2;
                    
                    % Subtract second sine wave and its first harmonic
                    myfit2 = fittype(sprintf('a1*sin(%f*x+c1)+a2*sin(%f*x+c2)',...
                        2*pi*fit_freq2,4*pi*fit_freq2));
                    ffit2 = fit(x,fxf1,myfit2,'Startpoint',[0,0,0,0]);
                    fxf2 = fxf1 - ffit2(x);
                                        
                    % Print on display and on the csv
                    fprintf('[%s-%s-%0.1fmps-%d]Sines filtered: %g,%g,%g,%g hz \n',...
                        experiments{i},trials{j},windspeed(k),l,...
                        fit_freq1,fit_freq1*2,fit_freq2,fit_freq2*2);
                    fprintf(fileID,'%g,%g,%g,%g\n',...
                        fit_freq1,fit_freq1*2,fit_freq2,fit_freq2*2);
                    
                    %% Combine data together                    
                    % Get pieces
                    perturbed_angle = CalculatedAngles.projAngles.(experiments{i}).(trials{j}).(windspeeds{k})...
                        (perturb_end(l)-window(1)+1:perturb_end(l),2);
                    % realign filtered data (add the perturbed mean value
                    % to the changes as all the changes are with respect to
                    % the mean)
                    % perturbed_mean = nanmean(perturbed_angle);
                    notch_filt = fxn;
                    fit_filt1 = [current(1);current(1)+cumsum(fxf1)];
                    fit_filt2 = [current(1);current(1)+cumsum(fxf2)];
                    
                    if complete_response_mode
                        % Combine them
                        current = [perturbed_angle;current(1:window(2))];
                        notch_filt = [perturbed_angle;notch_filt(1:window(2))];
                        fit_filt1 = [perturbed_angle;fit_filt1(1:window(2))];
                        fit_filt2 = [perturbed_angle;fit_filt2(1:window(2))];
                        % create separate time array for change values
                        time_change = (0:length(change)-1)'./sam_freq;
                    else
                        % Combine them
                        current = current(1:window(1)+window(2));
                        notch_filt = notch_filt(1:window(1)+window(2));
                        fit_filt1 = fit_filt1(1:window(1)+window(2));
                        fit_filt2 = fit_filt2(1:window(1)+window(2));
                        
                        % create separate time array for change values
                        time_change = (0:length(change)-1)'./sam_freq;
                    end

                    
                    %% Find start point of return
                    % Define which filtered trajectory to use for
                    % start_point selection
                    fit_filt = fit_filt2;
                    
                    % Find the mean antennal in th perturbed state and the
                    % steady state
                    % perturbed mean
                    perturbed = fit_filt(1:window(1));
                    perturbed_mean = nanmean(perturbed);
                    perturbed_std = nanstd(perturbed);
                    % unperturbed (steady state) mean
                    unperturbed = fit_filt(window(3)+window(1)+1:end);
                    unperturbed_mean = nanmean(unperturbed);
                    unperturbed_std = nanstd(unperturbed);
                    % return response
                    return_response = fit_filt(window(1)+1:window(3)+window(1));
                    
                    % Check if perturbed mean > unperturbed mean (as
                    % it should traditionally be). Additionally check if
                    % the perturbed and unperturbed mean are significantly
                    % apart using the mean_separation cutoff
                    if (perturbed_mean > unperturbed_mean) && ...
                            (mean_separation_cutoff*perturbed_mean > unperturbed_mean)
                        
                        % Try multiple ways of defining the start point.
                        
                        % 1) Standard deviation based startpoint selection
                        lower_cutoff = perturbed_mean-std_cutoff_start_point*perturbed_std;
                        
                        % Is this selection necessary? - Include it later,
                        % if necessary
                        %if (lower_cutoff > unperturbed_mean)
                            % Find the first point and the last point that
                            % cross the lower cutoff in the return response
                            % (works! double checked algo on a fake
                            % sinusoid + exponential waveform)                            
                            start_point_std_1 = find(return_response<lower_cutoff,1,'first')+window(1);
                            start_point_std_2 = find(diff(return_response<lower_cutoff)==1,1,'last')+window(1)+1;
                            perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = true;
                            
                            % handle empty start_point_mean_2
                            if isempty(start_point_std_2)
                                % The first point itself is past the
                                % cutoff.
                                start_point_std_2 = start_point_std_1;
                            end
                            
                        %else
                            % Most probably not perturbed
                        %    perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = false;
                        %end
                        
                        % 2) Mean based startpoint selection
                        cutoff = perturbed_mean - ...
                            (perturbed_mean - unperturbed_mean)*mean_cutoff;
                        
                        % First point which crosses the cutoff
                        start_point_mean_1 = find(return_response<cutoff,1,'first')+window(1);
                        % Last point which crosses the cutoff
                        start_point_mean_2 = find(diff(return_response<cutoff)==1,1,'last')+window(1)+1;
                        
                        % handle empty start_point_mean_2
                        if isempty(start_point_mean_2)
                            % The first point itself is past the
                            % cutoff.
                            start_point_mean_2 = start_point_mean_1;
                        end
                        
                        % Note that the start point based mean is only
                        % valid if perturbStatus is set to true above (in
                        % the std based start point detection)
                        
                    elseif (unperturbed_mean > perturbed_mean) && ...
                            (mean_separation_cutoff*unperturbed_mean > perturbed_mean)
                        % Handles the not so intuitive case of unperturbed
                        % mean being greater than the perturbed mean, and
                        % qualifying the mean separation cutoff
                        
                        % Clean perturbation status is set to false and the
                        % alternate perturbation status is evaluated here.
                        perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = false;
                        
                        % Try multiple ways of defining the start point.
                        % 1) Standard deviation based startpoint selection
                        upper_cutoff = perturbed_mean+std_cutoff_start_point*perturbed_std;
                        
                        % Is this selection necessary? - Include it later,
                        % if necessary
                        % if (upper_cutoff < unperturbed_mean)
                            start_point_std_1 = find(return_response>upper_cutoff,1,'first')+window(1);
                            start_point_std_2 = find(diff(return_response>upper_cutoff)==1,1,'last')+window(1)+1;
                            altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = true;
                            
                            % handle empty start_point_mean_2
                            if isempty(start_point_std_2)
                                % The first point itself is past the
                                % cutoff.
                                start_point_std_2 = start_point_std_1;
                            end
                       % else
                            % Most probably not perturbed
                       %     altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = false;
                       % end
                        
                        % 2) Mean based startpoint selection
                        cutoff = perturbed_mean + ...
                            (unperturbed_mean - perturbed_mean)*mean_cutoff;
                        
                        % First point which crosses the cutoff
                        start_point_mean_1 = find(return_response>cutoff,1,'first')+window(1);
                        % Last point which crosses the cutoff
                        start_point_mean_2 = find(diff(return_response>cutoff)==1,1,'last')+window(1)+1;
                        
                        % handle empty start_point_mean_2
                        if isempty(start_point_mean_2)
                            % The first point itself is past the
                            % cutoff.
                            start_point_mean_2 = start_point_mean_1;
                        end
                        
                        % Note that the start point based mean is only
                        % valid if altPerturbStatus is set to true above (in
                        % the std based start point detection)
                        
                    else
                        % Most probably no perturbation
                        perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = false;
                        altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) = false;
                    end
                    
                    %% Misc calculations (for plotting)
                    % Obtain cutoffs
                    perturbed_upper_cutoff = perturbed_mean+std_cutoff_start_point*perturbed_std;   % Same as upper cutoff
                    perturbed_lower_cutoff = perturbed_mean-std_cutoff_start_point*perturbed_std;   % Same as lower cutoff
                    unperturbed_upper_cutoff = unperturbed_mean+std_cutoff_start_point*unperturbed_std;
                    unperturbed_lower_cutoff = unperturbed_mean-std_cutoff_start_point*unperturbed_std;
                    
                    %% Plot Data
                    if plotData
                        
                        % Plot timeseries, means and start points
                        subplot(3,1,1),hold on;
                        title('Return curves')
                        
                        % Draw time series
                        plot(time,current,'LineWidth',2.0);
                        plot(time,notch_filt,'LineWidth',2.0);
                        plot(time,fit_filt1,'LineWidth',2.0);
                        plot(time,fit_filt2,'LineWidth',2.0);
                        legend(legend_names,'Box','off','Location','best');
                        
                        % Draw perturbed lines
                        plot([0 time(window(1))],[perturbed_mean,perturbed_mean],'--r','LineWidth',2.0);
                        plot([time(window(1)) time(window(2))],[perturbed_mean,perturbed_mean],'--r');
                        plot([0 time(end)],[perturbed_upper_cutoff,perturbed_upper_cutoff],':r');
                        plot([0 time(end)],[perturbed_lower_cutoff,perturbed_lower_cutoff],':r');
                        
                        % Draw unperturbed lines
                        plot([time(window(3)+window(1)) time(end)],[unperturbed_mean,unperturbed_mean],...
                            '--g','LineWidth',2.0);
                        plot([0 time(window(3)+window(1))],[unperturbed_mean,unperturbed_mean],'--g');
                        plot([0 time(end)],[unperturbed_upper_cutoff,unperturbed_upper_cutoff],':g');
                        plot([0 time(end)],[unperturbed_lower_cutoff,unperturbed_lower_cutoff],':g');
                        
                        if (~perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l)) && ...
                                (~altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l))
                            text(time(window(1)+1),fit_filt(window(1)+1), ...
                                '\leftarrow random noise elicited startpoint. No perturbation occurred here');
                        elseif (altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l))
                            text(time(window(1)+1),fit_filt(window(1)+1), ...
                                '\leftarrow alternate perturbation occurred here');
                        else
                            % Draw start points
                            plot(time(start_point_std_1),fit_filt(start_point_std_1),'o',...
                                'MarkerSize',10,'Color',rgb('navy blue'),'LineWidth',2.0);
                            plot(time(start_point_std_2),fit_filt(start_point_std_2),'d',...
                                'MarkerSize',10,'Color',rgb('navy blue'),'LineWidth',2.0);
                            plot(time(start_point_mean_1),fit_filt(start_point_mean_1),'o',...
                                'MarkerSize',10,'Color',rgb('purple'),'LineWidth',2.0);
                            plot(time(start_point_mean_2),fit_filt(start_point_mean_2),'d',...
                                'MarkerSize',10,'Color',rgb('purple'),'LineWidth',2.0);
                        end
                        
                        % Plot derivative of time series
                        subplot(3,1,2),hold on;
                        title('Derivative of return curves');
                        plot(time_change,change,'LineWidth',2.0);
                        %plot(time_change,fxn,'LineWidth',2.0);
                        plot(time_change,fxf1,'LineWidth',2.0);
                        plot(time_change,fxf2,'LineWidth',2.0);
                        legend(legend_names([1,3:4]),...
                            'Box','off','Location','best');
                        
                        % Plot fourier series spectrum
                        subplot(3,1,3),hold on;
                        fline = [-sam_freq/2,sam_freq/2];
                        plot(f,fsig);
                        % Fourier based filtering
                        plot(fline,[mean(fsig),mean(fsig)]);
                        plot(fline,[upper,upper]);
                        plot(f,nfsig);
                        plot(sort(freqs),fsig(ismember(f,freqs)),'+g');
                        %plot(f(indd(selected)),lmval(selected),'+g');
                        % Fit based filtering
                        plot(fit_freq1,fit_freq_sig1,'sk'); % First sine
                        [~,findd] = min(abs(f-2*fit_freq1)); % Its harmonic
                        plot(f(findd),fsig(findd),'sm'); % Plot first harmonic
                        plot(fit_freq2,fit_freq_sig2,'dk'); % Second sine
                        [~,findd] = min(abs(f-2*fit_freq2)); % Its harmonic
                        plot(f(findd),fsig(findd),'dm'); % Second harmonic
                        legend({'Raw','Mean','Threshold',...
                            'Thesholded','Selected',...
                            'Sine-1','Sine-1-Har-1','Sine-2','Sine-2-Har-1'},...
                            'Box','off','Location','best');
                        
                        % Overll Title
                        figtitle(strjoin(strsplit(sprintf('%s-%s-%s-%d',...
                            experiments{i},trials{j},windspeeds{k},l),'_'),'-'),...
                            'fontweight','bold');
                        
                        % save into plots
                        export_fig(fullfile(outdirExp,...
                            sprintf('%s_%s_%d',trials{j},windspeeds{k},l)),...
                            '-transparent',h);
                        export_fig(fullfile(outDir,sprintf('%s',experiments{i})),...
                            '-pdf','-append','-transparent',h);
                        
                        clf(h);
                    end
                    
                    %% Save Extracted curves
                    
                    % Save returncurve if only perturbed
                    if (perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l) || ...
                            altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l))
                        
                        % Preferred startpoint
                        pref_start_point = start_point_mean_2;
                        
                        % Generate return curve and the associated set
                        % point
                        returnCurveData = [fit_filt(pref_start_point:end),...
                            repmat(unperturbed_mean,length(fit_filt)-pref_start_point+1,1)];
                        
                        if (perturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l))
                            % Save to return curve structure
                            returnCurves.(experiments{i}).(trials{j}).(windspeeds{k}){l,1} = ...
                                returnCurveData;
                            
                        elseif altPerturbStatus.(experiments{i}).(trials{j}).(windspeeds{k})(l)
                            % Save to alternate return curve structure
                            altReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k}){l,1} = ...
                                returnCurveData;
                        end

                        % save start points
                        startPoints.(experiments{i}).(trials{j}).(windspeeds{k}).Std_first(l,1) = ...
                            iff(isempty(start_point_std_1),NaN,start_point_std_1);
                        startPoints.(experiments{i}).(trials{j}).(windspeeds{k}).Std_last(l,1) = ...
                            iff(isempty(start_point_std_2),NaN,start_point_std_2);
                        startPoints.(experiments{i}).(trials{j}).(windspeeds{k}).Mean_first(l,1) = ...
                            start_point_mean_1;
                        startPoints.(experiments{i}).(trials{j}).(windspeeds{k}).Mean_last(l,1) = ...
                            start_point_mean_2;
                        
                    end
                    
                    % save raw, filtered return curves into struct
                    rawReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k})(:,l) = current;
                    filtReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k})(:,l) = notch_filt;
                    returnCurveVars.(experiments{i}).(trials{j}).(windspeeds{k}).sinefit1(:,l) = fit_filt1;
                    returnCurveVars.(experiments{i}).(trials{j}).(windspeeds{k}).sinefit2(:,l) = fit_filt2;
                    returnCurveVars.(experiments{i}).(trials{j}).(windspeeds{k}).selectedfit(:,l) = fit_filt;
                    
                    
                end
            end            
        end
        % Experiment done, close fileID
        fclose(fileID);
    end
    
    close(h);
    save(fullfile(inputDir,...
        sprintf('ReturnCurves_%s_%s_%s.mat',...
        iff(filter_vals(mf),'Filtered','Unfiltered'),...
        angle_name{mf},treatment_name{mf})),...
        'rawReturnCurves','filtReturnCurves','startPoints',...
        'returnCurves','returnCurveVars','perturbStatus');
    
    warning('on','MATLAB:legend:IgnoringExtraEntries');
    
end