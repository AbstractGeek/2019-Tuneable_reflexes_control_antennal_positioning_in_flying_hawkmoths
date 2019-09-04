% Analyze the intrinsic setpoint data using control theory
clear, close all;
% Set defaults
trials = 5;
csv_dir = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Simulation/Antennal-motor-system-neural-network/Simulated-Data/Intrinsic-setpoint';
matfolder = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Simulation/Antennal-motor-system-neural-network/Control-Theoretic-Analysis/mat files';
outfolder = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Simulation/Antennal-motor-system-neural-network/Control-Theoretic-Analysis/';
% extraction defaults
perturb_end = [151, 251, 351, 451, 551];
step_size = 49;
mean_window = [300, 50]; % Total size 300, 50 mean unperturbed
mean_cutoff = 0.25;
perturbations = [0,25,75,100];
% other
angle_name = {'SimPosIn'};
filter_vals = [false];
wsNames = {'ws_0_mps'
    'ws_1_5mps'
    'ws_2_5mps'
    'ws_4_mps'};    % Mimicing perturbations as airflow (for now)
stepNames = arrayfun(@(x) sprintf('Step-%g',x),[1,2,3,4,5],...
    'UniformOutput',false);

%% Extract data (mimic extractAntennalReturnCurves function)
treatment_name = dir(csv_dir);
treatment_name = {treatment_name([treatment_name.isdir]).name}';
treatment_name(ismember(treatment_name, {'.','..'})) = [];

for trt = 1:length(treatment_name)
    curr_trt = strrep(strrep(treatment_name{trt},'.','_'),'-','_');
    % Create data structure
    returnCurves.(curr_trt) = struct;
    perturbStatus.(curr_trt) = struct;
    rawData.(curr_trt) = struct;

    % File handling
    % Obtain csv files
    csv_files = dir(fullfile(csv_dir, treatment_name{trt}, '*.csv'));
    csv_files = {csv_files(:).name}';
    csv_files(ismember(csv_files, {'Simulation-details.csv'}))=[];
    
    % Begin extraction
    time = (0:10:490)';
    new_time = (0:1:499)';
    for i=1:length(csv_files)
        curr_data = readtable(fullfile(csv_dir, treatment_name{trt}, csv_files{i}));
        % Obtain trial number and fake airflow name
        [~,tok] = regexp(csv_files{i},...
            'P-(\d*)_T-(\d*).csv', 'match','tokens');
        trial_no = str2double(tok{1}{2});
        windspeed = wsNames{ismember(perturbations, str2double(tok{1}{1}))};
        % Create data structures if necessary
        if ~isfield(returnCurves.(curr_trt), sprintf('Trial_%d',trial_no))
            returnCurves.(curr_trt).(sprintf('Trial_%d',trial_no)) = struct;
            for ws=1:length(wsNames)
                returnCurves.(curr_trt).(sprintf('Trial_%d',trial_no)).(wsNames{ws}) = cell(trials,1);
            end
            perturbStatus.(curr_trt).(sprintf('Trial_%d',trial_no)) = ...
                table(true(5,1),true(5,1),true(5,1),true(5,1),...
                'VariableNames',wsNames,'RowNames',stepNames);            
        end
        
        % Save rawData
        rawData.(curr_trt).(sprintf('Trial_%d',trial_no)).(windspeed) = curr_data;
        
        % Extract steps
        for j=1:length(perturb_end)
            perturbed_data = curr_data.position(perturb_end(j)-step_size:perturb_end(j));
            unperturbed_data = spline(time, ...
                curr_data.position(perturb_end(j):perturb_end(j)+step_size), new_time);
            % Obtain cutoff
            perturbed_mean = mean(perturbed_data);
            unperturbed_mean = mean(unperturbed_data(...
                mean_window(1)-mean_window(2)+1:mean_window(1)));
            cutoff = perturbed_mean - ...
                (perturbed_mean - unperturbed_mean)*mean_cutoff;
            % Obtain start point
            if (perturbed_mean > unperturbed_mean)
                start_point_mean = find(diff(unperturbed_data<cutoff)==1,1,'last')+1;
            elseif (perturbed_mean < unperturbed_mean)
                start_point_mean = find(diff(unperturbed_data>cutoff)==1,1,'last')+1;
            else
                error('Weird. No perturbation in a simulation?');
            end
            
            % Save dataset
            returnCurves.(curr_trt).(sprintf('Trial_%d',trial_no)).(windspeed){j,1} = ...
                [unperturbed_data(start_point_mean:mean_window(1)),...
                repmat(unperturbed_mean,length(unperturbed_data(start_point_mean:mean_window(1))),1)];
            
            
        end
        
    end
    
    % Save as matfile
    save(fullfile(matfolder,...
        sprintf('ReturnCurves_%s_%s_%s.mat',...
        iff(filter_vals(1),'Filtered','Unfiltered'),...
        angle_name{1},treatment_name{trt})),...
        'returnCurves','perturbStatus');
    clearvars returnCurves perturbStatus;
    
    % Save raw data
    save(fullfile(matfolder,...
        sprintf('RawData_%s_%s_%s.mat',...
        iff(filter_vals(1),'Filtered','Unfiltered'),...
        angle_name{1},treatment_name{trt})),...
        'rawData');    
    clearvars rawData;
end



%% Analyze the return curves (add it as batch jobs)

for i=1:length(treatment_name)
%     [ModelFits, ModelData, ExpFits] = analyzeReturnCurves(matfolder, outfolder,...
%         treatment_name, angle_name, filter_vals);
    batch('analyzeReturnCurves',3 ,{matfolder, outfolder,...
        treatment_name(i), angle_name, filter_vals});
end