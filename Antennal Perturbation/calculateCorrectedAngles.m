function [] = calculateCorrectedAngles(inputDir_pref, outputDir_pref,...
    correctedAngleDir_pref, dataset_names, angle_names, filter_vals)
% function [] = calculateCorrectedAngles(inputDir_pref, outputDir_pref,...
%    correctedAngleDir_pref, dataset_names, angle_names)
%
% Calculates left and right antennal angle, corrected angles and saves them
% as a mat file
% 
% Dinesh Natesan
% Last Modified: 24th Feb 2017

if (nargin < 1)
    inputDir_pref = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Digitized Data/Completed Digitization';
    outputDir_pref = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Analysis/mat files';
    correctedAngleDir_pref = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Analysis/Corrected Angles';
%     dataset_names = {'Control','JO','Control','JO'};
%     angle_names = {'Perturbed', 'Perturbed', 'Perturbed', 'Perturbed'};    
%     filter_vals = [false, false, true, true];
    dataset_names = {'Control','JO','Control','JO'};
    angle_names = {'HeadVect', 'HeadVect', 'HeadVect', 'HeadVect'};    
    filter_vals = [false, false, true, true];
%     dataset_names = {'Control','JO'};
%     angle_names = {'HeadVectAP', 'HeadVectAP'};    
%     filter_vals = [false, false];
end

initvars = who;
initvars = [initvars;'ds';'initvars'];

%% Calculate corrected headvect antennal angles
for ds = 1:length(dataset_names)
    % Load dataset
    Dataset = load(fullfile(inputDir_pref,strcat(dataset_names{ds},'.mat')));
    TreatmentType = dataset_names{ds};
    AngleType = angle_names{ds};
    FilterFlag = filter_vals(ds);
   
    % Begin
    windspeed = [0,1.5,2.5,4];
    control_start = [1500,2500,3500,4500,5500];
    mean_window = [301, 499];
    experiments = fieldnames(Dataset);
    % Output Strucures
    max_num_moths = 20;
    projAngles = struct;
    antAngles = struct;
    sphAngles = struct;
    debugVars = struct;
    
    % % corrected angle - from means
    curr_ind = 0;
    mean_lca_means = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    stder_lca_means = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    % Right corrected angle
    mean_rca_means = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    stder_rca_means = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    % % corrected angle - from distribution
    mean_lca_dist = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    stder_lca_dist = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    % Right corrected angle
    mean_rca_dist = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    stder_rca_dist = table(cell(max_num_moths,1),cell(max_num_moths,1),...
        NaN(max_num_moths,1),NaN(max_num_moths,1),NaN(max_num_moths,1),...
        NaN(max_num_moths,1),'VariableNames',...
        {'Moth','Trial','ws_0_mps','ws_1_5mps','ws_2_5mps','ws_4_mps'});
    % csv output folder
    csvfolder = fullfile(correctedAngleDir_pref,...
        sprintf('CorrectedAngle DataFiles (%s-%s)',AngleType,...
        iff(FilterFlag,'Filtered','Unfiltered')),TreatmentType);
    if ~isdir(csvfolder)
        mkdir(csvfolder);
    end
    
    % Start loop
    for i=1:length(experiments)
        trials = fieldnames(Dataset.(experiments{i}));
        for j=1:length(trials)
            windspeeds = fieldnames(Dataset.(experiments{i}).(trials{j}));
            curr_ind = curr_ind + 1;
            % set left metadata - means
            mean_lca_means.Moth{curr_ind} = experiments{i};
            mean_lca_means.Trial{curr_ind} = trials{j};
            stder_lca_means.Moth{curr_ind} = experiments{i};
            stder_lca_means.Trial{curr_ind} = trials{j};
            % set right metadata - means
            mean_rca_means.Moth{curr_ind} = experiments{i};
            mean_rca_means.Trial{curr_ind} = trials{j};
            stder_rca_means.Moth{curr_ind} = experiments{i};
            stder_rca_means.Trial{curr_ind} = trials{j};
            % set left metadata - dist
            mean_lca_dist.Moth{curr_ind} = experiments{i};
            mean_lca_dist.Trial{curr_ind} = trials{j};
            stder_lca_dist.Moth{curr_ind} = experiments{i};
            stder_lca_dist.Trial{curr_ind} = trials{j};
            % set right metadata - dist
            mean_rca_dist.Moth{curr_ind} = experiments{i};
            mean_rca_dist.Trial{curr_ind} = trials{j};
            stder_rca_dist.Moth{curr_ind} = experiments{i};
            stder_rca_dist.Trial{curr_ind} = trials{j};
                       
            % create calculated angle matrix
            left_corrected_angles = [];
            right_corrected_angles = [];
            
            for k=1:length(windspeeds)
                % Get antennal angle
                data = Dataset.(experiments{i}).(trials{j}).(windspeeds{k});
                
                if strcmp(AngleType,'HeadVect')
                    [projangle,debugvars] = calculatePerturbedAntennalAngleUsingHeadVect(data,FilterFlag);
                    % Save data
                    projAngles.(experiments{i}).(trials{j}).(windspeeds{k}) = projangle;
                    debugVars.(experiments{i}).(trials{j}).(windspeeds{k}) = debugvars;
                elseif strcmp(AngleType,'HeadVectAP')
                    [projangle,debugvars] = calculatePerturbedAntennalAngleUsingAntennaPlaneHeadVect(data,FilterFlag);
                    % Save data
                    projAngles.(experiments{i}).(trials{j}).(windspeeds{k}) = projangle;
                    debugVars.(experiments{i}).(trials{j}).(windspeeds{k}) = debugvars;
                elseif strcmp(AngleType,'HeadVectHP')
                    [projangle,debugvars] = calculatePerturbedAntennalAngleUsingHeadPlaneHeadVect(data,FilterFlag);
                    % Save data
                    projAngles.(experiments{i}).(trials{j}).(windspeeds{k}) = projangle;
                    debugVars.(experiments{i}).(trials{j}).(windspeeds{k}) = debugvars;
                else  % if perturbed angle
                    if k==1
                        [antangle,projangle,sphangle,refvect,projrefvect,debugvars] = ...
                            calculatePerturbedAntennalAngle(data);
                    else
                        [antangle,projangle,sphangle,refvect,projrefvect,debugvars] = ...
                            calculatePerturbedAntennalAngle(data,refvect,projrefvect);
                    end
                    % Save data
                    antAngles.(experiments{i}).(trials{j}).(windspeeds{k}) = antangle;
                    projAngles.(experiments{i}).(trials{j}).(windspeeds{k}) = projangle;
                    sphAngles.(experiments{i}).(trials{j}).(windspeeds{k}) = sphangle;
                    debugVars.(experiments{i}).(trials{j}).(windspeeds{k}) = debugvars;
                end
                
                % Angle used
                rightprojangle = projangle(:,1);
                leftprojangle = projangle(:,2);
                right_ca = [];
                left_ca = [];
                
                % Take individual means
                mean_left_ca = NaN(length(control_start),1);
                mean_right_ca = NaN(length(control_start),1);
                for l=1:length(control_start)
                    % Get corrected right antennal angle
                    mean_right_ca(l) = ...
                        nanmean(rightprojangle(control_start(l)+mean_window(1):control_start(l)+mean_window(2)));
                    temp = rightprojangle(control_start(l)+mean_window(1):control_start(l)+mean_window(2));
                    % corrected_angle = [corrected_angle;temp(~isnan(temp))];
                    right_ca = [right_ca;temp]; %#ok<*AGROW>
                    
                    mean_left_ca(l) = ...
                        nanmean(leftprojangle(control_start(l)+mean_window(1):control_start(l)+mean_window(2)));
                    temp = leftprojangle(control_start(l)+mean_window(1):control_start(l)+mean_window(2));
                    % corrected_angle = [corrected_angle;temp(~isnan(temp))];
                    left_ca = [left_ca;temp];
                end
                % Remove NaNs and Update corrected angle matrix
                curr_speed = floor(windspeed) == str2double(windspeeds{k}(4));
                right_ca = right_ca(~isnan(right_ca));
                right_corrected_angles = [right_corrected_angles;...
                    [windspeed(curr_speed).*ones(size(right_ca)),right_ca]];
                left_ca = left_ca(~isnan(left_ca));
                left_corrected_angles = [left_corrected_angles;...
                    [windspeed(curr_speed).*ones(size(left_ca)),left_ca]];
                
                % Take mean and stderr for trials - means
                mean_rca_means.(windspeeds{k})(curr_ind) = nanmean(mean_right_ca);
                stder_rca_means.(windspeeds{k})(curr_ind) = nanstd(mean_right_ca)./sqrt(length(control_start));
                mean_lca_means.(windspeeds{k})(curr_ind) = nanmean(mean_left_ca);
                stder_lca_means.(windspeeds{k})(curr_ind) = nanstd(mean_left_ca)./sqrt(length(control_start));
                 % Take mean and stderr for trials - dist
                mean_rca_dist.(windspeeds{k})(curr_ind) = nanmean(right_ca);
                stder_rca_dist.(windspeeds{k})(curr_ind) = nanstd(right_ca)./sqrt(length(right_ca));
                mean_lca_dist.(windspeeds{k})(curr_ind) = nanmean(left_ca);
                stder_lca_dist.(windspeeds{k})(curr_ind) = nanstd(left_ca)./sqrt(length(left_ca));
                
            end
            
            if j==1
                writeCSV(fullfile(csvfolder,sprintf('%s_RCA.csv',experiments{i})),{'windspeed','RCA(HeadVect)'},right_corrected_angles);
                writeCSV(fullfile(csvfolder,sprintf('%s_LCA.csv',experiments{i})),{'windspeed','LCA(HeadVect)'},left_corrected_angles);
            else
                writeCSV(fullfile(csvfolder,sprintf('%s_%s_RCA.csv',experiments{i},trials{j})),{'windspeed','RCA(HeadVect)'},right_corrected_angles);
                writeCSV(fullfile(csvfolder,sprintf('%s_%s_LCA.csv',experiments{i},trials{j})),{'windspeed','LCA(HeadVect)'},left_corrected_angles);
            end
        end
    end
    
    % Remove empty columns - means
    mean_rca_means(curr_ind+1:end,:) = [];   
    stder_rca_means(curr_ind+1:end,:) = [];  
    mean_lca_means(curr_ind+1:end,:) = [];    %#ok<*NASGU>
    stder_lca_means(curr_ind+1:end,:) = [];  
    % Remove empty columns - dist
    mean_rca_dist(curr_ind+1:end,:) = [];   
    stder_rca_dist(curr_ind+1:end,:) = [];  
    mean_lca_dist(curr_ind+1:end,:) = [];   
    stder_lca_dist(curr_ind+1:end,:) = [];  
    
    % Clean up data
    attributes = struct;
    attributes.AngleType = AngleType;
    attributes.TreatmentType = TreatmentType;
    misc = struct;
    misc.control_start = control_start;
    misc.mean_window = mean_window;
    misc.windspeed = windspeed;
    misc.antAngles = antAngles;
    misc.sphAngles = sphAngles;
    misc.debugVars = debugVars;     %#ok<*STRNU>
   
    % Clear unneccessary parameters
    clearvars antangle projangle sphangle sphangle angleUsed curr_speed;
    clearvars i j k l data curr_ind debugvars experiments left_ca right_ca;
    clearvars left_corrected_angles right_corrected_angles leftprojangle rightprojangle;
    clearvars max_num_moths mean_left_ca mean_right_ca temp trials windspeeds csvfolder;
    clearvars AngleType TreatmentType control_start mean_window windspeed antAngles sphAngles debugVars;

    % Get list to save
    currvars = who;
    tosave = currvars(~ismember(currvars, initvars));
    
    save(fullfile(outputDir_pref,...
        sprintf('CalculatedAngles_%s_%s_%s.mat',...
        iff(FilterFlag,'Filtered','Unfiltered'),...
        attributes.AngleType,attributes.TreatmentType)),...
        tosave{:});
    clearvars('-except',initvars{:});    
end

end