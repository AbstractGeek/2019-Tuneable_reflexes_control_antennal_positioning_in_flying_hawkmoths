function [] = calculateIAfromData()
% function [] = calculateIAfromData()
% Scans through the input folder (which is the output folder of the
% copyAntennaExperiment.m file), pulls in experiment data from the mat files,
% calculate the interantennal angles, plot them and save the processed data
% as a mat file into the same input directory.
%
% Convention
% Right antennae, Left antennae, Right antennal base, Left antennal base
%
% Dinesh Natesan, 12-09-2014
% Modified, 9-10-2014
%

% User Preferences
inputDir_pref = pwd;
% inputDir_pref = '/home/dinesh/Dropbox/Ncbs/Moth Projects/Johnston''s organ experiments/Experimental DataFiles';
debug = false;
plotData = true;

% Constants
wsNames = {'ws_0_mps'
    'ws_0_5mps'
    'ws_1_mps'
    'ws_1_5mps'
    'ws_2_mps'
    'ws_2_5mps'
    'ws_3_mps'
    'ws_3_5mps'
    'ws_4_mps'
    'ws_4_5mps'
    'ws_5_mps'};

% Treatments
groupNames = {{'control','normal'},'sham','_jo_'};
trtNames = {'control','sham','jo'};

% Ensure Names of the variables are same! And the order in the cell array
% are correct!
plotGroups = {{'raa','laa','iaa'}};

% Let's rock and roll?
% Get the experiment folder
expDir = uigetdir(inputDir_pref,'Experimental DataFiles');
resDir = fullfile(expDir,'CurrentPlots');
if ~isdir(resDir)
    mkdir(resDir);
end
% Obtain a list of subfolders
listOfContents = dir(expDir);
isubDir = [listOfContents(:).isdir];
dirList = {listOfContents(isubDir).name}';
dirList(ismember(dirList,{'.','..','ArchivedPlots','CurrentPlots'})) = [];
dirNum = length(dirList);
% Done

% Open a waitbar
h = waitbar((1-1)/(dirNum-1),strcat('Processing: Trial ',num2str(1),'-',num2str(dirNum),'...'));
pos = get(h,'Position');
set(h,'Position',[pos(1) pos(2)+2.*pos(4) pos(3) pos(4)]);

% Go into each folder, load matfile and calculate ia.
for i=1:dirNum
    % Update waitbar
    waitbar((i-1)/(dirNum-1),h,strcat('Processing: Trial ',num2str(i),'-',num2str(dirNum),'...'));
    % Obtain list of mat files
    currDir = fullfile(expDir,dirList{i});
    matFiles = dir(fullfile(currDir,'*.mat'));
    % Just keep the names. Ignore rest.
    matFiles = {matFiles(:).name}';
    
    % Begin Loading matData
    for j=1:length(matFiles)
        matData = load(fullfile(currDir,matFiles{j}));
        digitizationError = NaN(11,2);      % Left Std, Right Std
        currTrt = '';
        % Find which group the file belongs to
        for k=1:length(groupNames)
            if iscellstr(groupNames{k})
                for l=1:length(groupNames{k})
                    if ~isempty(strfind(matFiles{j},groupNames{k}{l})) && isempty(currTrt)
                        currTrt = trtNames{k};
                    elseif ~isempty(strfind(matFiles{j},groupNames{k}{l})) && ~isempty(currTrt)
                        error(strcat(matFiles{j},':Contains more than one treatment name!'));
                    end
                    
                end
            else
                if ~isempty(strfind(matFiles{j},groupNames{k})) && isempty(currTrt)
                    currTrt = trtNames{k};
                elseif ~isempty(strfind(matFiles{j},groupNames{k})) && ~isempty(currTrt)
                    error(strcat(matFiles{j},':Contains more than one treatment name!'));
                end
            end
            
        end
        
        % Calculate IA
        for k=1:11
            M = matData.(wsNames{k});
            % Get Digitization Error
            nanIndexes = isnan(M(:,1));
            antennaLength = [matrixNorm(M(~nanIndexes,1:3)-M(~nanIndexes,7:9)),matrixNorm(M(~nanIndexes,4:6)-M(~nanIndexes,10:12))];
            digitizationError(k,:) = std(antennaLength)./mean(antennaLength);
            % Calculate IA:  Digitization format
            % Right antennae, Left antennae, Right antennal base, Left
            % antennal base
            % antAngles = calculateAntennalAngleUsingCustomPlane(M);
            antAngles = calculateAntennalAngleUsingGlobalCoords(M);
            angleGroups = fieldnames(antAngles);
            
            % Initialize values of k=1;
            if k==1
                for l=1:length(angleGroups)
                    angleNames = fieldnames(antAngles.(angleGroups{l}));
                    nans = cell(length(angleNames),1);
                    nans(:) = {NaN(11,3)};
                    meanAntAngles.(angleGroups{l}) = cell2struct(nans,angleNames);
                    clear nans;
                end
            end
            
            % Set mean values for each of the variables
            for l=1:length(angleGroups)
                angleNames = fieldnames(antAngles.(angleGroups{l}));
                for m=1:length(angleNames)
                    meanAntAngles.(angleGroups{l}).(angleNames{m})(k,1) = (k-1)/2;
                    meanAntAngles.(angleGroups{l}).(angleNames{m})(k,2) = nanmean(antAngles.(angleGroups{l}).(angleNames{m}));
                    meanAntAngles.(angleGroups{l}).(angleNames{m})(k,3) = nanstd(antAngles.(angleGroups{l}).(angleNames{m}));
                end
            end
            
            % Save raw data
            rawData.(currTrt).(matFiles{j}(1:end-4)).(char(wsNames(k))) = antAngles;
            rawData.(currTrt).(matFiles{j}(1:end-4)).(char(wsNames(k))).digitizedPoints = M;
            clear M antAngles angleNames;
        end
        % Save iaData
        iaData.(currTrt).(matFiles{j}(1:end-4)) = meanAntAngles;
        angleGroups = fieldnames(meanAntAngles);
        
        % Normalize Data
        for l=1:length(angleGroups)
            angleNames = fieldnames(meanAntAngles.(angleGroups{l}));
            for m=1:length(angleNames)
                zerodAngle = meanAntAngles.(angleGroups{l}).(angleNames{m})(:,2)-min(meanAntAngles.(angleGroups{l}).(angleNames{m})(:,2));
                normAntAngles.(angleGroups{l}).(angleNames{m}) = [meanAntAngles.(angleGroups{l}).(angleNames{m})(:,1),...
                    [zerodAngle,meanAntAngles.(angleGroups{l}).(angleNames{m})(:,3)]./max(zerodAngle)];
            end
        end
        
        % Save normData
        normData.(currTrt).(matFiles{j}(1:end-4)) = normAntAngles;
        
        % Save digitizationData
        errorData.(matFiles{j}(1:end-4)) = digitizationError;
        
        clear matData zeroedAngle meanAntAngles normAntAngles;
    end
end
save(fullfile(resDir,'InterAntennalAngleData.mat'),'iaData','normData','rawData','errorData');
angData.iaData = iaData;
angData.normData = normData;
close(h);
clearvars -except angData debug expDir resDir plotGroups plotData;

% Saving data done - begin plotting
if plotData
    dataStructs = fieldnames(angData);
    for i=1:length(dataStructs)
        treatments = fieldnames(angData.(dataStructs{i}));
        for j=1:length(treatments)
            % Begin plotting every sample
            exptNames = fieldnames(angData.(dataStructs{i}).(treatments{j}));
            colors = distinguishable_colors(length(exptNames));
            [~,ind] = max(structfun(@(x)(length(fieldnames(x))),angData.(dataStructs{i}).(treatments{j})));
            angleGroups = fieldnames(angData.(dataStructs{i}).(treatments{j}).(exptNames{ind}));
            for k=1:length(angleGroups)
                angleNames = fieldnames(angData.(dataStructs{i}).(treatments{j}).(exptNames{ind}).(angleGroups{k}));
                legendNames = {};
                % Each antennal angle has its own plot
                h1=figure('Position',[1 1 900 1080]);
                                subplot(3,1,1),title(strcat(angleNames{1},'plot'));
                hold on;
%                 h2=figure('Position',[1 1 1280 720]);
                                subplot(3,1,2),title(strcat(angleNames{2},'plot'));
                hold on;
%                 h3=figure('Position',[1 1 1280 720]);
                                subplot(3,1,3),title(strcat(angleNames{3},'plot'));
                hold on;
                % Start loop
                for l=1:length(exptNames)
                    if ~isfield(angData.(dataStructs{i}).(treatments{j}).(exptNames{l}),(angleGroups{k}))                                            
                        continue;
                    end
                    legendNames = {legendNames{1:end},exptNames{l}};
                    subplot(3,1,1),errorbar(angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{1})(:,1),...
                        angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{1})(:,2),...
                        angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{1})(:,3),...
                        'color',colors(l,:),'LineStyle','-');
                    
                    subplot(3,1,2),errorbar(angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{2})(:,1),...
                        angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{2})(:,2),...
                        angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{2})(:,3),...
                        'color',colors(l,:),'LineStyle','-');
                    subplot(3,1,3),errorbar(angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{3})(:,1),...
                        angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{3})(:,2),...
                        angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{3})(:,3),...
                        'color',colors(l,:),'LineStyle','-');
                    %                     if debug
                    %                         text(angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{3})(end,1),...
                    %                             angData.(dataStructs{i}).(treatments{j}).(exptNames{l}).(angleGroups{k}).(angleNames{3})(end,2),strcat('   ',exptNames{l}));
                    %                     end
                end
                
%                 figure(h1),axis tight;
%                 figure(h2),axis tight;
%                 figure(h3),axis tight;               
                
                subplot(3,1,1),axis tight;
                subplot(3,1,2),axis tight;
                subplot(3,1,3),axis tight;      
                
%                 figure(h1),legend(exptNames),legend('hide');
%                 figure(h2),legend(exptNames),legend('hide');
%                 figure(h3),legend(exptNames),legend('hide');
                
                subplot(3,1,1),legend(legendNames),legend('hide');
                subplot(3,1,2),legend(legendNames),legend('hide');
                subplot(3,1,3),legend(legendNames),legend('hide');
                
                saveDir = fullfile(resDir,treatments{j});
                if ~isdir(saveDir)
                    mkdir(saveDir)
                end
                
                saveas(h1,fullfile(saveDir,strcat('AntennalAngles_',treatments{j},'_AngleGroup',angleGroups{k},'_',dataStructs{i},'.fig')));
                close(h1);
                
%                 saveas(h1,fullfile(saveDir,strcat('AntennalAngles_',treatments{j},'_AngleGroup',angleGroups{k},'_',angleNames{1},'_',dataStructs{i},'.fig')));
%                 close(h1);
%                 saveas(h2,fullfile(saveDir,strcat('AntennalAngles_',treatments{j},'_AngleGroup',angleGroups{k},'_',angleNames{2},'_',dataStructs{i},'.fig')));
%                 close(h2);
%                 saveas(h3,fullfile(saveDir,strcat('AntennalAngles_',treatments{j},'_AngleGroup',angleGroups{k},'_',angleNames{3},'_',dataStructs{i},'.fig')));
%                 close(h3);
                
            end
            
        end
    end
    close all;
end
% Bye!
end
