function performStatisticalAnalysisOnAntennalData
% Get the mat file that has all the Antennal data
[FileName,PathName,FilterIndex] = uigetfile('*.mat','Select Interantennal Data file (mat file)');
% if Cancel button is pressed
if (FilterIndex == 0)
    disp('Action cancelled, quiting');
    return;
end
% Enter and load mat files
cd (PathName);
load (FileName);

% Get the treatments to compare and analyze
treatmentNames = fieldnames(iaData);
[selection,ok] = listdlg('ListString',treatmentNames,'SelectionMode','multiple','Name','Select Treatment','PromptString','Select Treatments to Analyze');
% selection = [] or if Cancel button is pressed
if (ok == 0)
    disp('No Selection / Action cancelled, quiting');
    return;
end

% Sort RawData files first.
for i=1:length(selection)
    w = fieldnames(iaData.(treatmentNames{selection(i)}));
    for k=1:length(w)
        q = fieldnames(rawData.(w{k}));
        iaa = [];
        windspeeds = (0:0.5:5)';
        for j=1:length(q)
            temp = rawData.(w{k}).(q{j}).customplane.iaa;
            iaa = [iaa;[windspeeds(j).*ones(sum(~isnan(temp)),1),temp(~isnan(temp)),j.*ones(sum(~isnan(temp)),1)]];
        end
        normiaa = [iaa(:,1),normalizeVector(iaa(:,2)),iaa(:,3)];
        sortedRawData.(treatmentNames{selection(i)}).(w{k}) = iaa;
        sortedNormData.(treatmentNames{selection(i)}).(w{k}) = normiaa;
    end
end
clearvars w q iaa temp windspeeds;
save('sortedRawData.mat','sortedRawData','sortedNormData');

% Sorting done. Begin statistical analysis
treatmentNames = fieldnames(sortedRawData);
windspeeds = cellstr(num2str((0:0.5:5)'));
% Output Folder
sigSavePath = fullfile(PathName,'BoxPlots');
if ~isdir(sigSavePath)
    mkdir(sigSavePath);
end
insigSavePath = fullfile(sigSavePath,'TheInsignificants');
if ~isdir(insigSavePath)
    mkdir(insigSavePath);
end

for i=1:length(treatmentNames)
    w = fieldnames(sortedRawData.(treatmentNames{i}));
    sigFolder = fullfile(sigSavePath,treatmentNames{i});
    if ~isdir(sigFolder)
        mkdir(sigFolder);
    end
    insigFolder = fullfile(insigSavePath,treatmentNames{i});
    if ~isdir(insigFolder)
        mkdir(insigFolder);
    end
    
    for j=1:length(w)
        % Perform Kruskal Wallis test and 
        statsData.(treatmentNames{i}).(w{j}).X = sortedRawData.(treatmentNames{i}).(w{j})(:,2);
        statsData.(treatmentNames{i}).(w{j}).Group = sortedRawData.(treatmentNames{i}).(w{j})(:,3);
        [statsData.(treatmentNames{i}).(w{j}).h,statsData.(treatmentNames{i}).(w{j}).P,...
            statsData.(treatmentNames{i}).(w{j}).stats,statsData.(treatmentNames{i}).(w{j}).stats_all]...
            = PerformStats(sortedRawData.(treatmentNames{i}).(w{j})(:,2),sortedRawData.(treatmentNames{i}).(w{j})(:,3));
        % Check results and save plots accordingly
        if strcmp(statsData.(treatmentNames{i}).(w{j}).stats.table{2,7},'Reject H0')
            fprintf('%s: %s \n',treatmentNames{i},w{j});
            boxplot(statsData.(treatmentNames{i}).(w{j}).X,statsData.(treatmentNames{i}).(w{j}).Group,'labels',windspeeds);
            ylabel('Interantennal angle (deg)');
            saveas(gcf,fullfile(sigFolder,strcat(w{j},'_BoxPlot')),'fig');
            close(gcf);
        else
            boxplot(statsData.(treatmentNames{i}).(w{j}).X,statsData.(treatmentNames{i}).(w{j}).Group,'labels',windspeeds);
            ylabel('Interantennal angle (deg)');
            saveas(gcf,fullfile(insigFolder,strcat(w{j},'_BoxPlot')),'fig');
            close(gcf);            
        end        
    end    
end

save('statsData.mat','statsData');

end