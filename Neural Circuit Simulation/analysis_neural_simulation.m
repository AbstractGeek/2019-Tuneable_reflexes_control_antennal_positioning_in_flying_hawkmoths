%% Defaults
sam_freq = 1000;
wsNames = {'ws_0_mps'
    'ws_1_5mps'
    'ws_2_5mps'
    'ws_4_mps'};
wsLabels = {'0 mps', '1.5 mps', '2.5 mps', '4 mps'};
windspeeds = [0,1.5,2.5,4];
modelnames = {'P','I','PI','II','PID','Exp'};

matfolder = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Simulation/Antennal-motor-system-neural-network/Control-Theoretic-Analysis/mat files';
outfolder = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Simulation/Antennal-motor-system-neural-network/Control-Theoretic-Analysis/analysis';

hfig = figure('Units','centimeters','Position',[2 2 20 15]);
%% Save raw data plots
mat_files =  dir(fullfile(matfolder,'RawData*.mat'));
savefolder = fullfile(outfolder,'plots');

if ~isdir(fullfile(savefolder, 'SimPosIn'))
    mkdir(fullfile(savefolder, 'SimPosIn'));
end

if ~isdir(fullfile(savefolder, 'SimPosEx'))
    mkdir(fullfile(savefolder, 'SimPosEx'));
end

for f=1:size(mat_files, 1)
    clearvars rawData;
    [~, tok] = regexp(mat_files(f).name,...
        'RawData_Unfiltered_(SimPos\w*)_(mKp-[\d \.]*).mat',...
        'match', 'tokens');
    load(fullfile(matfolder, mat_files(f).name));  
    
    experiments = fieldnames(rawData);
    
    for i=1:length(experiments)
        trials = fieldnames(rawData.(experiments{i}));
        for j=1:length(trials)
            hold on;
            for k=1:length(wsNames)
                curr_data = rawData.(experiments{i}).(trials{j}).(wsNames{k});
                plot(curr_data.time, curr_data.position, 'DisplayName', wsLabels{k});                
            end
            % save plot
            filename = sprintf('rawplot_%s_%s_Trial-%d.eps', tok{1}{1}, tok{1}{2}, j);
            print('-depsc2',fullfile(savefolder, tok{1}{1}, filename));
            filename = sprintf('rawplot_%s_%s_Trial-%d.fig', tok{1}{1}, tok{1}{2}, j);
            savefig(hfig, fullfile(savefolder, tok{1}{1}, filename));
            filename = sprintf('rawplot_%s_%s_Trial-%d.png', tok{1}{1}, tok{1}{2}, j);
            saveas(hfig, fullfile(savefolder, tok{1}{1}, filename));
            clf(hfig);        
        end
    end
    
end
%% Perform stats
mat_files =  dir(fullfile(matfolder,'ModelFits*.mat'));
savefolder = fullfile(outfolder,'stats');

if ~isdir(fullfile(savefolder, 'SimPosIn'))
    mkdir(fullfile(savefolder, 'SimPosIn'));
end

if ~isdir(fullfile(savefolder, 'SimPosEx'))
    mkdir(fullfile(savefolder, 'SimPosEx'));
end


for f=1:size(mat_files, 1)
    clearvars ModelFits returnCurves perturbStatus;
    [~, tok] = regexp(mat_files(f).name,...
        'ModelFits_Unfiltered_(SimPos\w*)_(mKp-[\d \.]*).mat',...
        'match', 'tokens');
    load(fullfile(matfolder, mat_files(f).name));    
    
    % adjR2 box plot
    experiments = fieldnames(ModelFits);
    adjr2_values = cell(5,1);
    
    for i=1:length(experiments)
        trials = fieldnames(ModelFits.(experiments{i}));
        for j=1:length(trials)
            for k=1:length(wsNames)
                for l=1:5
                    if isempty(ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l})
                        continue;
                    end
                    
                    adjr2_values{1} = [adjr2_values{1};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.adjR2{1}];
                    adjr2_values{2} = [adjr2_values{2};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.adjR2{2}];
                    adjr2_values{3} = [adjr2_values{3};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.adjR2{3}];
                    adjr2_values{4} = [adjr2_values{4};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.adjR2{4}];
                    adjr2_values{5} = [adjr2_values{5};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.adjR2{5}];
                    
                end
                
            end
        end
    end
    
    adjr2_values = cell2mat(adjr2_values');
    
    CategoricalScatterplot(adjr2_values,modelnames(1:end-1));
    axis square;
    filename = sprintf('adjR2_%s_%s.eps', tok{1}{1}, tok{1}{2});
    print('-depsc2',fullfile(savefolder, tok{1}{1}, filename));
    filename = sprintf('adjR2_%s_%s.fig', tok{1}{1}, tok{1}{2});
    savefig(hfig, fullfile(savefolder, tok{1}{1}, filename));
    clf(hfig);
    
    % nAIC box plot
    nAIC = cell(5,1);
    
    for i=1:length(experiments)
        trials = fieldnames(ModelFits.(experiments{i}));
        for j=1:length(trials)
            for k=1:length(wsNames)
                for l=1:5
                    if isempty(ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l})
                        continue;
                    end
                    
                    nAIC{1} = [nAIC{1};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.nAIC{1}];
                    nAIC{2} = [nAIC{2};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.nAIC{2}];
                    nAIC{3} = [nAIC{3};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.nAIC{3}];
                    nAIC{4} = [nAIC{4};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.nAIC{4}];
                    nAIC{5} = [nAIC{5};ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.nAIC{5}];
                    
                end
                
            end
        end
    end
    
    nAIC = cell2mat(nAIC');
    
    CategoricalScatterplot(nAIC,modelnames(1:end-1));
    axis square;    
    filename = sprintf('nAIC_%s_%s.eps', tok{1}{1}, tok{1}{2});
    print('-depsc2',fullfile(savefolder, tok{1}{1}, filename));
    filename = sprintf('nAIC_%s_%s.fig', tok{1}{1}, tok{1}{2});
    savefig(hfig, fullfile(savefolder, tok{1}{1}, filename));
    clf(hfig);
    
end
close(hfig);