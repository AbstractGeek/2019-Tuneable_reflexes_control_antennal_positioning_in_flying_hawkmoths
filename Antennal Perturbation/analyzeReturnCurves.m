function [ModelFits, ModelData, ExpFits] = analyzeReturnCurves(inputDir, outputDir,...
    treatment_name, angle_name, filter_vals)
%
% Extracts antennal return curves, subtracts the wing beat frequency from
% them and saves it for the next step.
%
% Dinesh Natesan
% 28th Feb 2016

% If no arguments, set default inputs
if (nargin < 1)
    inputDir = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Analysis/mat files';
    outputDir = '/home/dinesh/Dropbox/Science/Projects/Antennal-Positioning/Behavior/Perturbation Experiments/Analysis';
%     angle_name = {'HeadVect','HeadVect','HeadVect','HeadVect'};
%     treatment_name = {'Control','JO','Control','JO'};
%     filter_vals = [false, false, true, true];
    
    % temporary (till JO is digitized)
    angle_name = {'HeadVect','HeadVect'};
    treatment_name = {'Control','JO'};
    filter_vals = [false, false];

end

% Set misc default parameters
folder_name =  'returnCurves (ModelFits)';
experiment_exceptions = {};  % Temp. Till they are digitized.

% Plotting related defaults
figure_size = [1 21 29.7 21];       % A4 size
wsNames = {'ws_0_mps'
    'ws_1_5mps'
    'ws_2_5mps'
    'ws_4_mps'};
stepNames = arrayfun(@(x) sprintf('Step-%g',x),[1,2,3,4,5],...
    'UniformOutput',false);
plotData = 1;

% Switch off warning
warning('off','Ident:estimation:transientDataCorrection');  % Switch off warning temporarily.

% Feedback Model related constants
sam_freq = 1000;
InitialState = 'estimate';
modelpzs = [0,0; 1,0; 1,1; 1,1; 2,0; 2,2]; % P, I, PI, II, PID
modelnames = {'P','I','PI', 'PD','II','PID','Exp'};
modelNo = length(modelnames);
modelNamesValue = {'P [X0 Kp]','I [X0 Ki]','PI [X0 Kp Ki]', 'PD [X0 Kp Kd]', 'II [X0 Ki]','PID [X0 X''0 Kp Ki Kd', 'Exp [a b Tau(ms)]'};
l_modelfns = {'l_proportional','l_integrator','l_proportional_integral','l_proportional_differential','l_double_integrator','l_proportional_integral_differential'};
% nl_modelfns = {'nl_proportional','nl_integrator','nl_proportional_integral','nl_proportional_integral_differential'};

% Log file formatting defaults
colgroups = {'<','','','','','>'};
seperator = '|||||';

% Obtain parameters not to be cleared
initvars = who;
initvars = [initvars;'mf';'initvars'];

for mf=1:length(angle_name)
    % Load return curves
    clearvars('-except',initvars{:})
    load(fullfile(inputDir,...
        sprintf('ReturnCurves_%s_%s_%s.mat',...
        iff(filter_vals(mf),'Filtered','Unfiltered'),...
        angle_name{mf},treatment_name{mf})));
    
    %% fit models to return curves    
    % Out directory
    outDir = fullfile(outputDir,folder_name, sprintf('%s (%s-%s)',...
        treatment_name{mf},iff(filter_vals(mf),'Filtered','Unfiltered'),...
        angle_name{mf}));
    if ~isdir(outDir)
        mkdir(outDir);
    end
    
    h = figure('Units','centimeters','Position',figure_size,'Visible','off');
    
    % Create and initiate fit and value org files
    % Fit
    orgfit = fullfile(outDir,sprintf('%s_fit.org',treatment_name{mf}));
    fidfit = fopen(orgfit,'w');
    fprintf(fidfit, '#+OPTIONS: toc:nil num:nil tags:nil \n');
    fprintf(fidfit, '#+STARTUP: align \n');
    fprintf(fidfit, '#+LATEX_CLASS: article \n');
    fprintf(fidfit, '#+LATEX_CLASS_OPTIONS: [a3paper] \n');
    fprintf(fidfit, '#+LATEX_HEADER: \\usepackage[margin=0.1in, landscape]{geometry} \n');
    % Values
    orgvalue = fullfile(outDir,sprintf('%s_value.org',treatment_name{mf}));
    fidvalue = fopen(orgvalue,'w');
    fprintf(fidvalue, '#+OPTIONS: toc:nil num:nil tags:nil \n');
    fprintf(fidvalue, '#+STARTUP: align \n');
    fprintf(fidvalue, '#+LATEX_CLASS: article \n');
    fprintf(fidvalue, '#+LATEX_CLASS_OPTIONS: [a3paper] \n');
    fprintf(fidvalue, '#+LATEX_HEADER: \\usepackage[margin=0.1in, landscape]{geometry} \n');
    
    % create necessary strings
    headerStringFit = strjoin(repmat(modelnames,1,4),' | ');
    columnStringFit = strjoin(repmat(colgroups,1,4),' | ');
    headerStringValue = strjoin(repmat(modelNamesValue,1,2),' | ');
    columnStringValue = strjoin(repmat(colgroups,1,2),' | ');    
    
    % Start!
    experiments = fieldnames(returnCurves);
    for i=1:length(experiments)
        
        % check if the experiment is an exception
        if any(strcmp(experiment_exceptions,experiments{i}))
            continue;
        end
        
        trials = fieldnames(returnCurves.(experiments{i}));
        outdirExp = fullfile(outDir,(experiments{i}));
        if ~isdir(outdirExp)
            mkdir(outdirExp);
        end
        
        % Create Org files
        orgfile = fullfile(outDir,sprintf('%s_log.org',experiments{i}));
        fidfile = fopen(orgfile,'w');
        fprintf(fidfile,'#+OPTIONS: toc:nil num:nil tags:nil \n');
        
        % Enter experiment names into org
        fprintf(fidfit,'* %s\n',strjoin(strsplit(experiments{i},'_'),'-'));
        fprintf(fidvalue,'* %s\n',strjoin(strsplit(experiments{i},'_'),'-'));
        
        
        for j=1:length(trials)
            
            windspeeds = fieldnames(returnCurves.(experiments{i}).(trials{j}));
            % Enter trial details into org file
            fprintf(fidfile,'* %s\n',strjoin(strsplit(trials{j},'_'),'-'));
            fprintf(fidfit,'** %s\n',strjoin(strsplit(trials{j},'_'),'-'));
            fprintf(fidvalue,'** %s\n',strjoin(strsplit(trials{j},'_'),'-'));
            
            % Initialize necessary tables
            % Fit string table
            fitString = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            fitString(:,:) = {seperator};
            
            % R2 string table
            rsquareString = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            rsquareString(:,:) = {seperator};
            
            % R2 string table
            adjustedrsquareString = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            adjustedrsquareString(:,:) = {seperator};
            
            % Akaike information criterion table
            naicString = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            naicString(:,:) = {seperator};
        
            % Constants Value string
            valueString = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            valueString(:,:) = {seperator};
            
            % Table for fit data and models
            fitData = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            fittedModels = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            expGOF = cell2table(cell(length(stepNames),length(wsNames)),...
                'VariableNames',wsNames,'RowNames',stepNames);
            
            % Go through each windspeed and fit each return trajectory
            for k=1:length(windspeeds)
                
                % Method 1: Linear Model
                ws_data = returnCurves.(experiments{i}).(trials{j}).(windspeeds{k});
                
                % Method 2: Non-linear grey box fitting (Not used
                % currently)
                % ws_fitdata = filtReturnCurves.(experiments{i}).(trials{j}).(windspeeds{k});
                
                for l=1:length(ws_data)                 
                    
                    % Linear modeling
                    % Copy data
                    curr_data = ws_data{l,1};
                   
                    % Check if there was a valid perturbation
                    if isempty(curr_data)
                        
                        % No perturbation found
                        fprintf('Step-%d of %s-%s-%s: No perturbation found \n',...
                            l,experiments{i},trials{j},windspeeds{k});                        
                        
                        % Write empty fit values into string tables
                        fitString.(windspeeds{k}){l} = strjoin(...
                            arrayfun(@(x) sprintf('%0.4g',x),NaN(modelNo,1)',...
                            'UniformOutput',false),' | ');
                        valueString.(windspeeds{k}){l} = strjoin(...
                            arrayfun(@(x) sprintf('%0.4g',x),NaN(modelNo,1)',...
                            'UniformOutput',false),' | ');
                        
                        % Write rsquared values and aic values into string table                        
                        rsquareString.(windspeeds{k}){l} = strjoin(...
                            arrayfun(@(x) sprintf('%0.4g',x),NaN(modelNo,1)',...
                            'UniformOutput',false),' | ');
                        adjustedrsquareString.(windspeeds{k}){l} = strjoin(...
                            arrayfun(@(x) sprintf('%0.4g',x),NaN(modelNo,1)',...
                            'UniformOutput',false),' | ');
                        naicString.(windspeeds{k}){l} = strjoin(...
                            arrayfun(@(x) sprintf('%0.4g',x),NaN(modelNo,1)',...
                            'UniformOutput',false),' | ');
                        
                        % Continue the loop
                        continue;
                    end                    
                    
                    % Extract input, output, error and time from the
                    % dataset
                    output_data = curr_data(:,1);
                    input_data = curr_data(:,2);
                    error_data = input_data - output_data;
                    time_data = (0:length(error_data)-1)'./sam_freq;
                    
                    % Initialize necessary cell arrays to store the linear
                    % models and the obtained outputs
                    sys = cell(modelNo,1);
                    rsquare_temp = cell(modelNo,1);
                    naic_temp = cell(modelNo,1);
                    adjust_rsquare_temp = cell(modelNo,1);
                    tempString = cell(modelNo,1);
                    tempString(:) = {''};
                    fprintf(fidfile,'*** Step %g\n',l);
                    
                    % Fit each model to the return curve and obtain extract
                    % the outputs                    
                    for m=1:size(modelpzs,1) 
                        
                        % Linear Model fitting (default inital state =
                        % estimate)
                        [sys{m},data] = feval(l_modelfns{m},input_data,...
                            output_data, sam_freq, InitialState);
                        
                        % Print fit details into the org file
                        fprintf(fidfile,'**** Model: %s\n',modelnames{m});
                        tftext = evalc('sys{m}');
                        fprintf(fidfile,'Model Details: %s',tftext(8:end));
                        
                        % Save fit constants
                        tempString{m} = sprintf('[%s %s]',...
                            num2str(sys{m}.Report.Parameters.X0','%0.2f '),...
                            num2str(sys{m}.Report.Parameters.ParVector','%0.2f '));
                        
                        % Save rsquared value
                        temp_y = compare(data,sys{m});
                        [rsquare_temp{m},adjust_rsquare_temp{m}] = rsquared(...
                            output_data, temp_y.OutputData,...
                            length(sys{m}.Report.Parameters.X0) + ...
                            length(sys{m}.Report.Parameters.ParVector));
                        
                        % Save naic value (Akaike information criterion
                        % table)
                        naic_temp{m} = sys{m}.Report.Fit.nAIC;
                        
                    end
                    
                    % Try an exponential fit of the error
                    [xData, yData] = prepareCurveData( time_data, error_data );
                    % Set up fittype and options.
                    ft = fittype( 'exp1' );
                    opts = fitoptions( 'Method', 'NonlinearLeastSquares');
                    opts.Display = 'Off';
                    opts.StartPoint = [0,0];
                    % Fit model to data.
                    [sys{end}, gof] = fit( xData, yData, ft, opts );
                    
                    tempString{end} = sprintf('[%s %s %s]',...
                        num2str(sys{end}.a,'%0.2f '),...
                        num2str(sys{end}.b,'%0.2f '),...
                        num2str(1000/sys{end}.b,'%0.2f '));
                    rsquare_temp{end} = gof.rsquare;                  
                    adjust_rsquare_temp{end} = gof.adjrsquare;
                    naic_temp{end} = NaN;
                    
                    % Print model fit complete for step
                    fprintf('Step-%d of %s-%s-%s: Model Fit Complete \n',...
                        l,experiments{i},trials{j},windspeeds{k});
                    
                    if plotData
                        % Plot the fitted models                        
                        subplot(2,1,1);
                        compare(data,sys{1},sys{2},sys{3},sys{4},sys{5},sys{6});
                        legend([{'Data'},modelnames(1:end-1)],'Box','off','Location','best');
                        
                        subplot(2,1,2);
                        plot( sys{end}, xData, yData );
                        legend('Error Data', 'Exponential Fit', 'Location', 'best' );
                        
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
                    
                    % Save model outputs and fit values
                    [y_compare,fit_compare,x0_compare] = ...
                        compare(data,sys{1},sys{2},sys{3},sys{4},sys{5},sys{6});                    
                    y_compare{length(y_compare)+1} = sys{end}(time_data);
                    fit_compare{length(y_compare)} = 100 * goodnessOfFit(y_compare{end}, yData, 'NRMSE');
                    x0_compare{length(y_compare)} = [];
                    
                    % Write fit values into string table
                    fitString.(windspeeds{k}){l} = strjoin(...
                        cellfun(@(x) sprintf('%0.4g',x),fit_compare',...
                        'UniformOutput',false),' | ');
                    valueString.(windspeeds{k}){l} = strjoin(tempString,'|');
                    % Write rsquared values and aic values into string table
                    rsquareString.(windspeeds{k}){l} = strjoin(...
                        cellfun(@(x) sprintf('%0.4g',x),rsquare_temp',...
                        'UniformOutput',false),' | ');
                    adjustedrsquareString.(windspeeds{k}){l} = strjoin(...
                       cellfun(@(x) sprintf('%0.4g',x),adjust_rsquare_temp',...
                       'UniformOutput',false),' | ');
                    naicString.(windspeeds{k}){l} = strjoin(...
                        cellfun(@(x) sprintf('%0.4g',x),naic_temp',...
                        'UniformOutput',false),' | ');                    
                    
                    % Save the data
                    fitData.(windspeeds{k}){l} = table(input_data, output_data,...
                        time_data,error_data);
                    fittedModels.(windspeeds{k}){l} = table(sys,y_compare,...
                        fit_compare,x0_compare,rsquare_temp,...
                        adjust_rsquare_temp,naic_temp,...
                        'VariableNames',{'Sys','Y','Fit','X0','R2','adjR2','nAIC'},...
                        'RowNames',modelnames);
                    expGOF.(windspeeds{k}){l} = gof;
                    
                end
                
            end
            
            % Save model fit data into the main structure
            ModelFits.(experiments{i}).(trials{j}) = fittedModels;
            ModelData.(experiments{i}).(trials{j}) = fitData;
            ExpFits.(experiments{i}).(trials{j}) = expGOF;
            
            % Write into the goodness of fit parameters to file
            % Fit percentage (RSE)
            fprintf(fidfit,'Fit Percentage\n');
            WriteTableAsOrgFile(fidfit,fitString,headerStringFit,columnStringFit);
            fprintf(fidfit,'\n');
            % R2
            fprintf(fidfit,'Coefficient of Determination (R2)\n');
            WriteTableAsOrgFile(fidfit,rsquareString,headerStringFit,columnStringFit);
            fprintf(fidfit,'\n');
            % Adjusted R2
            fprintf(fidfit,'Adjusted Coefficient of Determination (adjR2)\n');
            WriteTableAsOrgFile(fidfit,adjustedrsquareString,headerStringFit,columnStringFit);
            fprintf(fidfit,'\n');
            % nAIC
            fprintf(fidfit,'Normalized Akaike Information Criterion\n');
            WriteTableAsOrgFile(fidfit,naicString,headerStringFit,columnStringFit);
            fprintf(fidfit,'\n');
            
            % Write the value file as two tables
            WriteTableAsOrgFile(fidvalue,valueString(:,1:2),headerStringValue,columnStringValue);
            fprintf(fidvalue,'\n');
            WriteTableAsOrgFile(fidvalue,valueString(:,3:4),headerStringValue,columnStringValue);
            fprintf(fidvalue,'\n');
            
        end
        
        % Add a new page command for latex
        fprintf(fidfit,'#+LATEX: \\newpage \n');
        fprintf(fidvalue,'#+LATEX: \\newpage \n');
        
        % Keep saving the mat file (so as not to lose it)
        save(fullfile(inputDir,...
        sprintf('ModelFits_%s_%s_%s.mat',...
        iff(filter_vals(mf),'Filtered','Unfiltered'),...
        angle_name{mf},treatment_name{mf})),...
        'ModelFits','ModelData','ExpFits','perturbStatus');
    
    end
    
    close(h);        
    
end

% Switch on warning again
warning('on','Ident:estimation:transientDataCorrection');

end
