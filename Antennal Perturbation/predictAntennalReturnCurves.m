function []= predictAntennalReturnCurves(inputDir, outputDir,...
    treatment_name, angle_name, filter_vals)

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
folder_name =  'returnCurves (Predicted)';

% Plotting related defaults
figure_size = [1 21 29.7 21];       % A4 size
colors = rgb('red', 'green', 'blue', 'purple', 'brown');
wsNames = {'ws_0_mps'
    'ws_1_5mps'
    'ws_2_5mps'
    'ws_4_mps'};
stepNames = arrayfun(@(x) sprintf('Step-%g',x),[1,2,3,4,5],...
    'UniformOutput',false);
stepNamesPrediction = arrayfun(@(x,y) sprintf('[%0.1f m/s] Step-%g',x ,y),...
    [0.*ones(5,1);1.5.*ones(5,1);2.5.*ones(5,1);4.*ones(5,1)]', ...
    repmat([1,2,3,4,5],1,4), 'UniformOutput',false)';
parameterValNames = arrayfun(@(x) sprintf('Parameters_%s',x{1}(4:end)),...
    wsNames, 'UniformOutput',false)';

plotData = 1;

% Switch off warning
warning('off','Ident:estimation:transientDataCorrection');  % Switch off warning temporarily.

% Feedback Model related constants
sam_freq = 1000;
InitialState = 'Estimate';
l_modelfns = {'l_proportional','l_integrator','l_proportional_integral', 'l_proportional_differential','l_double_integrator','l_proportional_integral_differential'};
modelnames = {'P','I','PI','PD','II','PID'};

%Create figure
h = figure('Units','centimeters','Position',figure_size,'Visible','off');
hold on;

% Obtain parameters not to be cleared
initvars = who;
initvars = [initvars;'mf';'initvars'];

for mf=1:length(angle_name)
    
    % Load return curves
    clearvars('-except',initvars{:})
    load(fullfile(inputDir,...
        sprintf('ModelFits_%s_%s_%s.mat',...
        iff(filter_vals(mf),'Filtered','Unfiltered'),...
        angle_name{mf},treatment_name{mf})));
    
    % Out directory
    outDir = fullfile(outputDir,folder_name, sprintf('%s (%s-%s)',...
        treatment_name{mf},...
        iff(filter_vals(mf),'Filtered','Unfiltered'),angle_name{mf}));
    if ~isdir(outDir)
        mkdir(outDir);
    end
    
    % Obtain experiments
    experiments = fieldnames(ModelFits);        
    
    % Parameter initialization
    modelparameters = struct;
    modelparametercovar = struct;
    predicted_models = struct;
    predicted_fit = struct;
    predicted_R2 = struct;
    predicted_adjR2 = struct;
    predicted_nAIC = struct;
    
    for i=1:length(experiments)
        
        % Get trials
        trials = fieldnames(ModelFits.(experiments{i}));
        
        % Create Experiment directories
        outdirExp = fullfile(outDir,(experiments{i}));
        if ~isdir(outdirExp)
            mkdir(outdirExp);
        end
        
        for j=1:length(trials)
            
            for ml = 1:length(l_modelfns)
                
                % Extract model constants
                modelparameters.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(cell(6,1),cell(6,1),cell(6,1),cell(6,1),...
                    'VariableNames',wsNames,'RowNames',[stepNames,{'Step-Mean'}]);
                modelparametercovar.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(cell(5,1),cell(5,1),cell(5,1),cell(5,1),...
                    'VariableNames',wsNames,'RowNames',stepNames);
                % Made tables for predicted output
                predicted_models.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(cell(20,1),cell(20,1),cell(20,1),cell(20,1),...
                    'VariableNames',parameterValNames,'RowNames',stepNamesPrediction);
                predicted_fit.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(NaN(20,1),NaN(20,1),NaN(20,1),NaN(20,1),...
                    'VariableNames',parameterValNames,'RowNames',stepNamesPrediction);
                predicted_R2.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(NaN(20,1),NaN(20,1),NaN(20,1),NaN(20,1),...
                    'VariableNames',parameterValNames,'RowNames',stepNamesPrediction);
                predicted_adjR2.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(NaN(20,1),NaN(20,1),NaN(20,1),NaN(20,1),...
                    'VariableNames',parameterValNames,'RowNames',stepNamesPrediction);
                predicted_nAIC.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}) = ...
                    table(NaN(20,1),NaN(20,1),NaN(20,1),NaN(20,1),...
                    'VariableNames',parameterValNames,'RowNames',stepNamesPrediction);   
                
                % Run through windspeeds
                for k=1:length(wsNames)
                    
                    % If the windspeed is empty, skip the loop
                    if sum(perturbStatus.(experiments{i}).(trials{j}).(wsNames{k})) == 0
                        fprintf('[%s] %s-%s-%s: No perturbation found (0 out of 5). Skipping model prediction for %s \n',...
                            modelnames{ml},experiments{i},trials{j},wsNames{k},wsNames{k});
                        continue;
                    end
                    
                    % Extract model parameters for the prediction step
                    for l=1:length(stepNames)
                        if perturbStatus.(experiments{i}).(trials{j}).(wsNames{k})(l)
                            current_model = ModelFits.(experiments{i}).(trials{j}).(wsNames{k}){l}.Sys{ml};
                            modelparameters.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(wsNames{k}){l} = current_model.Report.Parameters.ParVector;
                            modelparametercovar.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(wsNames{k}){l} = current_model.Report.Parameters.FreeParCovariance;
                        end
                    end
                    
                    % Get the mean for this airflow
                    parameter_mean = nanmean(cell2mat(modelparameters.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(wsNames{k})'),2);
                    modelparameters.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(wsNames{k}){end} = parameter_mean;
                    
                    % Predict all antennal curves using this parameter mean
                    % of each linear model.
                    for m = 1:length(wsNames)
                                                
                        for l=1:length(stepNames)
                            
                            % If no perturbation, print and continue.
                            if ~perturbStatus.(experiments{i}).(trials{j}).(wsNames{m})(l)
                                continue;
                            end
                            
                            % Obtain necessary time series data
                            time_data = ModelData.(experiments{i}).(trials{j}).(wsNames{m}){l}.time_data;
                            input_data = ModelData.(experiments{i}).(trials{j}).(wsNames{m}){l}.input_data;
                            output_data = ModelData.(experiments{i}).(trials{j}).(wsNames{m}){l}.output_data;
                            
                            % Try to predict the return curve using the
                            % parameter mean for current airflow (outer
                            % loop)
                            [lgr,data] = feval(l_modelfns{ml},input_data, output_data,...
                                sam_freq, InitialState, parameter_mean);
                            [model_output,fit_percent,X0] = compare(data,lgr);
                            model_output = model_output.OutputData;
                            
                            % Sort and save model output into necessary tables                            
                            predicted_models.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(parameterValNames{k}){l+(m-1)*5} = ...
                                struct('Fit_Model', lgr, 'Model_Output', model_output,...
                                'Fit_Percent', fit_percent, 'InitialState', X0);
                            % Fit percentage (RSE)
                            predicted_fit.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(parameterValNames{k})(l+(m-1)*5) = fit_percent;
                            % R2 and AdjR2. The independent parameters are
                            % as before - number of parameters used to fit
                            % the dataset. (even though the main parameter
                            % vector is constant, the presence of a fixed
                            % parameter will sometimes increase R2).
                            [predicted_R2.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(parameterValNames{k})(l+(m-1)*5),...
                                predicted_adjR2.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(parameterValNames{k})(l+(m-1)*5)]= ...
                                rsquared(output_data, model_output, ...
                                length(X0)+length(parameter_mean));
                            predicted_nAIC.(l_modelfns{ml}(3:end)).(experiments{i}).(trials{j}).(parameterValNames{k})(l+(m-1)*5) = ...
                                lgr.Report.Fit.nAIC;
                            
                            % Plot data
                            if plotData
                                plot(time_data, output_data, 'Color', colors(l,:), ...
                                    'DisplayName', sprintf('Step-%d',l));
                                plot(time_data, model_output, 'Color', colors(l,:), ...
                                    'LineStyle', '--', 'DisplayName',...
                                    sprintf('Fit-%d(%0.2f)',l,fit_percent));
                            end
                            
                        end
                        
                        % Print status
                        % Print model fit complete for step
                        fprintf('[%s] %s based model prediction complete for %s-%s \n',...
                            modelnames{ml},wsNames{k},experiments{i},trials{j});
                        
                        if plotData
                            % Add title
                            figtitle(strjoin(strsplit(sprintf('%s-%s:[%s] Predicted %s using parameters from %s',...
                                experiments{i},trials{j},modelnames{ml},wsNames{m},wsNames{k},l),'_'),'-'),...
                                'fontweight','bold');
                            
                            % Switch on Legend
                            legend show;
                            
                            % save into plots
                            export_fig(fullfile(outdirExp,...
                                sprintf('%s_%s_Predicted_%s_Parameters_%s',trials{j},modelnames{ml},wsNames{m},wsNames{k})),...
                                '-transparent',h);
                            export_fig(fullfile(outDir,sprintf('%s',experiments{i})),...
                                '-pdf','-append','-transparent',h);
                            
                            clf(h);
                            hold on;
                        end
                    end
                end
            end
        end
        
        
        save(fullfile(inputDir,...
            sprintf('ModelPredictions_all_%s_%s_%s.mat',...
            iff(filter_vals(mf),'Filtered','Unfiltered'),...
            angle_name{mf},treatment_name{mf})),...
            'modelparameters','modelparametercovar','predicted_models',...
            'predicted_fit','predicted_R2','predicted_adjR2','predicted_nAIC');
        % clearvars -except Ki_Value Ki_Covar Predicted_Models Predicted_Fit Predicted_R2;
        
    end
end

% Switch on warning again
warning('on','Ident:estimation:transientDataCorrection');

end