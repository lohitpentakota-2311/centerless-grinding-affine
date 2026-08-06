% % % % % % % % % % % % % %
% % Regularized Matrix Identification (05/06/2025)
% % % angular domain identification (26/06/2025) not updaTESD 
% % the code should work for angular and freq domain (05/08/2025)

function [bestlambda, bestTransMatrix, Sigma_Y, TestDataSet] = ridge_cv(Input, Output,AngularDomain_u,AngularDomain_y, lambdas, k, seed,Roundness,SimData,InpParV,LFData,Save,Visibility,AngularDomain, ID)

rng(seed)

NsamplePoints=3600/2;

% % % Feed term vector
Input= [Input;ones(1,size(Input,2))];

% Ensure column orientation
[n_features, ~] = size(Input);

[~, ~] = size(Output);

% First: split off for seperating test set
% cv1 = cvpartition(numSamples, 'HoldOut', 0.2);

% trainIdx1 = training(cv1);   % Logical indices for training & validation

trainIdx1 = ID.Train; 

% % % % Roundness
TrainRndMeas= Roundness.MeasFinalRnd(trainIdx1);

% testIdx= test(cv1); % Logical indices for testing

testIdx= ID.Test; 

% % % % Roundness
TestRndMeas= Roundness.MeasFinalRnd(testIdx);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if isequal(LFData,[])
    LFData.IterFinalProfile=[];
    LFData.MatrFinalProfile=[];
    LFData.IterFFT=[];
    LFData.LFMRR= [];
else
    if isequal(LFData.IterFinalProfile,[])
        LFData.IterFinalProfile=[];
        LFData.IterFFT=[];
    else
    LFData.IterFinalProfile=LFData.IterFinalProfile(idx);
    LFData.IterFFT=LFData.IterFFT(:,idx);
    end
    if isequal(LFData.MatrFinalProfile,[])
        LFData.MatrFinalProfile=[];
        LFData.MatrFFT=[];
    else
        % TrainLFData.MatrFinalProfile=LFData.MatrFinalProfile(trainIdx1);
        TrainLFData.MatrFinalProfile= LFData.MatrFinalProfile(1:numel(find(trainIdx1)));
        % TrainLFData.MatrFFT=LFData.MatrFFT(:,trainIdx1);
        TrainLFData.MatrFFT= LFData.MatrFFT(:,1:numel(find(trainIdx1)));
        % % % % % 
        % TestLFData.MatrFinalProfile=LFData.MatrFinalProfile(testIdx);
        TestLFData.MatrFinalProfile= LFData.MatrFinalProfile(numel(find(trainIdx1))+1:end);
        TestLFData.MatrFFT=LFData.MatrFFT(:,numel(find(trainIdx1))+1:end);
    end
    % % High fidelity section
    TrainSimData= SimData(1:numel(find(trainIdx1))); 
    
    TestSimData= SimData(numel(find(trainIdx1))+1:end);
    
    TrainLFData.LFMRR=  LFData.LFMRR(1:numel(find(trainIdx1)));   
    
    TestLFData.LFMRR=  LFData.LFMRR(numel(find(trainIdx1))+1:end);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % rr
% % splits for training and validation
xTrain_Valid = Input(:,trainIdx1);

yTrain_Valid= Output(:,trainIdx1);

% % % % 
FullTestDataSet.xTrain = AngularDomain_u(:,trainIdx1);

FullTestDataSet.yTrain = AngularDomain_y(:,trainIdx1);

% % % For final testing data set
TestDataSet.xTest = Input(:,testIdx);

TestDataSet.yTest = Output(:,testIdx);

% % % % 
FullTestDataSet.xTest = AngularDomain_u(:,testIdx);

FullTestDataSet.yTest = AngularDomain_y(:,testIdx);


cv2 = cvpartition(size(xTrain_Valid ,2), 'KFold', k);

errors = zeros(length(lambdas), 1);            % To store validation errors

errors_Train= zeros(length(lambdas), 1);

for li = 1:length(lambdas)

    lambda = lambdas(li);

    mse_total = 0; rmse_Traintotal= 0;

    for fold = 1:k

        trainIdx = training(cv2, fold);

        validIdx= test(cv2, fold);

        X_train = xTrain_Valid(:, trainIdx);

        Y_train = yTrain_Valid(:, trainIdx);

        % Y_train= [Y_train;ones(1,size(yTrain_Valid,2))];

        X_val   = xTrain_Valid(:, validIdx);

        Y_val   = yTrain_Valid(:, validIdx);

        % Ridge regression formula:
        % TransMatrix = Y_train * X_train' * inv(X_train * X_train' + lambda * I)
        D = blkdiag(eye(n_features-1), 0); 

        TransMatrix = Y_train * X_train' / (X_train * X_train' + lambda * D);

        % Predict on training set
        Y_pred_Train = TransMatrix * X_train;


        % Predict on validation set
        Y_pred = TransMatrix * X_val;

        % Compute mean squared error
        % mse = mean((Y_val(:) - Y_pred(:)).^2);

        % mse= compute_rmse(Y_val,Y_pred, false);

        mse= rmse(Y_val,Y_pred, 'all');

        % % rmse_Train= compute_rmse(Y_train,Y_pred_Train, false);

        rmse_Train= rmse(Y_train,Y_pred_Train, 'all');

        mse_total = mse_total + mse;

        rmse_Traintotal = rmse_Traintotal + rmse_Train;

        % % % % % 
        errors_per_fold(fold, li) = mse;
    end

    errors(li) = mse_total / k;  % Average error across folds

    errors_Train(li) = rmse_Traintotal / k;  % Average error across folds
end

% Find best lambda
% [~, bestIdx] = min(errors);

% % bestLambda = lambdas(bestIdx);
% % % Model selection 
meanValid = errors;                  % 1 x M
meanTrain = errors_Train;            % 1 x M

foldErrorsValid = errors_per_fold;   % K x M

seValid = std(foldErrorsValid, 0, 1) ./ sqrt(size(foldErrorsValid,1));

relDiff = abs(meanTrain - meanValid) ./ max(meanValid, eps);

% === Step 1: best validation lambda ===
[bestVal, bestIdx] = min(meanValid);

oneSE_threshold = bestVal + seValid(bestIdx);

% === Step 2: candidates within 1-SE === rule
candidates = find(meanValid <= oneSE_threshold);

% === Step 3: pick best candidate ===
[~, idxRel] = min(relDiff(candidates));

bestIdxFinal = candidates(idxRel);

% If multiple, pick largest lambda
sameRel = find(abs(relDiff(candidates) - relDiff(candidates(idxRel))) < 1e-8);
if numel(sameRel) > 1
    [~, maxLambdaIdx] = max(lambdas(candidates(sameRel)));
    bestIdxFinal = candidates(sameRel(maxLambdaIdx));
end

bestlambda = lambdas(bestIdxFinal);

fprintf('✅ Best lambda (1-SE + rel diff): λ = %.4f\n', bestlambda);

% % % % % % % % % logic
% Initial threshold setup
% % % % initialThreshold = 0.00001;
% % % % 
% % % % thresholdStep = 0.001;
% % % % 
% % % % maxThreshold = 4;
% % % % 
% % % % threshold = initialThreshold;
% % % % 
% % % % meanTrain = errors_Train;
% % % % 
% % % % meanValid = errors;
% % % % 
% % % % meanDiff = abs(meanTrain - meanValid);
% % % % 
% % % % % Step 2: Sort models by lowest mean validation RMSE
% % % % [valSorted, sortedIdx] = sort(meanValid);  % sorted values and their linear index
% % % % 
% % % % % Initialize outputs
% % % % bestlambda = NaN; % Default value if no model meets the condition
% % % % 
% % % % found = false;
% % % % 
% % % % while ~found && threshold <= maxThreshold
% % % % 
% % % %     for i = 1:length(valSorted)
% % % % 
% % % %         idx = sortedIdx(i);
% % % % 
% % % %         if meanDiff(idx) <= threshold
% % % % 
% % % %             bestlambda =  lambdas(idx);
% % % % 
% % % %             % Step 3: Get the 10-fold RMSE values for this (V1, V2)
% % % % 
% % % %             % trainFolds = squeeze(foldErrorsTrain(r, c, :));  % 10×1
% % % %             %
% % % %             % validFolds = squeeze(foldErrorsValid(r, c, :));  % 10×1
% % % %             %
% % % %             % % Step 4: Find the fold with the smallest absolute difference
% % % %             % [~, bestFold] = min(abs(trainFolds - validFolds));  % returns index 1–10
% % % %             %
% % % %             found = true;
% % % %             break;
% % % %         end
% % % %     end
% % % %     if ~found
% % % %         threshold = threshold + thresholdStep;
% % % %     end
% % % % end
% % % % 
% % % % % Display results
% % % % if isnan(bestlambda)
% % % %     warning('⚠️ No suitable model found even with max threshold of %.2f.', maxThreshold);
% % % % else
% % % %     fprintf('✅ Best lambda found: λ = %.4f\n', bestlambda);
% % % % end

% [min_error, idx_min] = min(errors);
% 
% lambda_min_error = lambdas(idx_min);
% 
% errors_StdDev = std(errors_per_fold, 0, 1);
% % Get the standard deviation of errors at the minimum
% se_at_min = errors_StdDev(idx_min);
% 
% % Calculate the target error for the 1-SE rule
% target_error_1se = min_error + se_at_min;
% 
% idx_1se = find(errors <= target_error_1se, 1, 'last');
% 
% % Lambda selected by 1-SE rule
% lambda_1se = lambdas(idx_1se);
% 
% bestlambda= lambda_1se;

figure;
% loglog(lambdas, errors, '-o',LineWidth=1.5);
errorbar(lambdas, meanValid, seValid','-o',...
    'LineWidth',1.8,'MarkerSize',6,'CapSize',8,...
    'DisplayName','$\mathrm{Validation\ RMSE}$');
hold on;
plot(lambdas, errors_Train, '-o',LineWidth=1.5);
% Plot vertical dashed line at best lambda
xline(bestlambda, '--r', 'LineWidth', 1.5, ...
    'Label', sprintf('Best \\lambda = %.2e', bestlambda), ...
    'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
% Final plot settings
xlabel('$\eta$ (high to low)', 'Interpreter', 'latex', 'FontSize', 16); %
ylabel('RMSE', 'Interpreter', 'latex',FontSize=16);
title('Ridge Parameter Selection via Cross-Validation', 'Interpreter', 'latex',FontSize=16);
set(gca, 'XDir', 'reverse');           % Show X axis from high to low
legend('Validation Error', 'Training Error', 'Best Lambda', 'Location', 'best', 'Interpreter', 'latex');
% Improve ticks and font
% set(gca, 'FontSize', 17, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
set(gca,'XScale','log','YScale','log','XDir','reverse',...
        'FontSize',17,'LineWidth',1.2,...
        'TickLabelInterpreter','latex');

% Test on full dataset with best lambda
% % excluding intercept
D = blkdiag(eye(n_features-1), 0);

bestTransMatrix = yTrain_Valid * xTrain_Valid' / (xTrain_Valid * xTrain_Valid' + bestlambda * D);

% % % testing daat set
Y_predtest = bestTransMatrix * TestDataSet.xTest;
Sigma_y=Covariance_Calculation(xTrain_Valid,yTrain_Valid,0.2,bestlambda);
% % % --- Covariance Calculation ---
% % % --- Dimensions: Features (p=71), Observations (n=168) ---
% % [p, n] = size(xTrain_Valid); 
% % % 1. Components for the Ridge Inverse
% % XTX = xTrain_Valid * xTrain_Valid';    % (72x72)
% % I_reg = bestlambda * D;                % (72x72)
% % XtX_inv = inv(XTX + I_reg);            % (72x72)
% % % Ridge Hat Matrix (H)
% % % Dimensions: (168x71) * (72x72) * (72x168) = (168x168)
% % H = xTrain_Valid' * XtX_inv * xTrain_Valid;
% % % Effective Degrees of Freedom
% % % This will be a value between 0 and 72
% % df = trace(H);
% % % Residual Variance (sigma2)
% % % yTrain_Valid is (71x168), bestTransMatrix is (71x71)
% % resid = yTrain_Valid - (bestTransMatrix * xTrain_Valid);
% % sse = sum(resid(:).^2);
% % sigma2 = sse / (n - df);               % Residual variance
% % % Covariance of Ridge Estimator (Cov_beta)
% % % Dimensions: (72x72) * (72x72) * (72x72) = (72x72)
% % % This represents the uncertainty in 71x71 Transformation Matrix
% % Cov_beta = sigma2 * (XtX_inv * XTX * XtX_inv);
% % % Variance of each of the 71 coefficients
% % VarBeta = diag(Cov_beta);

% Add to output
% TestDataSet.sigma2 = sigma2;
% TestDataSet.Cov_beta = Cov_beta;
% TestDataSet.df = df;
% TestDataSet.VarBeta = VarBeta;


% % % Un normalization 
% Y_predtestUnNorm = Y_predtest.*NormFactors.ys+NormFactors.yc; 
% 
% % % actual data set 
% TestDataSet.yTestUnNorm = TestDataSet.yTest.*NormFactors.ys+NormFactors.yc; 
% 
% % % % X un normalization
% X_predtestUnNorm = TestDataSet.xTest(1:end-1,1:end).*NormFactors.us+NormFactors.uc; 
% 
% % % % rmse error
% TestError= rmse(TestDataSet.yTest,Y_predtest,'all');
% 
% % % Final data set
PredyTrain_Valid = bestTransMatrix * xTrain_Valid;
% 
% % % Un normalization 
% PredyTrain_ValidUnNorm = PredyTrain_Valid.*NormFactors.ys+NormFactors.yc; 
% 
% % % actual data set Un normalization 
% yTrain_ValidUnNorm = yTrain_Valid.*NormFactors.ys+NormFactors.yc; 
% 
% % % % X un normalization
% X_predtrainUnNorm = xTrain_Valid(1:end-1,1:end).*NormFactors.us+NormFactors.uc; 

if isequal(AngularDomain,'on')
% % % Rounndess Measurements
% % % Initialize
Visibility= 'on';
ActualTestSamples= zeros(1,size(TestDataSet.yTest,2));
PredTestSamples= zeros(1,size(Y_predtest,2));
% % % Test roundess
for TestRnd= 1:size(Y_predtest,2)
    ActualTestSamples(TestRnd)= roundness(TestDataSet.yTest(:,TestRnd));

    PredTestSamples(TestRnd)= roundness(Y_predtest(:,TestRnd));
end
% % % Plot
% figure('Visible',Visibility); plot(PredTestSamples,ActualTestSamples,'o');
% title('Test Samples'); xlabel('Predicted Samples'); ylabel('Actual Samples');
% lsline;
% % Annotate RMSE
% rmse_val = rmse(PredTestSamples,ActualTestSamples);
% text(0.05, 0.95, ...
%     ['\textbf{RMSE: }', num2str(rmse_val, '%.3f'), ' $\mu$m'], ...
%     'Units', 'normalized', ...
%     'FontSize', 18, ...
%     'Interpreter', 'latex', ...
%     'Color', 'k', ...
%     'HorizontalAlignment', 'left');
[~, RelErr_HF, ~] = computeAllErrors(ActualTestSamples, PredTestSamples, 0.3,'off');
figure('Visible','on','Position',[10 10 600 500]);
% Plot HF vs. MatrRnd
scatter(PredTestSamples, ActualTestSamples, 50, 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
hold on;
xlabel('Data Driven Process Model (Roundness)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('High-Fidelity Process Model (Roundness)', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison of Roundness Predictions', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'latex');
rmse_val = sqrt(mean((ActualTestSamples - PredTestSamples).^2));
lims = [0, max([ActualTestSamples(:); PredTestSamples(:)]) * 1.05];
axis([lims(1), lims(2), lims(1), lims(2)]); 
% Add least squares fit
lsline; 
h = findobj(gca,'Type','Line');
set(h(1), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); % Format lsline
% legend({...
%     ['(Relative Error: ', num2str(mean(RelErr_HF)*100, '%.2f'), '\%, RMSE: ', num2str(rmse_val, '%.4f'), '\%)'], ...
%     'Least Squares Fit'}, ...
%     'FontSize', 14, 'Location', 'southeast', 'Interpreter', 'latex');
legend({...
    ['$\mathrm{Relative\ Error}: ', num2str(mean(RelErr_HF)*100, '%.2f'), '\%$, ', ...
    '$\mathrm{RMSE}: ', num2str(rmse_val, '%.4f'), '\,\mu m$'], ...
    'Least Squares Fit'}, ...
    'FontSize', 14, 'Location', 'southeast', 'Interpreter', 'latex');
% legend boxoff
axis square; % Equal aspect ratio
% Improve ticks and font
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');

% % % Train
ActualTrainSamples= zeros(1,size(yTrain_Valid,2));
PredTrainSamples= zeros(1,size(PredyTrain_Valid,2));
% % % % % % % % % % % % % % % % % % % % % % % % % %
for TrainRnd= 1:size(yTrain_Valid,2)
    ActualTrainSamples(TrainRnd)= roundness(yTrain_Valid(:,TrainRnd));

    PredTrainSamples(TrainRnd)= roundness(PredyTrain_Valid(:,TrainRnd));
end
% % % Plot
% % % figure('Visible',Visibility); plot(PredTrainSamples,ActualTrainSamples,'o');
% % % title('Train Samples'); xlabel('Predicted Samples'); ylabel('Actual Samples');
% % % lsline;
% % % % Annotate RMSE
% % % rmse_val = rmse(PredTrainSamples,ActualTrainSamples);
% % % text(0.05, 0.95, ...
% % %     ['\textbf{RMSE: }', num2str(rmse_val, '%.3f'), ' $\mu$m'], ...
% % %     'Units', 'normalized', ...
% % %     'FontSize', 18, ...
% % %     'Interpreter', 'latex', ...
% % %     'Color', 'k', ...
% % %     'HorizontalAlignment', 'left');
[~, RelErr_HF, ~] = computeAllErrors(ActualTrainSamples, PredTrainSamples, 0.3,'off');
figure('Visible','on','Position',[10 10 600 500]);
% Plot HF vs. MatrRnd
scatter(PredTrainSamples, ActualTrainSamples, 50, 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
hold on;
xlabel('Data Driven Process Model (Roundness)', 'FontSize', 14, 'Interpreter', 'latex');
ylabel('High-Fidelity Process Model (Roundness)', 'FontSize', 14, 'Interpreter', 'latex');
title('Comparison of Roundness Predictions', 'FontSize', 15, 'FontWeight', 'bold', 'Interpreter', 'latex');
rmse_val = sqrt(mean((ActualTrainSamples - PredTrainSamples).^2));
lims = [0, max([ActualTrainSamples(:); PredTrainSamples(:)]) * 1.05];
axis([lims(1), lims(2), lims(1), lims(2)]); 
% Add least squares fit
lsline; 
h = findobj(gca,'Type','Line');
set(h(1), 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); % Format lsline
% legend({...
%     ['(Relative Error: ', num2str(mean(RelErr_HF)*100, '%.2f'), '\%, RMSE: ', num2str(rmse_val, '%.4f'), '\%)'], ...
%     'Least Squares Fit'}, ...
%     'FontSize', 14, 'Location', 'southeast', 'Interpreter', 'latex');
legend({...
    ['$\mathrm{Relative\ Error}: ', num2str(mean(RelErr_HF)*100, '%.2f'), '\%$, ', ...
    '$\mathrm{RMSE}: ', num2str(rmse_val, '%.4f'), '\,\mu m$'], ...
    'Least Squares Fit'}, ...
    'FontSize', 14, 'Location', 'southeast', 'Interpreter', 'latex');
% legend boxoff
axis square; % Equal aspect ratio
% Improve ticks and font
set(gca, 'FontSize', 14, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');

% % % % angular plots Training
figure('Visible',Visibility,'Position',[10 10 1500 900]);
subplot(1,2,1);
polarplot(X_predtestUnNorm(:,1)+10*mean(abs(X_predtestUnNorm(:,1))),'linewidth', 1.5);
title('Initial Experimental Profile', 'FontSize', 18, Interpreter='latex');
subplot(1,2,2);
PredProfile= Y_predtestUnNorm(:,1)-mean(Y_predtestUnNorm(:,1));
polarplot(PredProfile,'linewidth', 1.5);
hold on;
ActProfile= -(TestLFData.MatrFinalProfile{1,1}-mean(TestLFData.MatrFinalProfile{1,1}));
polarplot(ActProfile,'linewidth', 1.5)
hold on; 
ActProfile= -(TestLFData.MatrFinalProfile{1,1}-mean(TestLFData.MatrFinalProfile{1,1}));
polarplot(ActProfile,'linewidth', 1.5)
title('Comparison of Workpiece Profiles', 'FontSize', 18, Interpreter='latex');
legend('Estimated Profile (Data-Driven)', 'Affine Transformation', ...
    'Interpreter', 'latex', 'Location', 'southwest');
Minimum= min([-2*max(PredProfile)+3*min(PredProfile) -2*max(ActProfile)+3*min(ActProfile)]);
Maximum= max([max(PredProfile)*1 max(ActProfile)*1]);
rlim([real(Minimum) real(Maximum)]);
% Improve ticks and font
set(gca, 'FontSize', 18, 'LineWidth', 1.2, 'TickDir', 'out', 'Box', 'off');
else

% % % harmonic range 
Range=size(TestDataSet.xTest,1)-1;

% % % % % Transformation to angular domain
% % % training samples complex to angular domain transformation
% % roundness measuremnts using separate functions
% % % Initialize
% % % LOOP
for Train = 1:size(xTrain_Valid,2)

    Rmse_train(Train) = rmse(yTrain_Valid(:,Train),PredyTrain_Valid(:,Train));

    % rel_err on abs data (new approach, a hybrid error)
    % % % MEHE Error
    MEHEerrors_train(:,Train) = abs(abs(PredyTrain_Valid(:,Train))- abs(yTrain_Valid(:,Train))) ./ (0.5 + abs(abs(yTrain_Valid(:,Train)))); % Element-wise computation of the error

    RelErrabsdata_train(Train)  = max( MEHEerrors_train(:,Train)); % Maximum normalized error

    % Roundness estimation
    [Rndtrain(Train),Estimations]= freqToSpatialRoundness(yTrain_Valid(:,Train),NsamplePoints);

    % % % 
    [InitRndtrain(Train),InitEstimations]= freqToSpatialRoundness(xTrain_Valid(1:end-1,Train),NsamplePoints);

    % % storing geo prof
    TrainGeoProf(:,Train)= Estimations.RealGeoProf;

    InitTrainGeoProf(:,Train)= InitEstimations.RealGeoProf;

    % % Roundness estimation prediction
    [RndtrainPredicted(Train),Estimations]= freqToSpatialRoundness(PredyTrain_Valid(:,Train),NsamplePoints);

    % % storing geo prof predicted
    TrainGeoProfPredicted(:,Train)= Estimations.RealGeoProf;

    % Rel Err Roundness (adjueted rel errs, a new approach)
    RelErr_trainRnd(Train) =abs(Rndtrain(Train) - RndtrainPredicted(Train)) ./ (abs(Rndtrain(Train) ) + 0.5);

    % Rel Err on roundness before harmonic reduction (i.e.m,
    % measured comparison with clipped profile

    % Indices management
    % % %
    FinaltrainRnd_fold_noclipping= TrainRndMeas;

    RelErr_trainRnd_comparisonwithnoclipping(Train) = abs(FinaltrainRnd_fold_noclipping(Train) - RndtrainPredicted(Train)) ./ (abs(FinaltrainRnd_fold_noclipping(Train)) + 0.5);
end

% % % % % % % % % %
% % rms error
Rmse_valid= zeros(1,size(TestDataSet.xTest,2));
% % a hybrid error
MEHEerrors_valid= zeros(Range,size(TestDataSet.xTest,2));
% %
RelErrabsdata_valid= zeros(1,size(TestDataSet.xTest,2));
Rndvalid= zeros(1,size(TestDataSet.xTest,2));
validGeoProf= zeros(NsamplePoints,size(TestDataSet.xTest,2));
RndvalidPredicted= zeros(1,size(TestDataSet.xTest,2));
validGeoProfPredicted= zeros(NsamplePoints,size(TestDataSet.xTest,2));
RelErr_validRnd= zeros(1,size(TestDataSet.xTest,2));
RelErr_validRnd_comparisonwithnoclipping= zeros(1,size(TestDataSet.xTest,2));
% % Testing samples
for Test = 1:size(TestDataSet.xTest,2)

    Rmse_valid(Test) = rmse(TestDataSet.yTest(:,Test),Y_predtest(:,Test));

    % % % MEHE Error
    MEHEerrors_valid(:,Test) = abs(abs(Y_predtest(:,Test))- abs(TestDataSet.yTest(:,Test))) ./ (0.5 + abs(abs(TestDataSet.yTest(:,Test)))); % Element-wise computation of the error

    RelErrabsdata_valid(Test)  = max(MEHEerrors_valid(:,Test)); % Maximum normalized error

    %  Roundness
    [Rndvalid(Test),Estimations]= freqToSpatialRoundness(TestDataSet.yTest(:,Test),NsamplePoints);

    validGeoProf(:,Test)= Estimations.RealGeoProf;

    % % % % 
    [InitRndvalid(Test),InitEstimations]= freqToSpatialRoundness(TestDataSet.xTest(1:end-1,Test),NsamplePoints);

    InitvalidGeoProf(:,Test)= InitEstimations.RealGeoProf;

    [RndvalidPredicted(Test),Estimations]= freqToSpatialRoundness(Y_predtest(:,Test),NsamplePoints);

    validGeoProfPredicted(:,Test)= Estimations.RealGeoProf;

    % Rel Err Roundness (adjusted neew approach)
    % % % % RelErr_validRnd(Test) =norm(Rndvalid(Test)-RndvalidPredicted(Test)) / norm(Rndvalid(Test));
    RelErr_validRnd(Test) =abs(Rndvalid(Test) - RndvalidPredicted(Test)) ./ (abs(Rndvalid(Test)) + 1);

    % Rel Err on roundness before harmonic reduction (i.e.m,
    % measured
    FinalRnd= TestRndMeas;

    % Indices management per fold and the validation
    FinalvalidRnd_fold_noclipping= FinalRnd;

    % % % RelErr_validRnd_comparisonwithnoclipping(Test) =norm(FinalvalidRnd_fold_noclipping(Test)-RndvalidPredicted(Test)) / norm(FinalvalidRnd_fold_noclipping(Test));
    RelErr_validRnd_comparisonwithnoclipping(Test) = abs(FinalvalidRnd_fold_noclipping(Test) - RndvalidPredicted(Test)) ./ (abs(FinalvalidRnd_fold_noclipping(Test)) + 0.5);
end

% % % % % % % % % % % % % % %
% % train samples roundness
TrainingSamples= TrainRndMeas;
% % predicted train samples roundness
Predictedtrain_samples= RndtrainPredicted;

% % validation samples roundness
ValidSamples= TestRndMeas;
% % predicted validation samples roundness
Predictedvalid_samples= RndvalidPredicted;

Actual_TrainGeoProf = (TrainGeoProf);

% % Data Out
DataOut.Actual_TrainGeoProf= Actual_TrainGeoProf;

Pred_TrainGeoProf = (TrainGeoProfPredicted);

% % Data out prediction
DataOut.Pred_TrainGeoProf= Pred_TrainGeoProf;

% Valid (avg across all folds)
Actual_validGeoProf = (validGeoProf);

% % Data Out (validation)
DataOut.Actual_validGeoProf = Actual_validGeoProf;

Pred_validGeoProf = (validGeoProfPredicted);

% % Data Out (validation) prediction
DataOut.Pred_validGeoProf = Pred_validGeoProf;

% for Plot abs of average complex data across all folds
Actual_validabsProf = (TestDataSet.yTest);

Pred_validabsProf = (Y_predtest);

% % % train
Actual_trainabsProf = (yTrain_Valid);

Pred_trainabsProf = (PredyTrain_Valid);


% % Training dATA
OutPutData_DDtrain_Roundness=zeros(1,size(Actual_TrainGeoProf,2));
OutPutData_DDtrainPred_Roundness=zeros(1,size(Actual_TrainGeoProf,2));
IterRnd_train=zeros(1,size(Actual_TrainGeoProf,2));
MatrRnd_train=zeros(1,size(Actual_TrainGeoProf,2));

for TrainProf=1:size(Actual_TrainGeoProf,2)

    if isempty(LFData.LFMRR)
    else
        ax1=figure('Visible',Visibility,'Position',[10 10 1100 700]);

        % % Post Processing
        % % % % Theor MRR Estimation
        beta_s= asin(2*TrainSimData(TrainProf).MachProc.hw/(TrainSimData(TrainProf).dw+TrainSimData(TrainProf).MachProc.OperW.Diam));

        betaCW=asin(2*TrainSimData(TrainProf).MachProc.hw/(TrainSimData(TrainProf).MachProc.ContrW.Diam+TrainSimData(TrainProf).dw));

        Beta= beta_s+betaCW;

        Alpha= pi/2 - TrainSimData(TrainProf).MachProc.Gamma*pi/180 - beta_s;

        Radial_Infeed=  (InpParV(TrainProf, 3) / 1000); % mm

        Infeed= Radial_Infeed/(1-sin(Beta)/sin((Beta+Alpha))+sin(Alpha)/sin((Beta+Alpha))); % mm

        % % % % % MRR Estimation
        TheorMRR_REAL= Infeed*TrainSimData(TrainProf).Workpiece_rps*TrainSimData(TrainProf).MachProc.wplength*pi*TrainSimData(TrainProf).dw; % mm3/sec

        % % % MRR
        % % figure('Visible',Visibility,'Position',[10 10 1500 900]); subplot(1,2,1);
        ID= find(trainIdx1); ID= ID(TrainProf);

        Operation= ((ID));

        plot(TrainSimData(TrainProf).SimMRR' * TrainSimData(TrainProf).MachProc.wplength, 'DisplayName', sprintf('HF Simulation MRR (Avg: %.2f)', mean(TrainSimData(TrainProf).SimMRR((TrainSimData(TrainProf).Nrev/2+2)*SimData(TrainProf).Nwp:((TrainSimData(TrainProf).Nrev/2)+3)*TrainSimData(TrainProf).Nwp) * TrainSimData(TrainProf).MachProc.wplength)),LineWidth=1.5);

        hold on;
        if isequal(LFData.LFMRR,[])
        else
            plot(TrainLFData.LFMRR{1,TrainProf}, 'DisplayName', sprintf('Iterative MRR (Avg: %.2f)', mean(TrainLFData.LFMRR{1,TrainProf}((TrainSimData(TrainProf).Nrev/2+2)*TrainSimData(TrainProf).Nwp:((TrainSimData(TrainProf).Nrev/2)+3)*TrainSimData(TrainProf).Nwp))),LineWidth=1.5);
        end
        % plot(LFData.LFMRR{1, ID}, 'DisplayName', sprintf('HF Simulation MRR (Avg: %.2f)', mean(LFData.LFMRR{1, ID}((SimData(ID).Nrev/2+2)*SimData(ID).Nwp:((SimData(ID).Nrev/2)+3)*SimData(ID).Nwp))),LineWidth=1.5);
        yline(TheorMRR_REAL, 'r--', 'DisplayName', sprintf('TheorMRR REAL: %.2f', TheorMRR_REAL));
        xlabel('Steps', 'FontSize', 12);
        ylabel('mm³/sec', 'FontSize', 12);
        % title('MRR', 'FontSize', 12);
        title(['MRR Work-height (mm) = ', num2str(InpParV(ID,1)), ...
            ' | R-inf (\mu m) = ', num2str(InpParV(ID,3)), ...
            ' | Wp rps = ', num2str(SimData(ID).Workpiece_rps), ...
            ' | Wp Diameter = ', num2str(InpParV(ID,2)),' | Operation No = ', num2str(Operation),' | Gamma = ', num2str(InpParV(ID,5))], ...
            'FontSize', 12, 'Interpreter', 'tex');
        grid on;
        legend('show', 'FontSize', 12);
    end

    ax2= figure('Visible',Visibility,'Position',[10 10 1500 900]);

    if isequal(LFData.IterFinalProfile,[]) || isequal(LFData.MatrFinalProfile,[])
        subplot(1,2,1);
    else
    end
    % % % Initial Profile addind (23/08/2025)
    FullActualInitvalidGeoProfTrain= (FullTestDataSet.xTrain(1:end,TrainProf));
    polarplot(FullActualInitvalidGeoProfTrain, 'LineWidth', 2, 'Color', 'k');
    hold on; 
    ActualInitGeoProf= (InitTrainGeoProf(1:end,TrainProf));
    polarplot(ActualInitGeoProf, 'LineWidth', 1, 'Color', 'g');
    title('Initial Experimental Profile', 'FontSize', 18,Interpreter='latex');
    Minimum= min([-2*max(ActualInitGeoProf)+3*min(ActualInitGeoProf)]);
    Maximum= max([max(ActualInitGeoProf)*1]);
    rlim([real(Minimum) real(Maximum)]);
    set(gca, 'FontSize', 17);
    subplot(1,2,2);
    % HFProf=-SimData(ValidProf).SimFinalProfile+1*mean(abs(-SimData(ValidProf).SimFinalProfile));
    % ActualTrainGeoProf= (Actual_TrainGeoProf(1:end,TrainProf))+1*mean(abs((Actual_TrainGeoProf(1:end,TrainProf))));
    % ActualTrainGeoProf= (Actual_TrainGeoProf(1:end,TrainProf));
    ActualTrainGeoProf= (abs(FullTestDataSet.yTrain(1:end,TrainProf)))-mean(abs(FullTestDataSet.yTrain(1:end,TrainProf))); 
    % % storing the train roundness
    OutPutData_DDtrain_Roundness(TrainProf)= roundness(ActualTrainGeoProf);
    h1=polarplot(ActualTrainGeoProf, 'LineWidth', 2, 'Color', 'r');
    hold on;
    % PredTrainGeoProf= (Pred_TrainGeoProf(1:end,TrainProf))+1*mean(abs((Pred_TrainGeoProf(1:end,TrainProf))));
    PredTrainGeoProf= -(Pred_TrainGeoProf(1:end,TrainProf));
    % % storing predicted roundness
    OutPutData_DDtrainPred_Roundness(TrainProf)= roundness(PredTrainGeoProf);
    h2=polarplot(PredTrainGeoProf, '--', 'LineWidth', 2, 'Color', 'r');
    hold on; 
    % LFProf= LFData.IterFinalProfile(:,ID);
    if isequal(LFData.IterFinalProfile,[])
        LFProf=[];
    else
        CumDataMat.SpaceDomainProf= TrainLFData.IterFinalProfile{TrainProf};
        % CumDataMat.SpaceDomainProf= InitialProfiles(:,1:end);
        W= ones(size(CumDataMat.SpaceDomainProf,1),1);
        % % [Errors,DataOut] = NoFroundness(DataMat,num_harmonics_range,Plot,use_fixed_harmonics,Save,Case,W)
        [~,FreqDataOut] = NoFroundness(CumDataMat,Range(end),0,1,0,'LF Iter Profile',W);
        % LFProf= LFData.IterFinalProfile{ID}+1*mean(abs(LFData.IterFinalProfile{ID}));
        LFProf= FreqDataOut.reduced_harmonicGeoProf;
        % % storing iter rndness
        IterRnd_train(TrainProf)= roundness(LFProf);
        % % %
        h3=polarplot((LFProf),"linewidth",1);
        hold on;
    end
    % % Matr profile
    if isequal(LFData.MatrFinalProfile,[])
        MatrProf= [];
    else
        subplot(1,2,2);
        CumDataMat.SpaceDomainProf= TrainLFData.MatrFinalProfile{TrainProf};

        MatrProf= (CumDataMat.SpaceDomainProf-mean(CumDataMat.SpaceDomainProf)); 
        % CumDataMat.SpaceDomainProf= InitialProfiles(:,1:end);
        % W= ones(size(CumDataMat.SpaceDomainProf,1),1);
        % % % [Errors,DataOut]  = NoFroundness(DataMat,num_harmonics_range,Plot,use_fixed_harmonics,Save,Case,W)
        % [~,FreqDataOut] = NoFroundness(CumDataMat,Range(end),0,1,0,'LF Matr Profile',W);
        % MatrProf=  -(FreqDataOut.reduced_harmonicGeoProf);
        % % % storing iter rndness
        MatrRnd_train(TrainProf)= roundness(MatrProf);
        % % % MatrProf= LFData.MatrFinalProfile{ID}+1*mean(abs(LFData.MatrFinalProfile{ID}));
        % % % PredValidGeoProf= (Pred_validGeoProf(1:end,ValidProf))+1*mean(abs((Pred_validGeoProf(1:end,ValidProf))));
        h4=polarplot((MatrProf), 'LineWidth', 2,'Color', 'b' );
        grid on;
        % title('Transition Matrix Process Model: Training', 'FontSize', 18);
        % Minimum= min([-2*max(LFProf)+3*min(LFProf) -2*max(MatrProf)+3*min(MatrProf)]);
        % Maximum= max([max(LFProf)*1 max(MatrProf)*1]);
        % rlim([real(Minimum) real(Maximum)]);
        title('Profiles Comparison', 'FontSize', 18,Interpreter='latex');
        Minimum= min([-2*max(ActualTrainGeoProf)+3*min(ActualTrainGeoProf) -2*max(PredTrainGeoProf)+3*min( PredTrainGeoProf) -2*max(MatrProf)+3*min(MatrProf)]);
        Maximum= max([max(ActualTrainGeoProf)*1 max(PredTrainGeoProf)*1 max(MatrProf)*1]);
        rlim([real(Minimum) real(Maximum)]);
    end

    if isequal(LFData.IterFinalProfile,[]) &&  isequal(LFData.MatrFinalProfile,[])
        legend(['HF Rnd=',num2str(OutPutData_DDtrain_Roundness(TrainProf)),'\mum'],['Data-Driven- Rnd=',num2str(OutPutData_DDtrainPred_Roundness(TrainProf)),'\mum'] ...
            ,'FontSize', 14,Location='southoutside');
    elseif isequal(LFData.IterFinalProfile,[])
        % legend([h1,h2],{['HF Rnd=',num2str(OutPutData_DDtrain_Roundness(TrainProf)),'\mum'],['Data-Driven- Rnd=',num2str(OutPutData_DDtrainPred_Roundness(TrainProf)),'\mum']} ...
        %     ,'FontSize', 14,Location='southoutside');
        % legend([h4],{['LF Matr- Rnd=',num2str(MatrRnd_train(TrainProf)),'\mum'] },'FontSize', 14,Location='southoutside');
        legend('Ground truth (High Fidelity)','Prediction: Data-driven Affine Transformation', ...
            'Physics Based Affine Transformation', ...
            'Position',[0.662003786314762 0.519145581248412 0.165200003128052 0.0460858595491659], 'Interpreter', 'latex');
    elseif isequal(LFData.MatrFinalProfile,[])
        legend(['HF Rnd=',num2str(OutPutData_DDtrain_Roundness(TrainProf)),'\mum'],['Data-Driven- Rnd=',num2str(OutPutData_DDtrainPred_Roundness(TrainProf)),'\mum'] ...
            ,['LF Iter- Rnd=',num2str(IterRnd_train(TrainProf)),'\mum'] ,'FontSize', 14,Location='southoutside');
    end
    set(gca, 'FontSize', 17);
    if Save
        saveas(ax2, ['DD_TrainingGeoProf', num2str(TrainProf), '.png']);
        savefig(ax2,['DD_TrainingGeoProf', num2str(TrainProf), '.fig']);
        % Save the plot for training data
        if isequal(LFData.IterFinalProfile,[])  || isequal(LFData.MatrFinalProfile,[])
            saveas(ax1, ['DD_TrainingMRR', num2str(TrainProf), '.png']);
            savefig(ax1,['DD_TrainingMRR', num2str(TrainProf), '.fig']);
            % else
            %     saveas(gcf, ['Training Sample_MRR_GeoProf_LF_HF_DD', num2str(TrainProf), '_plot.png']);
            %     savefig(gcf,['Training Sample_MRR_GeoProf_LF_HF_DD', num2str(TrainProf), '_plot.fig']);
        end
    end
    hold off;
     % Pick random ID once
    if TrainProf == 1
        rng('shuffle'); 
        sampleID = randi([1,size(Actual_TrainGeoProf,2)]);
        disp(['Training Random Sample ID: ', num2str(sampleID)]);
    end
end

% Validation Samples
% % Initalize
OutPutData_DDvalidPred_Roundness= zeros(1,size(Actual_validabsProf,2));
OutPutData_DDvalid_Roundness= zeros(1,size(Actual_validabsProf,2));
IterRnd_valid= zeros(1,size(Actual_validabsProf,2));
MatrRnd_valid= zeros(1,size(Actual_validabsProf,2));

for ValidProf=1:size(Actual_validabsProf,2)

    if isempty(LFData.LFMRR)
    else
        ax1= figure('Visible',Visibility,'Position',[10 10 1100 700]);

        % % Post Processing
        % % % % Theor MRR Estimation
        beta_s= asin(2*TestSimData(ValidProf).MachProc.hw/(TestSimData(ValidProf).dw+TestSimData(ValidProf).MachProc.OperW.Diam));

        betaCW=asin(2*TestSimData(ValidProf).MachProc.hw/(TestSimData(ValidProf).MachProc.ContrW.Diam+TestSimData(ValidProf).dw));

        Beta= beta_s+betaCW;

        Alpha= pi/2 - TestSimData(ValidProf).MachProc.Gamma*pi/180 - beta_s;

        Radial_Infeed=  (InpParV(ValidProf, 3) / 1000); % mm

        Infeed= Radial_Infeed/(1-sin(Beta)/sin((Beta+Alpha))+sin(Alpha)/sin((Beta+Alpha))); % mm

        % % % % % MRR Estimation
        TheorMRR_REAL= Infeed*TestSimData(ValidProf).Workpiece_rps*TestSimData(ValidProf).MachProc.wplength*pi*TestSimData(ValidProf).dw; % mm3/sec

        % % % MRR
        % % figure('Visible',Visibility,'Position',[10 10 1500 900]); subplot(1,2,1);
        ID= find(testIdx); ID= ID(ValidProf);

        Operation= ((ID));

        plot(TestSimData(ValidProf).SimMRR'  * TestSimData(ValidProf).MachProc.wplength, 'DisplayName', sprintf('HF Simulation MRR (Avg: %.2f)', mean(TestSimData(ValidProf).SimMRR((TestSimData(ValidProf).Nrev/2+2)*TestSimData(ValidProf).Nwp:((TestSimData(ValidProf).Nrev/2)+3)*TestSimData(ValidProf).Nwp) * TestSimData(ValidProf).MachProc.wplength)),LineWidth=1.5);

        hold on;
        if isequal(LFData.LFMRR,[])
        else
            plot(TestLFData.LFMRR{1,ValidProf}, 'DisplayName', sprintf('Iterative MRR (Avg: %.2f)', mean(TestLFData.LFMRR{1,ValidProf}((TestSimData(ValidProf).Nrev/2+2)*TestSimData(ValidProf).Nwp:((TestSimData(ValidProf).Nrev/2)+3)*TestSimData(ValidProf).Nwp))),LineWidth=1.5);
        end
        % plot(LFData.LFMRR{1, ID}, 'DisplayName', sprintf('HF Simulation MRR (Avg: %.2f)', mean(LFData.LFMRR{1, ID}((SimData(ID).Nrev/2+2)*SimData(ID).Nwp:((SimData(ID).Nrev/2)+3)*SimData(ID).Nwp))),LineWidth=1.5);
        yline(TheorMRR_REAL, 'r--', 'DisplayName', sprintf('TheorMRR REAL: %.2f', TheorMRR_REAL));
        xlabel('Steps', 'FontSize', 12);
        ylabel('mm³/sec', 'FontSize', 12);
        % title('MRR', 'FontSize', 12);
        title(['MRR Work-height (mm) = ', num2str(InpParV(ValidProf,1)), ...
            ' | R-inf (\mu m) = ', num2str(InpParV(ValidProf,3)), ...
            ' | Wp rps = ', num2str(TestSimData(ValidProf).Workpiece_rps), ...
            ' | Wp Diameter = ', num2str(InpParV(ValidProf,2)),' | Operation No = ', num2str(Operation),' | Gamma = ', num2str(InpParV(ValidProf,5))], ...
            'FontSize', 12, 'Interpreter', 'tex');
        grid on;
        legend('show', 'FontSize', 12);

    end

    % subplot(1,2,2);

    ax2= figure('Visible',Visibility,'Position',[10 10 1500 900]);
    if isequal(LFData.IterFinalProfile,[])
        subplot(1,2,1);
    else
    end
    % % % Initial Profile addind (23/08/2025)
    FullActualInitvalidGeoProf= (FullTestDataSet.xTest(1:end,ValidProf));
    polarplot(FullActualInitvalidGeoProf, 'LineWidth', 2, 'Color', 'k');
    hold on; 
    ActualInitvalidGeoProf= (InitvalidGeoProf(1:end,ValidProf));
    polarplot(ActualInitvalidGeoProf, 'LineWidth', 1, 'Color', 'g');
    title('Initial Experimental Profile', 'FontSize', 18,Interpreter='latex');
    Minimum= min([-2*max(ActualInitvalidGeoProf)+3*min(ActualInitvalidGeoProf) -2*max(FullActualInitvalidGeoProf)+3*min(FullActualInitvalidGeoProf)]);
    Maximum= max([max(ActualInitvalidGeoProf)*1 max(FullActualInitvalidGeoProf)*1]);
    rlim([real(Minimum) real(Maximum)]);
    set(gca, 'FontSize', 17);
    subplot(1,2,2);
    % HFProf=-SimData(ValidProf).SimFinalProfile+1*mean(abs(-SimData(ValidProf).SimFinalProfile));
    % ActualValidGeoProf= (Actual_validGeoProf(1:end,ValidProf))+1*mean(abs((Actual_validGeoProf(1:end,ValidProf))));
    % ActualValidGeoProf= (Actual_validGeoProf(1:end,ValidProf));
    ActualValidGeoProf= (abs(FullTestDataSet.yTest(1:end,ValidProf)))-mean(abs(FullTestDataSet.yTest(1:end,ValidProf))); 
    % % storing the valid roundness
    OutPutData_DDvalid_Roundness(ValidProf)= roundness(ActualValidGeoProf);
    h1=polarplot(ActualValidGeoProf, 'LineWidth', 2, 'Color', 'r');
    hold on;
    % PredValidGeoProf= (Pred_validGeoProf(1:end,ValidProf))+1*mean(abs((Pred_validGeoProf(1:end,ValidProf))));
    PredValidGeoProf= -(Pred_validGeoProf(1:end,ValidProf));
    % % storing predicted roundness
    OutPutData_DDvalidPred_Roundness(ValidProf)= roundness(PredValidGeoProf);
    h2=polarplot(PredValidGeoProf, '--', 'LineWidth', 2, 'Color', 'g');
    hold on; 
    title('Profiles Comparison', 'FontSize', 18,Interpreter='latex');
    % Minimum= min([-2*max(ActualValidGeoProf)+3*min(ActualValidGeoProf) -2*max(PredValidGeoProf)+3*min(PredValidGeoProf)]);
    % Maximum= max([max(ActualValidGeoProf)*1 max(PredValidGeoProf)*1]);
    % rlim([real(Minimum) real(Maximum)]);
    % LFProf= LFData.IterFinalProfile(:,ID);
    if isequal(LFData.IterFinalProfile,[])
        LFProf=[];
    else
        CumDataMat.SpaceDomainProf= TestLFData.IterFinalProfile{ValidProf};
        % CumDataMat.SpaceDomainProf= InitialProfiles(:,1:end);
        W= ones(size(CumDataMat.SpaceDomainProf,1),1);
        % % [Errors,DataOut] = [Errors,DataOut] = NoFroundness(DataMat,num_harmonics_range,Plot,use_fixed_harmonics,Save,Case,W)
        [~,FreqDataOut] = NoFroundness(CumDataMat,Range(end),0,1,0,'LF Iter Profile',W);
        % LFProf= LFData.IterFinalProfile{ID}+1*mean(abs(LFData.IterFinalProfile{ID}));
        LFProf= FreqDataOut.reducedCentered_harmonicGeoProf;
        % % % % % [LFProf, ~] = PhaseShift(ActualValidGeoProf, LFProf);
        % % storing iter rndness
        IterRnd_valid(ValidProf)= roundness(LFProf);
        % % %
        h3=polarplot((LFProf),"linewidth",1.5);
        hold on;
    end
    % % Matr profile
    if isequal(LFData.MatrFinalProfile,[])
        MatrProf=[];
    else
        subplot(1,2,2);
        CumDataMat.SpaceDomainProf= TestLFData.MatrFinalProfile{ValidProf};

        MatrProf= (CumDataMat.SpaceDomainProf-mean(CumDataMat.SpaceDomainProf));
        % % CumDataMat.SpaceDomainProf= InitialProfiles(:,1:end);
        % W= ones(size(CumDataMat.SpaceDomainProf,1),1);
        % % % [Errors,DataOut] = NoFroundness(DataMat,num_harmonics_range,Plot,use_fixed_harmonics,Save,Case,W)
        % [~,FreqDataOut] = NoFroundness(CumDataMat,Range(end),0,1,0,'LF Matr Profile',W);
        % MatrProf=  -(FreqDataOut.reduced_harmonicGeoProf);
        % % % % Shifting the profiles
        % [MatrProf, ~] = PhaseShift(ActualValidGeoProf, MatrProf);
        % % storing iter rndness
        MatrRnd_valid(ValidProf)= roundness(MatrProf);
        % % % MatrProf= LFData.MatrFinalProfile{ID}+1*mean(abs(LFData.MatrFinalProfile{ID}));
        % % % PredValidGeoProf= (Pred_validGeoProf(1:end,ValidProf))+1*mean(abs((Pred_validGeoProf(1:end,ValidProf))));
        h4=polarplot((MatrProf), 'LineWidth', 2,'Color', 'b' );
        grid on;
        title('Profiles Comparison', 'FontSize', 18,Interpreter='latex');
        % Minimum= min([-2*max(LFProf)+3*min(LFProf) -2*max(MatrProf)+3*min(MatrProf)]);
        % Maximum= max([max(LFProf)*1 max(MatrProf)*1]);
        % rlim([real(Minimum) real(Maximum)]);
        Minimum= min([-2*max(ActualValidGeoProf)+3*min(ActualValidGeoProf) -2*max(PredValidGeoProf)+3*min(PredValidGeoProf) -2*max(MatrProf)+3*min(MatrProf)]);
        Maximum= max([max(ActualValidGeoProf)*1 max(PredValidGeoProf)*1 -2*max(MatrProf)+3*min(MatrProf)]);
        rlim([real(Minimum) real(Maximum)]);
    end
    if isequal(LFData.IterFinalProfile,[]) && isequal(LFData.MatrFinalProfile,[])
        legend(['HF Rnd=',num2str(OutPutData_DDvalid_Roundness(ValidProf)),'\mum'],['Data-Driven- Rnd=',num2str(OutPutData_DDvalidPred_Roundness(ValidProf)),'\mum'] ...
            ,'FontSize', 14,Location='southoutside');
    elseif isequal(LFData.IterFinalProfile,[])
        % legend([h1,h2],{['HF Rnd=',num2str(OutPutData_DDvalid_Roundness(ValidProf)),'\mum'],['Data-Driven- Rnd=',num2str(OutPutData_DDvalidPred_Roundness(ValidProf)),'\mum']} ...
        %     ,'FontSize', 14,Location='southoutside');
        % legend([h4],{['LF Matr- Rnd=',num2str(MatrRnd_valid(ValidProf)),'\mum'] },'FontSize', 14,Location='southoutside');
        legend('Ground truth (High Fidelity)','Prediction: Data-driven Affine Transformation', ...
            'Physics Based Affine Transformation', ...
            'Position',[0.662003786314762 0.519145581248412 0.165200003128052 0.0460858595491659], 'Interpreter', 'latex');
    elseif isequal(LFData.MatrFinalProfile,[])
        legend(['HF Rnd=',num2str(OutPutData_DDvalid_Roundness(ValidProf)),'\mum'],['Data-Driven- Rnd=',num2str(OutPutData_DDvalidPred_Roundness(ValidProf)),'\mum'] ...
            ,['LF Iter- Rnd=',num2str(IterRnd_valid(ValidProf)),'\mum'] ,'FontSize', 14,Location='southoutside');
    end
    set(gca, 'FontSize', 17);
    if Save
        saveas(ax2, ['DD_ValidationSampleGeoProf', num2str(ValidProf), '.png']);
        savefig(ax2,['DD_ValidationSampleGeoProf', num2str(ValidProf), '.fig']);
        % Save the plot for training data
        if isequal(LFData.IterFinalProfile,[])  || isequal(LFData.MatrFinalProfile,[])
            saveas(ax1, ['DD_ValidationSample', num2str(ValidProf), '.png']);
            savefig(ax1,['DD_ValidationSample', num2str(ValidProf), '.fig']);
            % % % else
            % % %     saveas(gcf, ['Validation Sample_MRR_GeoProf_LF_HF_DD', num2str(ValidProf), '_plot.png']);
            % % %     savefig(gcf,['Validation Sample_MRR_GeoProf_LF_HF_DD', num2str(ValidProf), '_plot.fig']);
        end
    end
    hold off;
    % Pick random ID once
    if ValidProf == 1
        rng('shuffle'); 
        sampleID = randi([1,size(Actual_validabsProf,2)]);
        disp(['Validation Random Sample ID: ', num2str(sampleID)]);
    end
end


% % % figure
[~, RelErr_train, ~] = computeAllErrors(TrainingSamples, Predictedtrain_samples, 0.5,'off');
[~, RelErr_valid, ~] = computeAllErrors(ValidSamples, Predictedvalid_samples, 0.5,'off');
% % % roundness error estimation with LF iter
if isequal(LFData.IterFinalProfile,[]) || isequal(LFData.MatrFinalProfile,[])
    [~, RelErr_train_LFIter, ~] = computeAllErrors(TrainingSamples,MatrRnd_train, 0.5,'off');
    [~, RelErr_valid_LFIter, ~] = computeAllErrors(ValidSamples,MatrRnd_valid, 0.5,'off');
end

if Visibility
    Rnds= figure('Visible',Visibility,'Position',[10 10 1000 500]);
    ax1=subplot(1,2,1);
    plot(Predictedtrain_samples,TrainingSamples, 'o',LineWidth=1);
    xlabel('Data-Driven Process Model',FontSize=12);
    ylabel('High-Fidelity Process Model',FontSize=12);
    title(['Training Samples: ', num2str(size(xTrain_Valid,2))],FontSize=12);
    legend({['Roundness Error (Err: ', num2str(mean(RelErr_train), '%.5f'), ')']}, ...
        'FontSize', 12, 'Location', 'best');
    lims = [min([TrainingSamples'; Predictedtrain_samples']) max([TrainingSamples'; Predictedtrain_samples'])]; % Calculate common limits
    axis([lims(1) lims(2) lims(1) lims(2)]);
    grid on;
    hold on;
    ax2=subplot(1,2,2);
    plot(MatrRnd_train,TrainingSamples,'o',LineWidth=1);
    xlabel('Transition Matrix Process Model',FontSize=12);
    ylabel('High Fidelity Process Model',FontSize=12);
    title(['Training Samples: ', num2str(size(xTrain_Valid,2))],FontSize=12);
    legend({['Roundness Error (Err: ', num2str(mean(RelErr_train_LFIter), '%.5f'), ')']}, ...
        'FontSize', 12, 'Location', 'best');
    % lsline
    % lims = [min([MatrRnd_train'; Predictedtrain_samples']) max([MatrRnd_train'; Predictedtrain_samples'])]; % Calculate common limits
    % axis([lims(1) lims(2) lims(1) lims(2)]);
    grid on;
    % lims = [min([ValidSamples'; Predictedvalid_samples']) max([ValidSamples'; Predictedvalid_samples'])]; % Calculate common limits
    % axis([lims(1) lims(2) lims(1) lims(2)]);
    linkaxes([ax1,ax2]);
    % % % % %
    if isequal(LFData.IterFinalProfile,[])  || isequal(LFData.MatrFinalProfile,[])
        LFRnds= figure('Visible',Visibility,'Position',[10 10 1000 500]);
        ax1=subplot(1,2,1);
        plot(Predictedvalid_samples,ValidSamples, '*',LineWidth=1);
        xlabel('Data Driven Process Model',FontSize=12);
        ylabel('High-Fidelity Process Model',FontSize=12);
        title(['Validation Samples: ', num2str(size(TestDataSet.xTest,2))],FontSize=12);
        grid on;
        % legend('HF vs. Matr Rnd','HF vs. Iter Rnd',FontSize=12
        legend({['Roundness Error (Err: ', num2str(mean(RelErr_valid), '%.5f'), ')']}, ...
            'FontSize', 12, 'Location', 'best');
        ax2=subplot(1,2,2);
        plot(MatrRnd_valid,ValidSamples,'*',LineWidth=1);
        xlabel('Transition Matrix Process Model',FontSize=12);
        ylabel('High Fidelity Process Model',FontSize=12);
        title(['Validation Samples: ', num2str(size(TestDataSet.xTest,2))],FontSize=12);
        grid on;
        % legend('HF vs. Matr Rnd','HF vs. Iter Rnd',FontSize=12
        legend({['Roundness Error (Err: ', num2str(mean(RelErr_valid_LFIter), '%.5f'), ')']}, ...
            'FontSize', 12, 'Location', 'best');
        % lsline
        % lims = [min([MatrRnd_valid'; Predictedvalid_samples']) max([MatrRnd_valid'; Predictedvalid_samples'])]; % Calculate common limits
        % axis([lims(1) lims(2) lims(1) lims(2)]);
        linkaxes([ax1 ax2]);
    end
end


% % % % % % % % % % % % % % % % % % % % % % % % %
% % plots
figure; bar(svd(bestTransMatrix));
xlabel('Singular Values',FontSize=12);
title('SVD Data Driven',FontSize=12);

% % % % training
x= round(linspace(0,length(Actual_trainabsProf(:,1)),length(Actual_trainabsProf(:,1))));

MEHEerrors_hfIter=zeros(numel(x)-2,size(Actual_trainabsProf,2)); MaxMEHE_HFIter=zeros(1,size(Actual_trainabsProf,2));

MEHEerrors_hfMatr=zeros(numel(x)-2,size(Actual_trainabsProf,2));  MaxMEHE_HFMatr=zeros(1,size(Actual_trainabsProf,2));

for trainProf=1:size(Actual_trainabsProf,2)
    figure('Visible',Visibility,'Position',[10 10 700 500]);
    if isempty(LFData.LFMRR)
    else
        % % % % %
        % [~, MeasFFT] = fftg(1/SimData(ValidProf).Nwp, SimData(ValidProf).InitialProfiles(1:SimData(ValidProf).Nwp,1), 0);
        plot(x(3:end),abs(xTrain_Valid(3:x(end),trainProf)),'LineWidth',1.5,'Color','k');
    end
    hold on; 
    plot(x(3:end),abs(Actual_trainabsProf(3:end,trainProf)), 'LineWidth', 2, 'Color', 'r');
    hold on;
    plot(x(3:end),abs(Pred_trainabsProf(3:end,trainProf)), '--', 'LineWidth', 2, 'Color', 'g');
    hold on;
    % ID= DataOut.ValidIndices(ValidProf,MinValidId);
    % if isequal(LFData.IterFFT,[])
    % else
    %     plot(x(3:end),abs(LFData.IterFFT(3:x(end), ID)),LineWidth=1);
    %     MEHEerrors_hfIter(:,ValidProf) = abs(abs(LFData.IterFFT(3:x(end), ID))- abs(Actual_validabsProf(3:end,ValidProf))) ./ (0.5 + abs(abs(Actual_validabsProf(3:end,ValidProf)))); % Element-wise computation of the error
    %     MaxMEHE_HFIter(ValidProf) = max(MEHEerrors_hfIter(:,ValidProf)); % Maximum normalized error
    %     hold on;
    % end 
    if isequal(LFData.MatrFinalProfile,[])
    else
        plot(x(3:end),abs(TrainLFData.MatrFFT(3:x(end), trainProf)),'LineWidth',1.5,'Color', 'b' );
        MEHEerrors_hfMatr(:,trainProf) = abs(abs(TrainLFData.MatrFFT(3:x(end), trainProf))- abs(Actual_trainabsProf(3:end,trainProf))) ./ (0.5 + abs(abs(Actual_trainabsProf(3:end,trainProf)))); % Element-wise computation of the error
        MaxMEHE_HFMatr(trainProf) = max(MEHEerrors_hfMatr(:,trainProf)); % Maximum normalized error
    end

    % % % % % % % % % % % % % % % % % % 
    if isequal(LFData.MatrFinalProfile,[])
    else
        % plot(x(3:end),abs(Pred_trainabsProf(3:x(end), trainProf)),'LineWidth',1.5,'Color', 'b' );
        MEHEerrors_hfDD(:,trainProf) = abs(abs(Pred_trainabsProf(3:x(end), trainProf))- abs(Actual_trainabsProf(3:end,trainProf))) ./ (0.5 + abs(abs(Actual_trainabsProf(3:end,trainProf)))); % Element-wise computation of the error
        MaxMEHE_HFdd(trainProf) = max(MEHEerrors_hfDD(:,trainProf)); % Maximum normalized error
    end
    % % % changing the direction (iter)
    % % err cal
    % rel_err on abs data ( new approach, a hybrid error)
    % % % MEHE Error
    % ID= DataOut.ValidIndices(ValidProf,MinValidId); Operation= idx(ID);
    % if isequal(LFData.IterFinalProfile,[])
    %     title(['FFT, Operation: ', num2str(Operation), ', HiFi vs. Data-Driven: ', num2str(RelErrabsdata_valid(ValidProf)), ...
    %         ', HF vs. Transition Matrix: ', num2str(MaxMEHE_HFMatr(ValidProf))], ...
    %         'Interpreter', 'none', 'FontSize', 12);
    %     legend('High-Fidelity', 'DataDriven','Transition Matrix','Initial Profile', 'FontSize', 12);
    % else
    %     % title(['FFT, Operation: ', num2str(Operation), ', HiFi vs. Data-Driven: ', num2str(RelErrabsdata_valid(ValidProf))], ...
    %     %     'Interpreter', 'none', 'FontSize', 12);
    %     % legend('High-Fidelity', 'Data-Driven', 'FontSize', 12);
    % end
    % grid on;
    axis tight;
    xlabel('Wave Per Turn (Lobes)',FontSize=12,Interpreter='latex');
    ylabel('Magnitude',FontSize=12,Interpreter='latex');
    title('Frequency Spectrum of the Profiles', 'FontSize', 18,Interpreter='latex');
    legend('Initial Experimental Profile','Ground truth (High Fidelity)','Prediction: Data-driven Affine Transformation', ...
            'Physics Based Affine Transformation', ...
            'Location','best', 'Interpreter', 'latex');
    % legend('Initial Experimental Profile','Ground truth (High Fidelity)','Prediction: Data-driven Affine Transformation', ...
    %         'Location','best', 'Interpreter', 'latex');
    set(gca, 'FontSize', 17);
    if Save
        % Save the plot for training data
        if isequal(LFData.IterFinalProfile,[])
            saveas(gcf, ['train Sample_FFT_HF_DD', num2str(trainProf), '_plot.png']);
            savefig(gcf,['train Sample_FFT_HF_DD', num2str(trainProf), '_plot.fig']);
        else
            saveas(gcf, ['train Sample_FFT_LF_HF_DD', num2str(trainProf), '_plot.png']);
            savefig(gcf,['train Sample_FFT_LF_HF_DD', num2str(trainProf), '_plot.fig']);
        end
    end
end

% % abs of complex data validation
% % Initialize
x= round(linspace(0,length(Actual_validabsProf(:,1)),length(Actual_validabsProf(:,1))));

MEHEerrors_hfIter=zeros(numel(x)-2,size(Actual_validabsProf,2)); MaxMEHE_HFIter=zeros(1,size(Actual_validabsProf,2));

MEHEerrors_hfMatr=zeros(numel(x)-2,size(Actual_validabsProf,2));  MaxMEHE_HFMatr=zeros(1,size(Actual_validabsProf,2));

for ValidProf=1:size(Actual_validabsProf,2)
    figure('Visible',Visibility,'Position',[10 10 700 500]);
    if isempty(LFData.LFMRR)
    else
        % % % % %
        % [~, MeasFFT] = fftg(1/SimData(ValidProf).Nwp, SimData(ValidProf).InitialProfiles(1:SimData(ValidProf).Nwp,1), 0);
        plot(x(3:end),abs(TestDataSet.xTest(3:x(end),ValidProf)),'LineWidth',1.5,'Color','k');
    end
    hold on; 
    plot(x(3:end),abs(Actual_validabsProf(3:end,ValidProf)), 'LineWidth', 3, 'Color', 'r');
    hold on;
    plot(x(3:end),abs(Pred_validabsProf(3:end,ValidProf)), '--', 'LineWidth', 3, 'Color', 'g');
    hold on;
    % ID= DataOut.ValidIndices(ValidProf,MinValidId);
    % if isequal(LFData.IterFFT,[])
    % else
    %     plot(x(3:end),abs(LFData.IterFFT(3:x(end), ID)),LineWidth=1);
    %     MEHEerrors_hfIter(:,ValidProf) = abs(abs(LFData.IterFFT(3:x(end), ID))- abs(Actual_validabsProf(3:end,ValidProf))) ./ (0.5 + abs(abs(Actual_validabsProf(3:end,ValidProf)))); % Element-wise computation of the error
    %     MaxMEHE_HFIter(ValidProf) = max(MEHEerrors_hfIter(:,ValidProf)); % Maximum normalized error
    %     hold on;
    % end 
    if isequal(LFData.MatrFinalProfile,[])
    else
        plot(x(3:end),abs(TestLFData.MatrFFT(3:x(end), ValidProf)),'LineWidth',1.5,'Color', 'b' );
        MEHEerrors_hfMatr(:,ValidProf) = abs(abs(TestLFData.MatrFFT(3:x(end), ValidProf))- abs(Actual_validabsProf(3:end,ValidProf))) ./ (0.5 + abs(abs(Actual_validabsProf(3:end,ValidProf)))); % Element-wise computation of the error
        MaxMEHE_HFMatr(ValidProf) = max(MEHEerrors_hfMatr(:,ValidProf)); % Maximum normalized error
    end

    if isequal(LFData.MatrFinalProfile,[])
    else
        % plot(x(3:end),abs(Pred_validabsProf(3:x(end), ValidProf)),'LineWidth',1.5,'Color', 'b' );
        MEHEerrors_hfDD(:,ValidProf) = abs(abs(Pred_validabsProf(3:x(end), ValidProf))- abs(Actual_validabsProf(3:end,ValidProf))) ./ (0.5 + abs(abs(Actual_trainabsProf(3:end,ValidProf)))); % Element-wise computation of the error
        MaxMEHE_HFdd(ValidProf) = max(MEHEerrors_hfDD(:,ValidProf)); % Maximum normalized error
    end
    % % % changing the direction (iter)
    % % err cal
    % rel_err on abs data ( new approach, a hybrid error)
    % % % MEHE Error
    % ID= DataOut.ValidIndices(ValidProf,MinValidId); Operation= idx(ID);
    % if isequal(LFData.IterFinalProfile,[])
    %     title(['FFT, Operation: ', num2str(Operation), ', HiFi vs. Data-Driven: ', num2str(RelErrabsdata_valid(ValidProf)), ...
    %         ', HF vs. Transition Matrix: ', num2str(MaxMEHE_HFMatr(ValidProf))], ...
    %         'Interpreter', 'none', 'FontSize', 12);
    %     legend('High-Fidelity', 'DataDriven','Transition Matrix','Initial Profile', 'FontSize', 12);
    % else
    %     % title(['FFT, Operation: ', num2str(Operation), ', HiFi vs. Data-Driven: ', num2str(RelErrabsdata_valid(ValidProf))], ...
    %     %     'Interpreter', 'none', 'FontSize', 12);
    %     % legend('High-Fidelity', 'Data-Driven', 'FontSize', 12);
    % end
    % grid on;
    axis tight;
    xlabel('Wave Per Turn (Lobes)',FontSize=12,Interpreter='latex');
    ylabel('Magnitude',FontSize=12,Interpreter='latex');
    title('Frequency Spectrum of the Profiles', 'FontSize', 18,Interpreter='latex');
    % legend('Initial Experimental Profile','Ground truth (High Fidelity)','Prediction: Data-driven Affine Transformation', ...
    %         'Physics Based Affine Transformation', ...
    %         'Location','best', 'Interpreter', 'latex');
    legend('Initial Experimental Profile','Ground truth (High Fidelity)','Prediction: Data-driven Affine Transformation', ...
            'Location','best', 'Interpreter', 'latex');
    set(gca, 'FontSize', 17);
    if Save
        % Save the plot for training data
        if isequal(LFData.IterFinalProfile,[])
            saveas(gcf, ['Validation Sample_FFT_HF_DD', num2str(ValidProf), '_plot.png']);
            savefig(gcf,['Validation Sample_FFT_HF_DD', num2str(ValidProf), '_plot.fig']);
        else
            saveas(gcf, ['Validation Sample_FFT_LF_HF_DD', num2str(ValidProf), '_plot.png']);
            savefig(gcf,['Validation Sample_FFT_LF_HF_DD', num2str(ValidProf), '_plot.fig']);
        end
    end
end


% Optionally plot validation error
% Reverse the order for plotting high to low
% % % lambdas_plot = lambdas(end:-1:1);
% % % errors_plot = errors(end:-1:1);
end

end

% Custom RMSE function
function e = compute_rmse(y_true, y_pred, isAngular)
if isAngular
    diff = angle(exp(1j * (y_true - y_pred)));  % Wrap difference to [-pi, pi]
    e = sqrt(mean(diff(:).^2));
else
    e = sqrt(mean((y_true(:) - y_pred(:)).^2));
end
end
