 % % eliminating the harmonics (operating wheel) using weights 19/07
 % % Indicator: W; 

 % % nothing changes changed cell 2 matrix for roundness with 3600 points
 % actual measurement to use for data driven model (01/21/2025)

 % % % line 98, array to cell initialization (05/05/2025)
 % % % added Root mean sqaure err for profiles reconstructio

function [Errors,DataOut] = NoFroundness(DataMat,num_harmonics_range,Plot,use_fixed_harmonics,Save,Case,W)

% Determine the number of harmonics to use (condition to deal with
% harmonic iterative or fixed number)
if use_fixed_harmonics
    harmonics_to_use = num_harmonics_range(1); % Use the first element as the fixed number of harmonics
else
    harmonics_to_use = num_harmonics_range;
end

% Analysis
% Store
% Preallocate  results for different harmonics
num_harmonics_length = length(harmonics_to_use);
num_profiles = size(DataMat(1).SpaceDomainProf, 2);
Rmse = zeros(length(num_harmonics_range), num_profiles);
Rel_ERR= zeros(length(num_harmonics_range), num_profiles);
MeasRnd = zeros(1, num_profiles);
ClipRnd = zeros(1, num_profiles);

% Iterate over different numbers of harmonics
for i = 1:num_harmonics_length
    if use_fixed_harmonics
        num_harmonics = harmonics_to_use; % Use the fixed number of harmonics
    else
        num_harmonics = harmonics_to_use(i); % Iterate through the range of harmonics
    end
    % Initialie
    CumFFTtransform_clipped=  zeros(size(DataMat(1).SpaceDomainProf, 1), num_profiles);
    CumNorm_FFTtransform_clipped=  zeros(size(DataMat(1).SpaceDomainProf, 1), num_profiles);
    Cumreduced_harmonicsGeoProf=  zeros(size(DataMat(1).SpaceDomainProf, 1), num_profiles);
    CumCentredreduced_harmonicsGeoProf=  zeros(size(DataMat(1).SpaceDomainProf, 1), num_profiles);
    AvgData= zeros(1,num_profiles);
    %  iterative over profiles
    for GeoProf= 1:size(DataMat(1).SpaceDomainProf,2)
        % Geo Prof
        Data= DataMat(1).SpaceDomainProf(:,GeoProf);
        AvgData(GeoProf)= mean(Data);
        % Actual Rnd Meas Exp
        MeasRnd(i,GeoProf)= roundness(Data);
        FFTtransform= fft(Data-mean(Data)); NsmpePoints= size(FFTtransform,1);
        % Normalize the fft (ref: GB, Rountine fftg)
        NormCoeffs = ones(NsmpePoints,1)*2/NsmpePoints; NormCoeffs(1)=NormCoeffs(1)/2;
        NormFFTtransform = NormCoeffs.*(FFTtransform);
        % Clipping at selected number of harmonics
        FFTtransform_clipped = FFTtransform;
        FFTtransform_clipped(num_harmonics+2:end-num_harmonics) = 0; % Clip harmonics to 0 excluding other harmonics
        % Normalized FFT Transform
        NormFFTtransform(num_harmonics+2:end-num_harmonics) = 0;
        % % eliminating the harmonics (operating wheel) using weights
        NormFFTtransform= NormFFTtransform.*W; 
        % Perform inverse FFT on Unnormalized FFT to get reliable estimate
        % of rouindness
        invFFTtransform_clipped = ifft(FFTtransform_clipped); % Inverse FFT
        % Verify that the result is real
        if ~isreal(invFFTtransform_clipped)
            error('RealGeoProf is not real. Check the input for proper conjugate symmetry.');
        end
       
        % Rnd Based on Clip harmonics
        ClipRnd(i,GeoProf)= roundness(invFFTtransform_clipped);
        % rms & rel err (based on norm) on roundness Meas on original
        % harmonics ex. 7200 and cLIPPED seleted harmonis
        Rmse(i, GeoProf)= rmse(MeasRnd(i,GeoProf),ClipRnd(i,GeoProf));
        % Rel_ERR(i, GeoProf)= norm(MeasRnd(i,GeoProf)-ClipRnd(i,GeoProf)) / norm(MeasRnd(i,GeoProf));
        Rel_ERR(i, GeoProf)= abs(MeasRnd(i,GeoProf) - ClipRnd(i,GeoProf)) ./ (abs(MeasRnd(i,GeoProf) ) + 0.5);
        % Cumulative FFT symmetric and anti symmetyric
        CumFFTtransform_clipped(:,GeoProf) =  FFTtransform_clipped;
        CumNorm_FFTtransform_clipped(:,GeoProf) = NormFFTtransform;
        Cumreduced_harmonicsGeoProf(:,GeoProf)= invFFTtransform_clipped;
        % Centering
        OneTurn = 2*pi*(1:NsmpePoints)'/NsmpePoints; % angular coordinate
        CosOneTurn = cos(OneTurn); SinOneTurn = sin(OneTurn);
        Profile =invFFTtransform_clipped;
        CosOffset = CosOneTurn'*Profile/NsmpePoints*2; % "X dir" offset
        SinOffset = SinOneTurn'*Profile/NsmpePoints*2; % "Y dir" offset
        % centered data:
        % meanProfile = mean(Profile);
        CumCentredreduced_harmonicsGeoProf(:,GeoProf) = (Profile)-SinOneTurn*SinOffset-CosOneTurn*CosOffset;
    end
    % FFT cLIPED transformed data
    DataOut.CumFFT_symmetric{i}= CumFFTtransform_clipped;
    % Spectra Separation Positive and Negative Freq with Normalized
    DataOut.Normalized_positive_half = CumNorm_FFTtransform_clipped(1:num_harmonics,1:end);
    DataOut.Normalized_negative_half = conj(DataOut.Normalized_positive_half(end:-1:2,1:end));
    % Geo prof uncentered and Centered
    DataOut.reduced_harmonicGeoProf = Cumreduced_harmonicsGeoProf;
    DataOut.reducedCentered_harmonicGeoProf = CumCentredreduced_harmonicsGeoProf;
    % Errors
    Errors.Rmse{i}= Rmse;
    Errors.Rel_ERR{i}= Rel_ERR;
    % Roudness
    DataOut.CumClipRnd{i}= ClipRnd;
    DataOut.CumMeasRnd{i}= MeasRnd;
    % Avg removal
    DataOut.AvgRemoved{i}= AvgData; 
    % %  variance
    % Cumvar_Profharmonic(i,:)= var(CumFFTtransform_clipped(1:end,1:end));
    % % mean spectra
    % Cumavg_spectra{i}= mean(CumFFTtransform_clipped,2);
end
% Plot RMSE values for different numbers of harmonics
% % rmsfig=figure('visible',Plot,'Position',[10 10 1700 1200]);
% % hold on;
% % % Iterate over different numbers of harmonics
% % for i = 1:length(num_harmonics_range)
% %     % Plot RMSE values for each number of harmonics
% %     plot(1:size(DataMat(1).SpaceDomainProf, 2), Rmse(i, :), '-o', 'DisplayName', sprintf('Harmonics: %d, Mean: %.2f, var: %.2f', num_harmonics_range(i), mean(Rmse(i, :)),var(Rmse(i, :))));
% % end
% % % Plot mean RMSE values
% % % Compute mean RMSE for each number of harmonics
% % hold off;
% % % Add labels and legend
% % xlabel('Spatial Domain Profiles','FontSize',12);
% % ylabel('RMSE µm (log scale)',FontSize=12);
% % title('Filtered Final Rnd Meas & Clip Harmonic (RMSE) vs. Spatial Domain Profiles',FontSize=12);
% % legend('Location', 'best',FontSize=12); grid on; axis tight;
% Set y-axis scale to logarithmic
% set(gca, 'YScale', 'log');
% Rel ERR
if Plot
relerrfig=figure('visible',Plot,'Position',[10 10 1700 1200]);
bar(mean(Rel_ERR, 2));
xlabel('Clipped at different Harmonics for Profiles Reconstruction', 'FontSize', 12);
xticks(1:numel(num_harmonics_range)); % Assuming num_harmonics_range is a vector
xticklabels(num2str(num_harmonics_range(:))); % Converting numeric values to strings
grid on;
title(['Roundness Rel Err, Case: ', Case], 'FontSize', 12); legend(['Total Profiles: ', num2str(size(DataMat(1).SpaceDomainProf,2))],FontSize=14);
ylabel('Mean Rel Err for each clipped harmonics %', 'FontSize', 12);

rmsfig=figure('visible',Plot,'Position',[10 10 1700 1200]);
bar(mean(Rmse, 2));
xlabel('Clipped at different Harmonics for Profiles Reconstruction', 'FontSize', 12);
xticks(1:numel(num_harmonics_range)); % Assuming num_harmonics_range is a vector
xticklabels(num2str(num_harmonics_range(:))); % Converting numeric values to strings
grid on;
title(['Roundness Rmse Err, Case: ', Case], 'FontSize', 12); legend(['Total Profiles: ', num2str(size(DataMat(1).SpaceDomainProf,2))],FontSize=14);
ylabel('Mean Rmse Err for each clipped harmonics', 'FontSize', 12);
% set(gca, 'YScale', 'log');
% var
% % figure('visible',Plot,'Position',[10 10 1700 1200]);
% % bar(mean(Cumvar_Profharmonic,2));
% % xlabel('Clipped at different Harmonics for Profiles Reconstruction', 'FontSize', 12);
% % xticks(1:numel(num_harmonics_range)); % Assuming num_harmonics_range is a vector
% % xticklabels(num2str(num_harmonics_range(:))); % Converting numeric values to strings
% % grid on;
% % title('Roundness Rel Err', 'FontSize', 12); legend(['Total Profiles: ', num2str(size(DataMat(1).SpaceDomainProf,2))],FontSize=14);
% % ylabel('Mean Rel Err for each clipped harmonics %', 'FontSize', 12);
% Box Plot on matrix of positive frreq
boxplotfig=figure('visible',Plot,'Position',[10 10 1700 1200]);
x= round(linspace(0,num_harmonics_range(end),num_harmonics_range(end)));
boxplot(abs(DataOut.Normalized_positive_half(3:end,1:end))',x(3:end));
grid on;
xlabel('WPT (0,1) Neglected starts from 2',FontSize=12);
% title(['For Harmonics: ', num2str(size(DataOut.Normalized_positive_half(3:end,1:end)',2))],FontSize=14);
title(['For Harmonics: ', num2str(size(DataOut.Normalized_positive_half(3:end,1:end)', 2)), ', Case=', Case], 'FontSize', 14);
grid on;
else
end 
if Save
    % Save the plot & fig file
    % RelErr_Profile reconstruction
    saveas(relerrfig, sprintf(['RelErr_Profile reconstruction at Clipped Harmonics' ...
        ',Case=%s.png'],Case),'png');
    saveas(relerrfig, sprintf(['RelErr_Profile reconstruction at Clipped Harmonics' ...
        ',Case=%s.fig'],Case),'fig');
    % % % % % % % % % % % % 
    saveas(rmsfig, sprintf(['RmseErr_Profile reconstruction at Clipped Harmonics' ...
        ',Case=%s.png'],Case),'png');
    saveas(rmsfig, sprintf(['RmseErr_Profile reconstruction at Clipped Harmonics' ...
        ',Case=%s.fig'],Case),'fig');
    % boxplotfig
    saveas(boxplotfig, sprintf(['Harmonic Contribution' ...
        ',Case=%s.png'],Case),'png');
    saveas(boxplotfig, sprintf(['Harmonic Contribution' ...
        ',Case=%s.fig'],Case),'fig');
end