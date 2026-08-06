     
% % Name of the excel file you want to load (Data Set Explain),
% EX. Tests= 260:282 (With Respect to the DoE);
% Profiles='Final' or 'Initial'; ProfileNo= [1:3] or Single Profile Selection;
% % DataSetType= Experimental or Simulation or Synthetic
% % ElabPar=  Settings for Figure Plotting Exp Meas & DataSet Type with
% Directory for the Meas Files

function [DataMat,Predictors,AnalysisPar,DataOut]= dataprocessing(FileNameExcel,Tests,Profiles,ProfileNo,ElabPar)

% % Input Par (DataSet type, & Experimental Meas Files Ditrectory)
DataSetType= ElabPar.DataSet;
ExperimentMeas= ElabPar.MeasFiles;

% % % Expt. Setup parameters for anova
if strcmp(DataSetType, 'Simulation') || strcmp(DataSetType, 'Experiments')
    DataDir = ExperimentMeas;
    RotAllSubDir = 'MISURAZIONE ROTONDITà TUTTE\';
    % % % % ******************** CREARE FUNZIONE A PARTE PER CARICARE DATI DA EXCEL E DEFINIRE VARIABILI
    LastDoErow = 309; % last row number in the DoE sheet of the excel data file
    PlotF = true; %
    SaveFF = false; % save figures flag
    if ~(exist([DataDir FileNameExcel],'file')==2) % file exists?
        error([[DataDir FileNameExcel] ': does not exist'])
    else
        FID=fopen([DataDir FileNameExcel]);
        if FID==-1
            error([[DataDir FileNameExcel] ': cannot be opened'])
        else
            fclose(FID);
        end
    end
    % % % % Reads the DoE data:
    OptsInp = detectImportOptions([DataDir FileNameExcel],'Sheet','DoE','VariableNamesRange','A3','DataRange','A4'); % ,'DataRange',sprintf('A4:CZ%i',LastDoErow)
    OptsInp = setvartype(OptsInp, {'provaEseguita_','NomeFileProfiloIniziale','NomeFileProfiloFinale','fareAnalisiProfilo_Iniziale', ...
        'fareAnalisiProfilo_Finale'},'char');  % or 'datetime' or 'double'
    DoET=readtable([DataDir FileNameExcel],OptsInp,'UseExcel',true); % 'UseExcel',true to import result of excel formula
    NTests = size(DoET,1);  hw=zeros(size(DoET.hw_mm));
    wp_rps= zeros(size(DoET.x_w_rps_)); CW_rpm=zeros(size(DoET.s1CW_rpm));
    Dur_s= zeros(size(DoET.Duration_s_)); InfeedD=zeros(size(DoET.DiametralInfeed_mm_));
    feedvel_mic_s= zeros(size(DoET.s1FeedVel_micron_s)); beta=zeros(size(DoET.x____));
    Nrev= zeros(size(DoET.Nrev));
    CumRoundMeas= zeros(1,3); % % Cumulative Round Measurements
    for IT = (Tests)
        if ~(isempty(DoET.fareAnalisiProfilo_Finale{IT})||strcmpi(DoET.fareAnalisiProfilo_Finale{IT},'no'))
            % % Machine Setup parameters
            ds = DoET.OWActualDiameter_DiametroMolaOperatrice(IT); % Diameter OP wheel in mm
            Predictors.hw(IT) = DoET.hw_mm(IT); % WP height (-2: stable; 20 unstable)
            Predictors.wp_rps(IT)= DoET.x_w_rps_(IT); % % % workpiece speed in rps
            Predictors.CW_rpm(IT)= DoET.s1CW_rpm(IT); % % % Control Wheel speed in rps
            Predictors.Dur_s(IT)=  DoET.Duration_s_(IT); % % % Duration
            Predictors.InfeedD(IT)= DoET.DiametralInfeed_mm_(IT)*1000; % % %  Diam.Infeed in Micron
        end
    end
    % % %  Hw, dur & InfeedD for anova
    Predictors.CW_rpm(:,~any(Predictors.CW_rpm,1))= [];
    Predictors.hw(:,~any(Predictors.hw,1))= [];
    Predictors.Dur_s(:,~any(Predictors.Dur_s,1))= [];
    Predictors.InfeedD(:,~any(Predictors.InfeedD,1))= [];
    Predictors.Dur_s(:,~any(Predictors.Dur_s,1))= [];
    % % %     Predictors.wp_rps(:,~any(Predictors.wp_rps,1))= [];
    % % % Parameters for Time Domain analysis
    timemax= 1/(min(Predictors.wp_rps(Predictors.wp_rps ~= 0))); % max time in sec
    Freqresolution= 1/timemax;
    Nolines= round(800/Freqresolution); AnalysisPar.lines= Nolines;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % Change of DataSet Simulation, Experimental or Synthetic Data Set
% % % % % SIM Data Set
if strcmp(DataSetType, 'Simulation')
    DataSetName= 'Res';
    disp('Simulation Data Set is used for analysis');
    load(DataSetName);
    % Identify empty elements in struct Res (Missing Sim Data)
    emptyElements = arrayfun(@(x) isempty(x.GeoProf), Res);
    % Remove empty elements
    Res(emptyElements) = [];
    NsamplePoints= length(Res(1).GeoProf{1,1})/2;
    % %     FFTNsamplePoints= length(Res(1).FFTProf{1,1}.value);
    FFTProf = zeros(1, numel(Res(1).GeoProf));
    numColumns = numel(Res(1).GeoProf);
    SimPlaneMeas= zeros(1,numel(Res(1).GeoProf));
    RndMeas= zeros(1,numel(Res(1).GeoProf));
    CumFreq= zeros(1,numel(Res(1).GeoProf));
    CumSTF= zeros(1,numel(Res(1).GeoProf));
    CumFilteredProfiles= zeros(1,numel(Res(1).GeoProf));
    % % % Sim Profile Centering Parameters
    OneTurn = 2*pi*(1:NsamplePoints)'/NsamplePoints;
    CosOneTurn = cos(OneTurn); SinOneTurn = sin(OneTurn);
    ElabPar.CutOffFreq= 50;
    for IT= 1:length(Res)
        for ICurve= ProfileNo
            % % %             % %                 FFTProf(1:FFTNsamplePoints, ICurve) = (Res(IT).FFTProf{ICurve}.value*1000 ...
            % % %             % %                      -mean(Res(IT).FFTProf{ICurve}.value*1000))';
            % % %             % %                 CosOffset = CosOneTurn'*(Res(IT).GeoProf{ICurve ...
            % % %             % %                     }*1000-(mean(Res(IT).GeoProf{ICurve})*1000))/NsamplePoints*2; % X offset
            % % %             % %                 SinOffset = SinOneTurn'*(Res(IT).GeoProf{ICurve ...
            % % %             % %                     }*1000-(mean(Res(IT).GeoProf{ICurve})*1000))/NsamplePoints*2; % Y offset
            % % %             % %                 SimPlaneMeas(1:NsamplePoints, ICurve) = (Res(IT).GeoProf{ICurve ...
            % % %             % %                     }*1000-(mean(Res(IT).GeoProf{ICurve})*1000))-SinOneTurn*SinOffset-CosOneTurn*CosOffset;
            % % % Downsample the Data
            SimMeas= downsample(Res(IT).GeoProf{ICurve},2);
            NsamplePoints= length(SimMeas);
            CosOffset = CosOneTurn'*(SimMeas*1000-(mean(SimMeas)*1000))/NsamplePoints*2; % X offset
            SinOffset = SinOneTurn'*(SimMeas*1000-(mean(SimMeas)*1000))/NsamplePoints*2; % Y offset
            SimPlaneMeas(1:NsamplePoints, ICurve) = (SimMeas*1000-(mean(SimMeas)*1000))-SinOneTurn*SinOffset-CosOneTurn*CosOffset;
            RndMeas(:,ICurve) = roundness(SimPlaneMeas(1:NsamplePoints, ICurve));
            [FilteredProfiles] = (GaussianLinFilter(SimPlaneMeas(1:NsamplePoints, ICurve), ElabPar.CutOffFreq, 0));
            DTurn = (0:NsamplePoints)'/NsamplePoints; % angular increment
            [~,MeasFFT] = fftg(DTurn,SimPlaneMeas(:,ICurve),0);
            FFTProf(1:Nolines,ICurve) = MeasFFT(1:Nolines);
            Predictors.wp_rps(:,~any(Predictors.wp_rps,1))= [];
            DT= 1/(Predictors.wp_rps(IT)*NsamplePoints);
            Timevect= (1:NsamplePoints)*DT;
            [Freq, STF] = fftg(Timevect,SimPlaneMeas(:,ICurve),(round((timemax-Timevect(end))/DT)));
            CumFreq(1:Nolines,ICurve) = Freq(1:Nolines);
            CumSTF(1:Nolines,ICurve) = STF(1:Nolines);
            CumFilteredProfiles(1:length(FilteredProfiles),ICurve)= FilteredProfiles;
        end
        % Calculate the column range to insert numArray
        colStart = (IT - 1) * numColumns + 1;
        colEnd = IT * numColumns;
        % Store numArray in the 2D matrix
        CumFFT(1:Nolines, colStart:colEnd) = FFTProf;
        CumSimPlaneMeas(1:NsamplePoints,colStart:colEnd)= SimPlaneMeas;
        CumRndMeas(colStart:colEnd,:)= RndMeas';
        CumProfFreq(1:Nolines,colStart:colEnd)= CumFreq;
        CumProfSTF(1:Nolines,colStart:colEnd)= CumSTF;
        CumFiltPlaneMeas(:,colStart:colEnd)= CumFilteredProfiles;
    end
    CumProfFreq(:,~any(CumProfFreq,1))= []; CumProfSTF(:,~any(CumProfSTF,1))= [];
    CumSimPlaneMeas(:,~any(CumSimPlaneMeas,1))= []; CumFFT(:,~any(CumFFT,1))= [];
    % % % Angular Domain
    ComplxProf = CumFFT; RealFreqProf= (abs(CumFFT));
    % % % Freq in WPT AngularDomain
    Freq_Xaxis= linspace(0,size(RealFreqProf,1),size(RealFreqProf,1));
    AngularDomain.Domain= 'Angular Domain';
    AngularDomain.Freq_Xaxis= repmat(Freq_Xaxis(3:end)',1,size(RealFreqProf,2));
    % abs and Complex Prof
    AngularDomain.RealFreqProf= (RealFreqProf(3:end,:));      AngularDomain.ComplxProf= ComplxProf(3:end,:);
    % % % Time Domain abs and frewq
    CumProfFreq= CumProfFreq(2:end-1,:);    CumProfSTF= (CumProfSTF(2:end-1,:));
    TimeDomain.Domain= 'Time Domain';
    TimeDomain.Freq_Xaxis=  CumProfFreq;    TimeDomain.RealFreqProf= abs(CumProfSTF);     TimeDomain.ComplxProf=  CumProfSTF;
    % % % Geo Prof for alignemnt
    AngularDomain.SpaceDomainProf= CumSimPlaneMeas; AngularDomain.RndMeas=  CumRndMeas;
    AngularDomain.FiltSpaceDomainProf= CumFiltPlaneMeas;
    TimeDomain.SpaceDomainProf= [];  TimeDomain.RndMeas= [];
    TimeDomain.FiltSpaceDomainProf= [];
    % % % % % Synthetic  Data Set
elseif strcmp(DataSetType, 'SyntheticDataSet')
    ElabPar.Sig= ['C:\Users\pentakota\C.N.R. STIIMA\DASM - Documents\Papers' ...
        '\23 Profiles Princ Comp Analysis\SyntheticDatafigure\'];
    disp('Synthetic Data Set is used for analysis');
    % % % Synthetic Data Section
    n = 100;   % % % % Number of profiles
    p = 7200;  % % % % Number of samples per profile
    % (true or false for AR(2) Noiss)
    % (Gaussian white process) for noise generation
    % Initialize variables
    y = zeros(p, n); % Matrix to store profile data
    cj1 = 0.4;  %  Auto-regressive Model coefficients cj1
    cj2 = 0.3;  %  Auto-regressive Model coefficients cj2
    % % % % Setup
    noise = zeros(p,1);  % Initialize AR(2) noise sequence
    % Amplitude variance and correlation fotr 2 & 3 Hormonics
    sigma_2_2 = 5;   % Variance of amplitude for second harmonic
    sigma_2_3 = 5;   % Variance of amplitude for third harmonic
    if isequal(ConstPhase,1)
        phase= rand()*2*pi;
        Phaserho= 1;
    end
    % % % % Generate profile data
    for j = 1:n
        if isequal(Harmonics,[2]) || isequal(Harmonics, [3])
            amplitude = mean_amp + amp_var * randn(); % % % % One Hormonic
        elseif isequal(Harmonics, [2,3]) || isequal(Harmonics, [3,16])
            amplitude = mvnrnd([4 2],[sigma_2_2^2 sigma_2_2 ...
                *sigma_2_3*rho; sigma_2_2*sigma_2_3*rho sigma_2_3^2])
            % % % % Multivariate Normal Distribution=
            % % % Co-variance Matrix + mean amplitude
            Cumamplitudes(j,:)=amplitude;
        end
        if isequal(UniformPhase,1)
            clear phase;
            phase= rand()*2*pi; Phaserho= 0; % %  Uniform Random phase
        elseif isequal(NormDistrPhase,1)
            clear phase; phase= sqrt(Phasevariance)*randn();
            Phaserho= 0;
            % % Normal Distribution Phase
        end
        % Harmonic components
        for i = (Harmonics) % % % %  hasrmonic
            % Changing Indices for the Profile Creation (Not optimal)
            if isequal(i, 2) || isequal(i, 3)
                i0=3;
            elseif isequal(i, 16)
                i0= 15;
            end
            % Formulation
            Ang = linspace(2.5, 2*pi, p); % % % % Ang Co-ordinates
            % % % % Equation (2) for 3 Hormonics
            y(:,j) = y(:,j) + (amplitude(i-i0+1) * (cos(i*Ang-phase)))';
            % Add noise to the profile: On choice white noise or AR(2)
            % Noise
            if use_white_noise
                noise = sqrt(noise_variance) * randn(p,1);  % Generate white noise
            else
                WhiteNoise = sqrt(noise_variance)*randn(p,1);
                for k = 3:p
                    noise(k) = cj1 * noise(k-1) + cj2 * noise(k-2) + WhiteNoise(k);  % Generate AR(2) noise
                end
            end
        end
        y(:,j) = y(:,j) + (noise);  % Add noise to profile j
    end
    DTurn = (0:length(y))'/length(y); % angular increment
    [~,FFT] = fftg(DTurn,y,0);
    phase = (angle(FFT));
    ComplxProf = FFT; RealFreqProf= abs(FFT);
    Freq_wpt= (0:size(RealFreqProf,1)-1);
    % % % Space Domain DataSet
    SpaceDomainProf= y;
    % % % Experimental Measures Data Set
elseif strcmp(DataSetType, 'Experiments')
    for IT = (Tests) % % % 1:length(DoET.NEsecuzione)
        fprintf('analysing roundness of test: %s\n',DoET.Id{IT});
        % % If requested, Initial Profiles
        if strcmp(Profiles, 'Initial') || strcmp(Profiles, 'Initial&Final')
            disp('Experimental Measured Initial Profiles are Selected');
            ElabPar.Sig= ['C:\Users\pentakota\C.N.R. STIIMA\DASM - Documents\Papers' ...
                '\23 Profiles Princ Comp Analysis\ExperimentDatafigure\Meas. DoE NEsecuzione 201-400 FinalProfiles\Meas. Cyl Init. Profiles DoE G201-400\'];
            ElabPar.FigureSubDir = ElabPar.Sig;
            if ~(isempty(DoET.fareAnalisiProfilo_Iniziale{IT})||strcmpi(DoET.fareAnalisiProfilo_Iniziale{IT},'no')) && ~(isempty(DoET.fareAnalisiCN{IT})||strcmpi(DoET.fareAnalisiCN{IT},'no')) && ~(isempty(DoET.fareAnalisiProfilo_Finale{IT})||strcmpi(DoET.fareAnalisiProfilo_Finale{IT},'no'))
                [~,FFT_UFInit,~,DataOutPrInit] = ElabRotondimetroAll([DataDir RotAllSubDir DoET ...
                    .NomeFileProfiloIniziale{IT} '.csv'],ElabPar.RangeMean,ElabPar);
                % %
                DataOut.InitialProfileName(IT)= DoET.NomeFileProfiloIniziale(IT);
                % DataMat.InitialProfileExist(IT)=
                % Initialize a numeric array to store the values
                FFTProf = zeros(1, numel(FFT_UFInit)); numColumns = numel(FFT_UFInit);
                PlaneMeas= zeros(1,numel(FFT_UFInit)); PlaneMeasUnCentr= zeros(1,numel(FFT_UFInit));
                FiltPlaneMeas= zeros(1,numel(FFT_UFInit));
                RndMeas= zeros(1,numel(FFT_UFInit));
                ExpCylPlaneMeas= zeros(1,numel(FFT_UFInit));
                % Extract the values from the struct and store them in the numeric array
                for i = ProfileNo
                    if IT<201 % % % Managing Tests (less than 50 contains observations with 3600)
                        % % % PlaneMeas(1:length(DataOutPrInit(i).DataUF),i) =  (DataOutPrInit(i) ...
                        % % %     .DataUF-mean(DataOutPrInit(i).DataUF));
                        PlaneMeas(1:length(DataOutPrInit(i).DataUF),i) =  (DataOutPrInit(i) ...
                            .DataUF);
                        % % % Uncentr (05/06/2025)
                        PlaneMeasUnCentr(1:length(DataOutPrInit(i).DataUF),i) =  (DataOutPrInit(i). ...
                            DataUFUnCentred);
                        % % % Filtered Data
                        % FiltPlaneMeas(1:length(DataOutPrInit(i).DataUF),i) =  (DataOutPrInit(i) ...
                        %     .Data-mean(DataOutPrInit(i).Data));
                        FiltPlaneMeas(1:length(DataOutPrInit(i).DataUF),i) =  (DataOutPrInit(i) ...
                            .Data);
                        RndMeas(:,i)= roundness(PlaneMeas(1:length(DataOutPrInit(i).DataUF),i));
                        FFTProf(1:length(FFT_UFInit(i).FFT),i) = FFT_UFInit(i).FFT;
                        % % %
                        ExpCylPlaneMeas(1:length(DataOutPrInit(i).DataUF),i)= 0;
                    else
                        warning('Data Downsampled');
                        % FFTProf(1:length(FFT_UFInit(i).FFT)/2,i) = downsample(FFT_UFInit(i).FFT-mean(FFT_UFInit(i).FFT),2);
                        PlaneMeas(1:length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered,2);
                        % % % Uncentr (05/06/2025)
                        PlaneMeasUnCentr(1:length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i) = downsample(DataOutPrInit.PlaneMeas.Profile(i) ...
                            .Profile,2);
                        % % % % % % % % % % % %
                        RndMeas(:,i)= roundness(PlaneMeas(1:length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i));
                        % % % Filtered Data
                        FiltPlaneMeas(1:length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).CenteredFilteredProfile)/2,i) =downsample(DataOutPrInit. ...
                            PlaneMeas.Profile(i).CenteredFilteredProfile,2);
                        DTurn = (0:length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2)'/length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2; % angular increment
                        [~,MeasFFT] = fftg(DTurn,PlaneMeas(:,i),0);
                        FFTProf(1:size(MeasFFT,1),i) = MeasFFT;
                        % % % Init Cycl meas
                        ExpCylPlaneMeas(1:length(DataOutPrInit. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrInit. ...
                            CylPlaneMeas.Profile(i).Profile,2);
                    end
                    % DT= 1/(Predictors.wp_rps(IT)*length(PlaneMeas));
                    % Timevect= [1:length(PlaneMeas)]*DT;
                    % [Freq, STF] = fftg(Timevect,PlaneMeas(:,i),(round((timemax-Timevect(end))/DT)));
                    % CumFreq(1:Nolines,i) = [];
                    % CumSTF(1:Nolines,i) = [];
                end
                % Calculate the column range to insert numArray
                colStart = (IT - 1) * numColumns + 1;
                colEnd = IT * numColumns;

                % Store numArray in the 2D matrix
                CumFFT(:, colStart:colEnd) = FFTProf;
                CumPlaneMeas(:,colStart:colEnd)= PlaneMeas;
                % % % % UnCentr
                CumPlaneMeasUnCentr(:,colStart:colEnd)= PlaneMeasUnCentr;
                % % % Cylindrical meas (Un processed)
                CumExpCylPlaneMeas(:,colStart:colEnd)= ExpCylPlaneMeas;
                CumRndMeas(colStart:colEnd,:)= RndMeas';
                CumFiltPlaneMeas(:,colStart:colEnd)= FiltPlaneMeas;
                CumProfFreq= [];
                CumProfSTF= [];
            end
        end
        % % If requested, Final Profiles
        if strcmp(Profiles, 'Final') || strcmp(Profiles, 'Initial&Final')
            disp('Experimental Measured Final Profiles are Selected');
            ElabPar.Sig= ['C:\Users\pentakota\C.N.R. STIIMA\DASM - Documents\Papers\23 Profiles Princ Comp Analysis\ExperimentDatafigure' ...
                '\Meas. DoE NEsecuzione 1-50 FinalProfiles\Meas. DoE 1-50 Angular & Time Domain 3 Profiles\'];
            if ~(isempty(DoET.fareAnalisiProfilo_Finale{IT})||strcmpi(DoET.fareAnalisiProfilo_Finale{IT},'no')) && ~(isempty(DoET.fareAnalisiCN{IT})||strcmpi(DoET.fareAnalisiCN{IT} ...
                    ,'no')) && ~(isempty(DoET.fareAnalisiProfilo_Iniziale{IT})||strcmpi(DoET.fareAnalisiProfilo_Iniziale{IT},'no'))
                [~,FFT_UFFin,~,DataOutPrFin] = ElabRotondimetroAll([DataDir RotAllSubDir DoET ...
                    .NomeFileProfiloFinale{IT} '.csv'],ElabPar.RangeMean,ElabPar);
                % %
                DataOut.FinalProfileName(IT)= DoET.NomeFileProfiloFinale(IT);
                % Initialize a numeric array to store the values
                FFTProf = zeros(1, numel(FFT_UFFin)); numColumns = numel(FFT_UFFin);
                PlaneMeas= zeros(1,numel(FFT_UFFin)); FiltPlaneMeas= zeros(1,numel(FFT_UFFin));
                RndMeas= zeros(1,numel(FFT_UFFin)); ExpCylPlaneMeas= zeros(1,numel(FFT_UFFin));
                CumFreq= zeros(1,numel(FFT_UFFin));
                CumSTF= zeros(1,numel(FFT_UFFin));
                % Extract the values from the struct and store them in the numeric array
                disp(['Experimental Measured Profile is Selected for the tests: ' , num2str(ProfileNo)]);
                for i = ProfileNo
                    if IT<51 % % % Managing Tests (less than 50 contains observations with 3600)
                        % % PlaneMeas(1:length(DataOutPrFin(i).DataUF),i) =  (DataOutPrFin(i) ...
                        % %     .DataUF-mean(DataOutPrFin(i).DataUF));
                        PlaneMeas(1:length(DataOutPrFin(i).DataUF),i) =  (DataOutPrFin(i) ...
                            .DataUF);
                        % % % Uncentr (05/06/2025)
                        PlaneMeasUnCentr(1:length(DataOutPrFin(i).DataUF),i) =  (DataOutPrFin(i). ...
                            DataUFUnCentred);
                        RndMeas(:,i)= roundness(PlaneMeas(1:length(DataOutPrFin(i).DataUF),i));
                        % % % Filtered Data
                        % FiltPlaneMeas(1:length(DataOutPrFin(i).DataUF),i) =  (DataOutPrFin(i) ...
                        %     .Data-mean(DataOutPrFin(i).Data));
                        FiltPlaneMeas(1:length(DataOutPrFin(i).DataUF),i) =  (DataOutPrFin(i) ...
                            .Data);
                        % % % FFT Data
                        FFTProf(1:Nolines,i) = (FFT_UFFin(i).FFT(1:Nolines)-mean(FFT_UFFin(i).FFT(1:Nolines)));
                        % % %
                        ExpCylPlaneMeas(1:length(DataOutPrFin(i).DataUF),i)= 0;
                        % % % some deviations from excel files are also
                        % including filtered deviations from monzesi
                        % included in else if condition
                    elseif IT > 51 && IT < 201
                        warning('Data Downsampled, Meas with Filetered deviations from Monzesi');
                        % % PlaneMeas(1:length(DataOutPrFin. ...
                        % %     PlaneMeas.Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrFin. ...
                        % %     PlaneMeas.Profile(i).ProfileCentered-(DataOutPrFin. ...
                        % %     PlaneMeas.Profile(i).MeanData),2);
                        PlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).ProfileCentered,2);
                        RndMeas(:,i)= roundness(PlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).ProfileCentered)/2,i));
                        % % % Filtered Data
                        FiltPlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).CenteredFilteredProfile)/2,i) =downsample(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).CenteredFilteredProfile,2);
                        DTurn = (0:length(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).ProfileCentered)/2)'/length(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).ProfileCentered)/2; % angular increment
                        [~,MeasFFT] = fftg(DTurn,PlaneMeas(:,i),0);
                        FFTProf(1:Nolines,i) = MeasFFT(1:Nolines);
                        % % % Final Cycl meas
                        ExpCylPlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas(2).Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrFin. ...
                            CylPlaneMeas(2).Profile(i).Profile,2);
                    else
                        % % % % FFTProf(1:length(FFT_UFFin(i).FFT)/2,i) = downsample(FFT_UFFin(i).FFT-mean(FFT_UFFin(i).FFT),2);
                        % % % Downsample the Data
                        warning('Data Downsampled');
                        % % PlaneMeas(1:length(DataOutPrFin. ...
                        % %     PlaneMeas.Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrFin. ...
                        % %     PlaneMeas.Profile(i).ProfileCentered-(DataOutPrFin. ...
                        % %     PlaneMeas.Profile(i).MeanData),2);
                        PlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered,2);
                        PlaneMeasUnCentr(1:length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i) = downsample(DataOutPrFin.PlaneMeas.Profile(i) ...
                            .Profile,2);
                        RndMeas(:,i)= roundness(PlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i));
                        % % % Filtered Data
                        FiltPlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).CenteredFilteredProfile)/2,i) =downsample(DataOutPrFin. ...
                            PlaneMeas.Profile(i).CenteredFilteredProfile,2);
                        DTurn = (0:length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2)'/length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2; % angular increment
                        [~,MeasFFT] = fftg(DTurn,PlaneMeas(:,i),0);
                        FFTProf(1:Nolines,i) = MeasFFT(1:Nolines);
                        % % % Final Cycl meas
                        ExpCylPlaneMeas(1:length(DataOutPrFin. ...
                            PlaneMeas.Profile(i).ProfileCentered)/2,i) =downsample(DataOutPrFin. ...
                            CylPlaneMeas.Profile(i).Profile,2);
                    end
                    DT= 1/(Predictors.wp_rps(IT)*length(PlaneMeas));
                    Timevect= (1:length(PlaneMeas))*DT;
                    [Freq, STF] = fftg(Timevect,PlaneMeas(:,i),(round((timemax-Timevect(end))/DT)));
                    CumFreq(1:Nolines,i) = Freq(1:Nolines);
                    CumSTF(1:Nolines,i) = STF(1:Nolines);
                end
                % Calculate the column range to insert numArray
                colStart = (IT - 1) * numColumns + 1;
                colEnd = IT * numColumns;

                % Store numArray in the 2D matrix
                CumFFT(:, colStart:colEnd) = FFTProf;
                CumPlaneMeas(:,colStart:colEnd)= PlaneMeas;
                % % % % UnCentr
                CumPlaneMeasUnCentr(:,colStart:colEnd)= PlaneMeasUnCentr;
                % % % Cylindrical meas (Un processed)
                CumExpCylPlaneMeas(:,colStart:colEnd)= ExpCylPlaneMeas;
                CumRndMeas(colStart:colEnd,:)= RndMeas';
                CumFiltPlaneMeas(:,colStart:colEnd)= FiltPlaneMeas;
                CumProfFreq(:,colStart:colEnd)= CumFreq;
                CumProfSTF(:,colStart:colEnd)= CumSTF;
            end
        end
    end
    % % % Angular Domain
    ComplxProf = CumFFT; RealFreqProf= (abs(CumFFT));  CumRndMeas(~any(CumRndMeas,2),:)= [];
    RealFreqProf(:,~any(RealFreqProf,1))= [];  ComplxProf(:,~any(ComplxProf,1))= [];
    AngularDomain.Domain= 'Angular Domain';
    % % % Freq in WPT Angular Domain
    Freq_Xaxis= linspace(0,size(RealFreqProf,1),size(RealFreqProf,1));
    AngularDomain.Freq_Xaxis= repmat(Freq_Xaxis(3:end)',1,size(RealFreqProf,2));
    % abs and Complex Prof Angular Domain
    AngularDomain.RealFreqProf= (RealFreqProf(3:end,:));
    AngularDomain.ComplxProf= ComplxProf(3:end,:);
    % % % Time Domain
    CumProfFreq(:,~any(CumProfFreq,1))= []; CumProfSTF(:,~any(CumProfSTF,1))= [];
    % % % Time Domain abs and frewq
    TimeDomain.Domain= 'Time Domain';
    CumProfFreq= CumProfFreq(2:end-1,:);    CumProfSTF= (CumProfSTF(2:end-1,:));
    TimeDomain.Freq_Xaxis=  CumProfFreq;   TimeDomain.RealFreqProf= abs(CumProfSTF);
    TimeDomain.ComplxProf=  CumProfSTF;
    % % % Space Domain Profile for alignment
    CumPlaneMeas(:,~any(CumPlaneMeas,1))= [];  CumPlaneMeasUnCentr(:,~any(CumPlaneMeasUnCentr,1))= [];
    CumFiltPlaneMeas(:,~any(CumFiltPlaneMeas,1))= [];
    AngularDomain.SpaceDomainProf= CumPlaneMeas;
    AngularDomain.SpaceDomainProfUnCentr= CumPlaneMeasUnCentr;
    % % Cyl Pmeas
    CumExpCylPlaneMeas(:,~any(CumExpCylPlaneMeas,1))= [];
    AngularDomain.ExpCylSpaceDomainProf= CumExpCylPlaneMeas;
    AngularDomain.FiltSpaceDomainProf= CumFiltPlaneMeas; AngularDomain.RndMeas= CumRndMeas;
    TimeDomain.SpaceDomainProf= []; TimeDomain.SpaceDomainProfUnCentr= [];
    TimeDomain.FiltSpaceDomainProf= [];  TimeDomain.RndMeas= [];
    TimeDomain.ExpCylSpaceDomainProf= [];
    % % % % % Synthetic  Data Set
elseif strcmp(DataSetType, 'RotatedProf')
    ElabPar.Sig= ['C:\Users\pentakota\C.N.R. STIIMA\DASM - Documents\Papers' ...
        '\23 Profiles Princ Comp Analysis\SyntheticDatafigure' ...
        '\Rotated Profiles Data Set\'];
    disp('Profile Rotated Data Set is used for analysis');
    % Initialize a cell array to store the profiles
    numProfiles = 100;
    profileSize = 7200;
    % Initialize a cell array to store the profiles
    profiles = zeros(profileSize,numProfiles);
    % Define a constant phase for rotation
    constantPhase = 0; % No phase shift
    % constant amplitude
    constantAmplitude = 2;
    % profile = rand(profileSize, 1) * (2 * pi / profileSize);
    i= 2; % ex. add 2nd harmonic
    % Formulation
    Ang = linspace(0, 2*pi, profileSize); % % % % Ang Co-ordinates
    % % % % Equation
    profile = (constantAmplitude * (cos(i*Ang-constantPhase))') + 1 * randn(profileSize,1) ;
    Phaserho=0; rho=0; Harmonics=[2 0];
    % Circularly shift the profile by a random amount
    for i = 1:numProfiles
        % Generate a random profile for each iteration
        shiftAmount = randi(profileSize);
        profileRotation = circshift(profile, shiftAmount);
        % Store the rotated profile
        profiles(:,i) = profileRotation;
    end
    DTurn = (0:length(profiles))'/length(profiles); % angular increment
    [~,FFT] = fftg(DTurn,profiles,0);
    phase = (angle(FFT));
    ComplxProf = FFT; RealFreqProf= abs(FFT);
    Freq_wpt= (0:size(RealFreqProf,1)-1);
    % % % Space Domain DataSet
    SpaceDomainProf= profiles;
end
% % % % ******************** END OF DATA GENERATION

DataMat= [AngularDomain,TimeDomain];

end