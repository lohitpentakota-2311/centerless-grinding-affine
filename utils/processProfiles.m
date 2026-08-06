%  Input, WP rps, Sim Res File, Initial Profiles, Tests as per Res File, Profile
%  Curves (Sections) Ex:
% Tests= 1:22; ProfileCurves=3; InitialProfiles= CumInitialProfileTest;

% Profile Offsets Preserving

% sparkout stage contact length output 04/07

% depth of cut evaluation 12/07 Simulation, for K factor estimation (K=
% set depth of cut/set depth of cut + x) assumption Rowe

% % added C estimation (Roughing, SparkOut) and (verified with global contact length by sim) as output from Simulation (29/08) (I.e., for Roughing) (MARCO's Approach assumption of
% same material characteristics of wheel and workpiece)
% % added new input Infeed (i.e., Not a radial one)

% % (07/09) for grinding stiffness- added condition to deal with a_actual feed present in simulation
% results or not. if, not it will use theoreticcal formulation to estimate
% the value

% % (12/09) grind stiffness for spark out (i.e., ref value), (16/009)
% updated to deal with different type of tests (i.e., variable operating
% conditions)

% % (24/09) Output MRR of the Simulations, Fn of the simulations, (25/09)
% No of revolutions

% % 10/12 wrong definition of contact length for rays estimation
% (corrected)

% % added C coeff for CW! (22/03/2025) line 206

% % line 193 corrected 07/04/2025

% % Line 185 when avg contact estimation is very small (25/04/2025) for
% experimental data

function [ContactLength,ExternVar,CentredGeoProf,UnCentredGeoProf,FFT] = processProfiles(wplength,ds,dcw,WP_diameter,WP_rps,SparkOut_Time,Res,InitialProfiles,Tests,ProfileCurves,Eeq)
% Input Data
CumInitialProfileTest= InitialProfiles; 

% wplength=146; ds= 610; dcw= 355;

% wplength=50; ds= 450; dcw= 238;

% % Factor= 4.5; % a factor to consider a span of rays of contact in steady state

% Initialize output variables
% Simulation Frequency Domain
SimFinProfRealFFTTest = cell(1,length(Tests));
SimFinProfFFTTest = cell(1,length(Tests));
% Meas Frequency Domain
CumInitRealFreqProfTest = cell(1,length(Tests));
CumInitFFTProfTest = cell(1,length(Tests));
% Geo Meas, Sim
CumSimPlaneMeas_centered = cell(1,length(Tests));
CumMeas_centered = cell(1,length(Tests));
% UnCentred Profiles
CumMeas_uncentered = cell(1,length(Tests));
CumSimMeas_uncentered = cell(1,length(Tests));
% Sim FFT magnitudes, Lobes
CumSimFFTmag = cell(1,length(Tests));
CumSimlobes = cell(1,length(Tests));
% Filtered Measurements Sim, Meas
CumSimFiltRoundMeas = zeros(length(Tests), ProfileCurves);
CumMeasFiltRoundMeas = zeros(length(Tests), ProfileCurves);
% Sample points
NsamplePoints = length(Res(Tests(1)).GeoProf{1,1});
DTurn = 1/NsamplePoints;

% Contact length initialize
SimContactLength_allcurves= zeros(length(Tests),1);
SimContactLength_SparkOut= zeros(length(Tests),1);
AvgC= zeros(length(Tests),1);
AvgC_SparkOut= zeros(length(Tests),1);
AvgLg= zeros(length(Tests),1);
AvgLc= zeros(length(Tests),1);
AvgLc_SparkOut= zeros(length(Tests),1);
AvgInfeed= zeros(length(Tests),1);

% depth of cut and grinding stiffness
Avgks_allcurves= zeros(length(Tests),1);
AvgSetDepthofCut_allcurves= zeros(length(Tests),1);

% Profile Centering Parameters
OneTurn = 2*pi*(1:NsamplePoints)'/NsamplePoints;
CosOneTurn = cos(OneTurn);
SinOneTurn = sin(OneTurn);

% Loop over tests
for IT = Tests
    %  Initialize for Profile Sections
    CumInitRealFreqProf = zeros(NsamplePoints./2, 1);
    CumInitFFTProf = zeros(NsamplePoints./2, 1);
    CumSimFinFFTRealProf = zeros(NsamplePoints./2, 1);
    CumSimFinFFTProf = zeros(NsamplePoints./2, 1);
    FFTmag = zeros(ProfileCurves, 1);
    lobes = zeros(ProfileCurves, 1);
    Meas_centered  = zeros(NsamplePoints,ProfileCurves);
    SimPlaneMeas_centered  = zeros(NsamplePoints,ProfileCurves);
    SimMeas_uncentered = zeros(NsamplePoints,ProfileCurves);
    Meas_uncentered = zeros(NsamplePoints,ProfileCurves);
    CumMeasFiltRound= zeros(1,ProfileCurves);
    CumSimFiltRound= zeros(1,ProfileCurves);
    SimContactLength_Curve= zeros(ProfileCurves,1);
    SparkOut_lc= zeros(ProfileCurves,1); SparkOut_Fn= zeros(ProfileCurves,1);
    SimFnmean= zeros(ProfileCurves,1); MRRMean= zeros(ProfileCurves,1);
    ActualInfeed= zeros(ProfileCurves,1);
    lg= zeros(ProfileCurves,1);
    C_allCurves= zeros(ProfileCurves,1);
    C_allCurves_SparkOut= zeros(ProfileCurves,1);
    l_f_SparkOut= zeros(ProfileCurves,1);
    l_f= zeros(ProfileCurves,1);
    MRRAlongWplen= zeros(ProfileCurves,1); SimInfVel= zeros(ProfileCurves,1);
    SimSetDepthofCut= zeros(ProfileCurves,1);
    % ActualInfeed= zeros(ProfileCurves,1);
    SimKs= zeros(ProfileCurves,1);
    for Icurve = 1:ProfileCurves
        % Sim Profile FFT
        SimFinFFTRealProf = abs(Res(IT).FFTProf{Icurve}.value(1:end)) * 1000;
        SimFinFFTProf = Res(IT).FFTProf{Icurve}.value(1:end) * 1000;
        % Lobes Extraction
        [SortFFTmod, Isort] = sort(abs(Res(IT).FFTProf{Icurve}.value(3:end)) * 1000, 'descend');
        NlobesRid = Res(IT).FFTProf{Icurve}.Nlobes(3:end);
        NlobesRid = NlobesRid(Isort);
        for Icomponents = 1:3
            FFTmag(Icurve, Icomponents) = round(SortFFTmod(Icomponents), 3);
            lobes(Icurve, Icomponents) = NlobesRid(Icomponents);
        end
        % Sim Profile Centering
        SimMeas= -(Res(IT).GeoProf{Icurve});
        % SimMeas_uncentered(1:NsamplePoints, Icurve)= -(Res(IT).GeoProf{Icurve}-mean(Res(IT).GeoProf{Icurve}))*1000;
        SimMeas_uncentered(1:NsamplePoints, Icurve)= -(Res(IT).GeoProf{Icurve})*1000;
        ExternVar.Profileoffsets.CosOffset(IT) = CosOneTurn' * (SimMeas * 1000 ) / NsamplePoints * 2;
        ExternVar.Profileoffsets.SinOffset(IT) = SinOneTurn' * (SimMeas * 1000 ) / NsamplePoints * 2;
        SimPlaneMeas_centered(1:NsamplePoints, Icurve) = (SimMeas * 1000 ) - SinOneTurn *  ExternVar.Profileoffsets.SinOffset(IT) - CosOneTurn *ExternVar.Profileoffsets.CosOffset(IT);
        % Measured Profile Centering
        % Conditional downsampling for InitMeas as per sim output profiles
        % Check the condition to deal with the HF simulation input of
        % initial profiles
        if exist('CumInitialProfileTest', 'var') && isfield(CumInitialProfileTest{IT}(Icurve), 'DataUF')
            % Use CumInitialProfileTest{IT}(Icurve).DataUF
            dataToCheck = CumInitialProfileTest{IT}(Icurve).DataUF;
        elseif exist('CumInitialProfileTest', 'var') && ~isempty(CumInitialProfileTest{1, IT})
            % Use CumInitialProfileTest{1, IT}
            dataToCheck = CumInitialProfileTest{1, IT};
        else
            error('Neither type of data sets is available.');
        end
        if size(SimMeas, 1) < size(dataToCheck, 1)
            DownsampleFactor= size(dataToCheck, 1)/ size(SimMeas, 1);
            InitMeas = downsample(dataToCheck, DownsampleFactor);
            disp('Initial Profiles are downsampled to match with Simulation Final Data Output');
        else
            InitMeas = dataToCheck;
        end
        Measdata = InitMeas;
        % % To del with singular case of wp_rps
        if isscalar(WP_rps)
            WP_Vel_grind = WP_rps*(2*pi);
            WorkpieceDiameter= WP_diameter;
        else
            WP_Vel_grind = WP_rps(IT)*(2*pi);
            WorkpieceDiameter= WP_diameter(IT);
        end
        % Meas_uncentered(1:NsamplePoints, Icurve) = InitMeas - mean(InitMeas);
        Meas_uncentered(1:NsamplePoints, Icurve) = InitMeas;
        % offsets are used as per the simulation data
        Meas_centered(1:NsamplePoints, Icurve) = (Measdata)-SinOneTurn*ExternVar.Profileoffsets.SinOffset(IT)-CosOneTurn*ExternVar.Profileoffsets.CosOffset(IT);
        % Measured Profile FFT
        [~, MeasFFT] = fftg(DTurn, (InitMeas - mean(InitMeas)), 0);
        CumInitRealFreqProf(:, Icurve) = abs(MeasFFT);
        CumInitFFTProf(:, Icurve) = MeasFFT;
        % Sim FFT Prof
        CumSimFinFFTRealProf(:,Icurve)= SimFinFFTRealProf;
        CumSimFinFFTProf(:,Icurve)= SimFinFFTProf;
        % Filtered Roundness Measurement
        MeasFilteredProfiles = GaussianLinFilter(Meas_centered(1:NsamplePoints, Icurve), 50, 0);
        CumMeasFiltRound(:, Icurve) = roundness(MeasFilteredProfiles);
        SimFilteredProfiles = GaussianLinFilter(SimPlaneMeas_centered(1:NsamplePoints, Icurve), 50, 0);
        CumSimFiltRound(:, Icurve) = roundness(SimFilteredProfiles);
        % Contact Length Estimation Simulation
        ExternVar.No_rev(IT)= length(Res(IT).InternVar{Icurve}{Icurve}.nRad)/NsamplePoints;
        % Factor= ExternVar.No_rev(IT)/3;
        Factor= 2.5;
        ExternVar.SparkOut_Rev(IT)= ((WP_rps(IT)*SparkOut_Time(IT)));
        % CentrRev=  (ExternVar.No_rev(IT)-ExternVar.SparkOut_Rev(IT))/2; 
        CentrRev=  (ExternVar.No_rev(IT))/2; 
        % if CentrRev<=Factor
        %     Factor= Factor-CentrRev;
        % else
        % end
        if (CentrRev - Factor) < 0
            warning('Check Data: Factor is too small for Avg Contact Est at Steady state');
        end
        minLimit = 1e-6; % Minimum allowed Factor value
        while (CentrRev - Factor) < 0 && Factor > minLimit
            Factor = Factor * 0.95;
        end
        SimContactLength_Curve(Icurve,:)= (mean(Res(IT).InternVar{Icurve}{Icurve}.nRad((CentrRev-Factor)*NsamplePoints:(CentrRev+Factor)*NsamplePoints))*((WorkpieceDiameter/2)*2*pi)/(NsamplePoints));
        % Contact Length Estimation Simulation for sparkout stage
        SparkOut_Stage= ((ExternVar.No_rev(IT)-ExternVar.SparkOut_Rev(IT))*NsamplePoints);
        % ks from simulation
        SimFnmean(Icurve,:)= mean(Res(IT).InternVar{Icurve}{Icurve}.Fn((CentrRev-Factor)*NsamplePoints:(CentrRev+Factor)*NsamplePoints)); % Mean Force from the simulation
        MRRMean(Icurve,:)= mean(Res(IT).InternVar{Icurve}{Icurve}.MRR((CentrRev-Factor)*NsamplePoints:(CentrRev+Factor)*NsamplePoints)); % Sim: MRR [mm2/s]
        WorkpieceTan_velocity= WP_Vel_grind*(WorkpieceDiameter/2); % Workpiece tangential Velocity [mm/s]
        % % checking if the simulation data contains a_actual feed data or
        % not
        if isfield(Res(IT).InternVar{Icurve}{Icurve}, 'a_actual')
            ActualInfeed(Icurve,:)= mean(Res(IT).InternVar{Icurve}{Icurve}.a_actual((CentrRev-Factor)*NsamplePoints:(CentrRev+Factor)*NsamplePoints));
        else
            ActualInfeed(Icurve,:)=  MRRMean(Icurve,:)/ WorkpieceTan_velocity; % [mm2/s] / [mm/s] = [mm]
        end
        SimKs(Icurve,:)= SimFnmean(Icurve,:)/ActualInfeed(Icurve,:); % [N/mm]
        % % C estimation (29/08) (I.e., Roughing)
        de = (1/(1/WorkpieceDiameter+1/ds));  % Equivalent Diameter (d_e)
        lcs = SimContactLength_Curve(Icurve,:); % Combined contact length given by simulation (in mm)
        lg(Icurve,:) = sqrt(ActualInfeed(Icurve,:)*de); % Geometric contact length (in mm)
        Fn_prime = SimFnmean(Icurve,:); % normal force (F_n') by simulation
        % Calculate the constant C
        C_allCurves(Icurve,:) = (lcs^2 - lg(Icurve,:)^2) / ((Fn_prime) * de);
        Rr(Icurve,:)= sqrt((C_allCurves(Icurve,:)*pi*wplength*Eeq)/(8));
        % % % % % for CW wheel (lg =0)
        de_cw = (1/(1/WorkpieceDiameter+1/dcw));
        CW_C_allCurves(Icurve,:) = (lcs^2 - 0) / ((Fn_prime) * de_cw);
        Rr_CW(Icurve,:)= sqrt((CW_C_allCurves(Icurve,:)*pi*wplength*Eeq)/(8));
        if  SparkOut_Time==0
            SparkOut_lc(Icurve,:)= 0; SparkOut_Fn= 0;
            C_allCurves_SparkOut(Icurve,:) =0;
            l_f_SparkOut(Icurve,:) =0;
            SimKs_sparkout(Icurve,:)=0;
        else
            SparkOut_lc(Icurve,:)= (mean(Res(IT).InternVar{Icurve}{Icurve}.nRad(SparkOut_Stage:end))*((WorkpieceDiameter/2)*2*pi)/(NsamplePoints));
            SparkOut_Fn(Icurve,:)= mean(Res(IT).InternVar{Icurve}{Icurve}.Fn(SparkOut_Stage:end));
            % % spark out C
            C_allCurves_SparkOut(Icurve,:) = (SparkOut_lc(Icurve,:)^2 - lg(Icurve,:)^2) / (SparkOut_Fn(Icurve,:) * de);
            Rr_SparkOut(Icurve,:)= sqrt((C_allCurves_SparkOut(Icurve,:)*pi*wplength*Eeq)/(8));
            % % % % % % % % % % % % % % % CW C (No feed)
            CW_C_allCurves_SparkOut(Icurve,:) = (SparkOut_lc(Icurve,:)^2 - 0) / (SparkOut_Fn(Icurve,:) * de_cw);
            Rr_CWSparkOut(Icurve,:)= sqrt((CW_C_allCurves_SparkOut(Icurve,:)*pi*wplength*Eeq)/(8));
            % % lc deflection contact length spark out
            l_f_SparkOut(Icurve,:) = sqrt(C_allCurves_SparkOut(Icurve,:) * SparkOut_Fn(Icurve,:) * de);
            % % grind stiffness
            SimKs_sparkout(Icurve,:)= SparkOut_Fn(Icurve,:)/ActualInfeed(Icurve,:); % [N/mm]
        end
        % % lc deflection contact length
        l_f(Icurve,:) = sqrt(C_allCurves(Icurve,:) * Fn_prime * de);
        % sim depth of cut
        WorkpieceTan= WP_Vel_grind*(WorkpieceDiameter/2); % Workpiece tangential Velocity [mm/s]
        % % MRRMean(Icurve,:)= mean(Res(IT).InternVar{Icurve}{Icurve}.MRR((CentrRev-Factor)*NsamplePoints):(CentrRev+Factor)*NsamplePoints); % Sim: MRR [mm2/s]
        MRRAlongWplen(Icurve,:)= MRRMean(Icurve,:)*wplength; % Sim: MRR [mm3/s]
        % Infeed velocity
        SimInfVel(Icurve,:)= (2*MRRAlongWplen(Icurve,:))/(pi*WorkpieceDiameter);
        % Depth of cut estimation
        SimSetDepthofCut(Icurve,:)= 1/2*pi*WorkpieceDiameter*(SimInfVel(Icurve,:)/WorkpieceTan);
    end
    % Store results
    % FFT Complex and absolute for sim & meas
    SimFinProfRealFFTTest{IT} = CumSimFinFFTRealProf;
    CumInitRealFreqProfTest{IT} = CumInitRealFreqProf;
    SimFinProfFFTTest{IT} = CumSimFinFFTProf;
    CumInitFFTProfTest{IT} = CumInitFFTProf;
    % Geo Prof Uncentred and Centred
    CumSimPlaneMeas_centered{IT} = SimPlaneMeas_centered;
    CumMeas_centered{IT} = Meas_centered;
    CumMeas_uncentered{IT} = Meas_uncentered;
    CumSimMeas_uncentered{IT} = SimMeas_uncentered;
    % FFT magnitude/lobes simulation
    CumSimFFTmag{IT} = FFTmag;
    CumSimlobes{IT} = lobes;
    % Filtered roundness Measurements
    CumSimFiltRoundMeas(IT, :) = CumSimFiltRound;
    CumMeasFiltRoundMeas(IT, :) = CumMeasFiltRound;
    % Contact Length
    SimContactLength_allcurves(IT)=mean(SimContactLength_Curve);
    SimContactLength_SparkOut(IT)=mean(SparkOut_lc);
    % AvgSimContactLength= mean(nonzeros( SimContactLength));
    % Ks from simulation
    Avgks_allcurves(IT)= mean(SimKs);
    % % sparkout ks
    Avgks_allcurves_sparkout(IT)= mean(SimKs_sparkout);
    AvgSetDepthofCut_allcurves(IT)= mean(SimSetDepthofCut);
    % % C est
    AvgC(IT)=mean(C_allCurves);
    % % % CW C
    CW_AvgC(IT)=mean(CW_C_allCurves);
    % % R est
    AvgRr(IT)=mean(Rr);
    % % % CW R
    CW_AvgRr(IT)=mean(Rr_CW);
    % % % % % % 
    AvgC_SparkOut(IT)=mean(C_allCurves_SparkOut);
    % % % % % 
    CW_AvgC_SparkOut(IT)=mean(CW_C_allCurves_SparkOut);
    % % % % sprak out
    % % R est
    AvgRr_SparkOut(IT)=mean(Rr_SparkOut);
    % % % CW R
    CW_AvgRrSparkOut(IT)=mean(Rr_CWSparkOut);
    % % geo cont length
    AvgLg(IT)= mean(lg);  AvgLc(IT)= mean(l_f);
    % % spark out
    AvgLc_SparkOut(IT)= mean(l_f_SparkOut);
    % % infeed estimation should be equal to theor equal
    AvgInfeed(IT)= mean(ActualInfeed);
    % % MRR Estimation
    AvgMRR(IT)= mean(MRRAlongWplen);
    % % % Fn
    AvgFn(IT)=mean(SimFnmean);
end
% return relevant results
FFT.u= cell2mat(CumInitFFTProfTest);
FFT.y= cell2mat(SimFinProfFFTTest);
CentredGeoProf.u = cell2mat(CumMeas_centered);
CentredGeoProf.y = cell2mat(CumSimPlaneMeas_centered);
UnCentredGeoProf.u = cell2mat(CumMeas_uncentered);
UnCentredGeoProf.y = cell2mat(CumSimMeas_uncentered);
% CONTACT length and Ks from Simulation
ExternVar.SimKsAvg= (nonzeros(Avgks_allcurves));
ExternVar.SimKsAvg_sparkout= (nonzeros(Avgks_allcurves_sparkout));
ExternVar.SimSetDepthofCut=nonzeros(AvgSetDepthofCut_allcurves);
ExternVar.Infeed= nonzeros(AvgInfeed);
ExternVar.SimMRR= nonzeros(AvgMRR);
ExternVar.SimFn=nonzeros(AvgFn);
% % C estimation with roughing and spark out
ContactLength.SimContactLenAvg=(nonzeros(SimContactLength_allcurves));
ContactLength.SimContactLength_SparkOutAvg=(nonzeros(SimContactLength_SparkOut));
ContactLength.C= (nonzeros(AvgC));
ContactLength.C_SparkOut= (nonzeros(AvgC_SparkOut));
% % % roughness factor
ContactLength.Rr= (nonzeros(AvgRr));
ContactLength.Rr_SparkOut= (nonzeros(AvgRr_SparkOut));
% % 
% % % % CW 
ContactLength.CW_C= (nonzeros(CW_AvgC));
ContactLength.CW_C_SparkOut= (nonzeros(CW_AvgC_SparkOut));
% % % roughness factor
ContactLength.CW_Rr= (nonzeros(CW_AvgRr));
ContactLength.CW_RrSparkOut= (nonzeros(CW_AvgRrSparkOut));
% % cont lengths
ContactLength.Cumlg=(nonzeros(AvgLg)); ContactLength.Cumlc=(nonzeros(AvgLc));
ContactLength.Cumlc_SparkOut=(nonzeros(AvgLc_SparkOut));
end