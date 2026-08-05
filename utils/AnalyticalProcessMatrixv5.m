% Transition Matrix Estimation
% % Cleaned & Updated to deal with stages 09/10

% % added a parameter for sparkout turns to eliminate the step (06/12)

% % to deal with non integer steps modified the iteration approach (18/12)

% % added DETACH counter for 2 stages 01/09/2025 (remaining cal are not changed)

% % Xo for roughing is modified and corrected (01/12/2025)

% % Contact filter defined for both wheels (02/17/2025). CW wheel lc takes
% into account only elastic defelction.

% % % removed filter coeff here, and K factor coeff on cw (19/03/2025)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [Estimations,InternVar] = AnalyticalProcessMatrixv5(Nwp,kFactors,hw,gamma,Feed,lcValues,dw,ds,dcw,turncalculation,ProcessSimulation,InitProfile,harmonic_range,PlotVisibility,Simulation,WP_rps,SparkoutTime,verbose,SparkOutturns,CW_lc)
DetachMent='on'; % for Iterative (On/off)
% Input Data Processinng
beta_s= asin(2*hw/(dw+ds));
betaCW=asin(2*hw/(dcw+dw));
Beta= beta_s+betaCW;
Alpha=pi/2-gamma*pi/180-(beta_s);

% % 2/05/2025
AngRes= 2*pi/Nwp; % angular position [rad]
Alpha= AngRes*round(Alpha/AngRes);
Beta= AngRes*round(Beta/AngRes);

% Grinding Contacts
K1 = sin(Beta)/sin((Alpha)+(Beta));
K2 = sin(Alpha)/sin((Alpha)+(Beta));
Col_Blade = round(Alpha/(2*pi)*Nwp);
Col_ContrWh = round((pi-Beta)/(2*pi)*Nwp);

% % % Contact Filter: Hasimoto
    function [BF,Nbuffer,Ntot,FilterOrder]= ContactFilter(EstNbuffer,Col_ContrWh,dw,lc,FilterOrder,Nwp,PlotVisibility,verbose)
        FilterOrder_FirstStage= FilterOrder;
        FreqNyq = floor(Nwp/2); % Nyquist frequency
        ZeroFreq = (pi*dw/lc); % 0 gain frequency (in N lobes per turn)
        if ZeroFreq>=FreqNyq
            error('wave filter cutting frequency too high! ZeroFreq = %f',ZeroFreq)
        end
        % Adjust FilterCoeff if ZeroFreq exceeds the Nyquist frequency
        % while ZeroFreq >= FreqNyq
        %     FilterCoeff = FilterCoeff * 0.95; % Reduce FilterCoeff by 10% (or adjust the factor as i needed)
        %     ZeroFreq = FilterCoeff * (pi * dw / lc);
        % end
        if verbose
            fprintf('Contact Filter FilterCoeff = %f, resulting in ZeroFreq = %f\n', FilterCoeff, ZeroFreq);
        end
        % % %
        FilterOrder = min(Nwp, 10+4*Nwp/ZeroFreq);
        % % % checking the size of the filter
        % Check if FilterOrder is smaller than FilterOrder_FirstStage (i.e.,
        % roughing)
        if FilterOrder <= FilterOrder_FirstStage
            % Assign FilterOrder to FilterOrder_FirstStage
            FilterOrder = FilterOrder_FirstStage;
            % disp('Filter order is set to FilterOrder_FirstStage.');
        else
            % Handle the case where FilterOrder is not smaller
            % disp('Filter order remains unchanged.');
        end
        FilterOrder = floor(FilterOrder/2)*2; % it must be even, to get an odd number of coefficients
        FiltFreqCos = [0 logspace(log10(0.05*ZeroFreq),log10(ZeroFreq),100)]; % profile sampling frequencies: focus on the cutting freq
        FiltAmplCos = 0.5*(1+(cos((lc/dw).*FiltFreqCos))); %  ampiezza =0 nella stop band (see TailFilt = .. later on)
        % add points with 0 gain, with the same "density" in frequency:
        PointsXFreq = length(FiltFreqCos)/FiltFreqCos(end);
        NZeroPoints = round((FreqNyq-FiltFreqCos(end))*PointsXFreq);
        AdditionalFreq = linspace(FiltFreqCos(end),FreqNyq,NZeroPoints); AdditionalFreq = AdditionalFreq(2:end);
        FiltFreq = [FiltFreqCos AdditionalFreq]; FiltAmpl = [FiltAmplCos zeros(1,NZeroPoints-1)];
        FiltFre01 = FiltFreq/FreqNyq; FiltFre01(end) = 1;
        % if strcmp(filter, 'Grinding wheel ')
        BF = fir2(FilterOrder,FiltFre01,FiltAmpl); ZG = sum(BF); BF= BF/ZG; % forces static gain=1
        % else
        % BF = fir2(FilterOrder,FiltFre01,FiltAmpl); BF= BF/max(abs(BF)); % forces max coeff =1
        % end
        FilterOrder = length(BF); % REM: fir2 can change the filter order
        Nbuffer = (FilterOrder-1)/2;
        Ntot = Nwp+Nbuffer; % total matrix size
        if Nwp-Col_ContrWh+EstNbuffer<Nbuffer
            error('no space for the defined contact filter')
        end
        % a buffer of CentraIndex elements must be used to avoid writing on
        % points that will be involved in subsequent filter calculations
        if isequal(PlotVisibility,'on')
            [Hfilt,Wfilt] = freqz(BF,1,100,Nwp); figure('Visible',PlotVisibility);
            ax1 = subplot(211); plot(Wfilt,abs(Hfilt)); title('contact filter response'); ylabel('modulus'); grid on;
            ax2 = subplot(212); plot(Wfilt,180/pi*unwrap(angle(Hfilt))); ylabel('phase [deg]'); grid on; xlabel('freq XX')
            linkaxes([ax1 ax2],'x');
        else
        end
    end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Matrix estimation
    function [G1step,Ntot,Nbuffer,FilterOrder]= TransitionMatrix(K1,K2,Col_Blade,Col_ContrWh,k,dw,lc,Nwp,PlotVisibility,FilterOrder,verbose,CW_lc)
        CW_ContactFilter='on';
        % % % % % % % % (02/17/2025)
        [BF,Nbuffer,Ntot,FilterOrder]= ContactFilter(0,Col_ContrWh,dw,lc,FilterOrder,Nwp,PlotVisibility,verbose);
        % % Control wheel Filter
        [BF_CW,Nbuffer_CW,~,~]= ContactFilter(Nbuffer,Col_ContrWh,dw,CW_lc,[],Nwp,PlotVisibility,verbose);
        % % % % Commented upto 152 line
        % % FilterOrder_FirstStage= FilterOrder;
        % % FreqNyq = floor(Nwp/2); % Nyquist frequency
        % % ZeroFreq =FilterCoeff*(pi*dw/lc); % 0 gain frequency (in N lobes per turn)
        % % if ZeroFreq>=FreqNyq
        % %     error('wave filter cutting frequency too high! ZeroFreq = %f',ZeroFreq)
        % % end
        % % % Adjust FilterCoeff if ZeroFreq exceeds the Nyquist frequency
        % % % while ZeroFreq >= FreqNyq
        % % %     FilterCoeff = FilterCoeff * 0.95; % Reduce FilterCoeff by 10% (or adjust the factor as i needed)
        % % %     ZeroFreq = FilterCoeff * (pi * dw / lc);
        % % % end
        % % if verbose
        % %     fprintf('Contact Filter FilterCoeff = %f, resulting in ZeroFreq = %f\n', FilterCoeff, ZeroFreq);
        % % end
        % % % % %
        % % FilterOrder = min(Nwp, 10+4*Nwp/ZeroFreq);
        % % % % % checking the size of the filter
        % % % Check if FilterOrder is smaller than FilterOrder_FirstStage (i.e.,
        % % % roughing)
        % % if FilterOrder <= FilterOrder_FirstStage
        % %     % Assign FilterOrder to FilterOrder_FirstStage
        % %     FilterOrder = FilterOrder_FirstStage;
        % %     % disp('Filter order is set to FilterOrder_FirstStage.');
        % % else
        % %     % Handle the case where FilterOrder is not smaller
        % %     % disp('Filter order remains unchanged.');
        % % end
        % % FilterOrder = floor(FilterOrder/2)*2; % it must be even, to get an odd number of coefficients
        % % FiltFreqCos = [0 logspace(log10(0.05*ZeroFreq),log10(ZeroFreq),100)]; % profile sampling frequencies: focus on the cutting freq
        % % FiltAmplCos = 0.5*(1+(cos((lc/dw).*FiltFreqCos))); %  ampiezza =0 nella stop band (see TailFilt = .. later on)
        % % % add points with 0 gain, with the same "density" in frequency:
        % % PointsXFreq = length(FiltFreqCos)/FiltFreqCos(end);
        % % NZeroPoints = round((FreqNyq-FiltFreqCos(end))*PointsXFreq);
        % % AdditionalFreq = linspace(FiltFreqCos(end),FreqNyq,NZeroPoints); AdditionalFreq = AdditionalFreq(2:end);
        % % FiltFreq = [FiltFreqCos AdditionalFreq]; FiltAmpl = [FiltAmplCos zeros(1,NZeroPoints-1)];
        % % FiltFre01 = FiltFreq/FreqNyq; FiltFre01(end) = 1;
        % % BF = fir2(FilterOrder,FiltFre01,FiltAmpl); ZG = sum(BF); BF= BF/ZG; % forces static gain=1
        % % FilterOrder = length(BF); % REM: fir2 can change the filter order
        % % Nbuffer = (FilterOrder-1)/2;
        % % Ntot = Nwp+Nbuffer; % total matrix size
        % % if Nwp-Col_ContrWh<Nbuffer
        % %     error('no space for the defined contact filter')
        % % end
        % % % a buffer of CentraIndex elements must be used to avoid writing on
        % % % points that will be involved in subsequent filter calculations
        % % if isequal(PlotVisibility,'on')
        % %     [Hfilt,Wfilt] = freqz(BF,1,100,Nwp); figure('Visible',PlotVisibility);
        % %     ax1 = subplot(211); plot(Wfilt,abs(Hfilt)); title('contact filter response'); ylabel('modulus'); grid on;
        % %     ax2 = subplot(212); plot(Wfilt,180/pi*unwrap(angle(Hfilt))); ylabel('phase [deg]'); grid on; xlabel('freq XX')
        % %     linkaxes([ax1 ax2],'x');
        % % else
        % % end
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
        % Transition Matrix and Eigenvalue Analysis
        G1step = diag(ones(1,Ntot-1),-1); % transition matrix including the wave filter
        G1step(1, Col_Blade) = k*K1;
        % % % Overlapping condition
        % Define the ranges
        % % Control wheel
        range_CW_start = Col_ContrWh - Nbuffer_CW;
        range_CW_end = Col_ContrWh + Nbuffer_CW;
        % % % grinding wheel
        range_Nbuffer_start = Nwp - Nbuffer;
        range_Nbuffer_end = Nwp + Nbuffer;
        % Check if the ranges overlap (i.e., it will throw a error)
        if range_CW_end >= range_Nbuffer_start && range_CW_start <= range_Nbuffer_end
            % If they overlap, throw an error
            error('Ranges overlap: Contact Filter Overlapped!');
        end
        % % % Filter Insertion
        % % G1step(1, Col_ContrWh) = -k*K2;
        if isequal(CW_ContactFilter,'on')
            % % % k_cw=k*kfactorcoeff_cw;
            G1step(1, Col_ContrWh+ (-Nbuffer_CW:Nbuffer_CW)) = -k*K2*BF_CW;
        else
            G1step(1, Col_ContrWh) = -k*K2;
        end
        G1step(1,Nwp + (-Nbuffer:Nbuffer)) = (1-k)*BF; % shift of the previous value + contact filter + contact compliance
        % G1step(1,Nwp) = (1-k);
    end

% % Transition Matrix different stages
% 1 step transition matrix for the different stages as required with the
% change of lc and k factor (i.e., to deal with spark out stage)
% initial filter order set to zero (ex. roughing)

% % % Initial Matrix computation with stage 1 K factor, lc
[G1step,Ntot,Nbuffer,~]= TransitionMatrix(K1,K2,Col_Blade,Col_ContrWh,kFactors(1),dw,lcValues(1),Nwp,'off',[],verbose,CW_lc);

% % if isequal(Stage,'Roughing')
EigVect1step_sparkout=[];
EigVal1step_sparkout=[]; G1step_sparkout=[];
% % else
% %     lc_sparkout= lc_sparkOut; % assumption on the rays of contact
% %     % k_sparkout= k_sparkout; % assumption (not changing)
% %     % % n = rev, e = size error calculation after roughing
% %     % n= 1;
% %     % e= (1-k)^n;
% %     % % Calculate k
% %     % k = 1 - e^(1/n);
% %     isequal(Stage,'Roughing+SparkOut');
% %     % spark out
% %     [G1step_sparkout,~,~]= TransitionMatrix(K1,K2,Col_Blade,Col_ContrWh,k_sparkout,dw,lc_sparkout,Nwp,PlotVisibility,FilterOrder);
% %     [EigVect1step_sparkout,EigVal1step_sparkout] = eig(G1step_sparkout);
% %     EigVal1step_sparkout = diag(EigVal1step_sparkout);
% %     [~,IVal] = sort(abs(angle(EigVal1step_sparkout))); % sorted by frequency
% %     EigVal1step_sparkout = EigVal1step_sparkout(IVal); EigVect1step_sparkout = EigVect1step_sparkout(:,IVal);
% % end


if isequal(Simulation,'off')
    [EigVect1step,EigVal1step] = eig(G1step);
    EigVal1step = diag(EigVal1step);
    [~,IVal] = sort(abs(angle(EigVal1step))); % sorted by frequency
    EigVal1step = EigVal1step(IVal); EigVect1step = EigVect1step(:,IVal);
    WavesXturn = angle(EigVal1step)/pi*Nwp/2;
    if isequal(PlotVisibility,'on')
        figure('Visible',PlotVisibility); plot(WavesXturn, abs(EigVal1step),'*');
        xlabel('waves per turn'); ylabel('abs(eigval)'); grid on; title('AbsEigVal1step WITH FILTER')
        % identifies the mode with the highest amplification:
        [~,IMax] = max(abs(EigVal1step));
        % Angle1step = unwrap(angle(EigVect1step(:,IMax)));
        % NWavesEigVect1step = abs(Angle1step(end)-Angle1step(1))/pi*0.5,
        figure('Visible',PlotVisibility,'Name','1stepOrthogEigenvect');
        imagesc(abs(EigVect1step'*EigVect1step)); colorbar; title('abs(EigVect1step''*EigVect1step)')
        figure('Visible',PlotVisibility,'Name','1step');
        plot3(1:length(EigVect1step(:,IMax)),real(EigVect1step(:,IMax)),imag(EigVect1step(:,IMax)));
        title('abs(EigVect1step(:,IMax))')
        xlabel('point'); ylabel('real eigvect'); zlabel('imag eigvect'); grid on;
    else
    end
else
    EigVal1step = [];
    EigVect1step = [];
    % WavesXturn = [];
end


% % Process Simulation, if Flagged (Not updated with feed term)
simulateProcess(ProcessSimulation, G1step, Nwp, Nbuffer, Col_Blade, Col_ContrWh,turncalculation,InitProfile,kFactors,dw,WP_rps,SparkoutTime,Feed,Ntot,Alpha,Beta)

% % disp(['Condition number for G1 Step: ', num2str(cond(G1step))]);
% % Function to compute final profiles with different stages
ProcessParameters.hw=hw; ProcessParameters.dcw= dcw; ProcessParameters.ds= ds; ProcessParameters.gamma= gamma;
[FinalProfiles,rInit,CentredData,~,AngularDomain,~,FeedTerm,MRRrate,RDToT_wear,TheorMRR,DataOut] = computeFinalProfile(beta_s,Beta,Alpha,dw,G1step,G1step_sparkout,EigVect1step,EigVal1step,turncalculation,Nwp,Ntot,Feed,InitProfile,beta_s,kFactors,Simulation,WP_rps,SparkoutTime,EigVect1step_sparkout,EigVal1step_sparkout,verbose,SparkOutturns,DetachMent,ProcessParameters);

if isequal(PlotVisibility,'on')
    % % To freq domain transformation
    % Initial Profiles
    Profile=rInit*1000; %  micron
    % GeoProf to Freq testing final profile
    DataMat.SpaceDomainProf= (FinalProfiles(1:end,1:end));
    W= ones(size(DataMat.SpaceDomainProf,1),1);
    [Errors,DataOut_Final] = NoFroundness(DataMat,harmonic_range,PlotVisibility,0,0,'Analytical Final Profile with Feed Effect',W);
    Estimations.FinalProfAvgRemoval= DataOut_Final.AvgRemoved{end};
    % Rel_err for harmonic reduction to range requested
    Estimations.RelErr_FinalProf_RoundnessFreqtoang= Errors.Rel_ERR{end};
    % INIT Prof
    DataMat.SpaceDomainProf= (Profile(1:end,1:end));
    [~,~] = NoFroundness(DataMat,harmonic_range,PlotVisibility,0,0,'Initial Profiles',W);
    % % % % Transformation to Freq Domain
    [~,~,H_Norm]= DFT(Ntot,0);
    TransMatrix_freqdomain= H_Norm'*(AngularDomain)*H_Norm;
    % % Eigen on freq domain
    [L_SV,Sing_val,R_SV]= eig(TransMatrix_freqdomain);
    Estimations.SingValues= Sing_val;
    Estimations.L_SV= L_SV;
    Estimations.R_SV= R_SV;
    % % % % % % % % % % % % % % % % % % % % %
    % % % % % Normalization of Profiles across fEATURES
    % [y,~,y_s]=normalize(DataOut_Final.Normalized_positive_half(1:end,1:end),2);
    % [u,~,u_s]=normalize(DataOut_Initial.Normalized_positive_half(1:end,1:end),2);
    % % % Transformation Matrix with the requested harmonic range based on
    % % final and initial on normalized Profiles
    % disp(['Condition number for inverse of u for Process Matrix Identification: ', num2str(cond(u))]);
    % u=[u;ones(1,size(u,2))];
    % TransMatrix_freqdomain= y*pinv(u);
    % TransMatrix_freqdomain_unNorm= diag(y_s)*(diag(u_s))^(-1)*TransMatrix_freqdomain;
    % % % % % % % % % % % % % % % % % %
    % Plot transformation matrix in Complex on the range requested
    figure; imagesc(abs(TransMatrix_freqdomain)); xlabel('Features','FontSize',12);
    ylabel('Features','FontSize',12); title('ProcessMatrix: abs of Complex','FontSize',12); colorbar;
    % % Eigen (Freq domAIN)
    figure; bar(abs(diag(Sing_val(2:harmonic_range(end),2:harmonic_range(end))))); xlabel('Eigen Values','FontSize',12);
    title('Freq Domain: Eigen values of Transformation Matrix from Order 1','FontSize',12); grid on; axis tight;
    % % Eigen on 1 step transition matrix (Angular Domain)
    figure; bar(abs(EigVal1step)); xlabel('Eigen Values','FontSize',12);
    title('Angular Domain: Eigen values of 1 step Transformation Matrix','FontSize',12); grid on; axis tight;
    % % % % % % % % % %
    Estimations.MaxEigVal1step = max(abs(EigVal1step));
    Estimations.WavesXturn = angle(EigVal1step)/pi*Nwp/2;
    % Estimations.TransMatrix_freqDomain_unNorm=TransMatrix_freqdomain_unNorm;
    % Centred Geo Profiles are based on the reduced harmonic range, constructed profiles
    Estimations.Harmonicsreduced_CentredData= DataOut_Final.reducedCentered_harmonicGeoProf;
    InternVar.TransMatrix_freqDomain= TransMatrix_freqdomain;
else
end

% figure('visible',PlotVisibility); bar(abs(EigValnstep)); xlabel('Eigen Values','FontSize',12);
% title('Angular Domain: Eigen values of N step Transformation Matrix','FontSize',12); grid on; axis tight;

% Simulation (If requested) not updated
    function simulateProcess(ProcessSimulation, G1step, Nwp, Nbuffer, Col_Blade, Col_ContrWh,turns,InitProfile,kFactors,dw,WP_rps,SparkoutTime,Feed,Ntot,~,~)
        % Function to simulate the process and generate a video of the evolving profile
        if ProcessSimulation
            No_Stages = size(nonzeros(kFactors ~= 0),1);
            if No_Stages==1
                SparkOutsteps=0; else
            end
            Nsteps = turns*Nwp;
            % WP_length= 146;
            DeltaW= pi*dw/Nwp;
            DT= (2*pi)/Nwp/(WP_rps*2*pi);
            % % SparkOut Steps
            SparkOutsteps= ceil(SparkoutTime*WP_rps)*Nwp;
            % % fprintf('G1step condition number: %f',cond(G1step))
            DeltaX_Stage1 = zeros(Ntot,1); X0 = zeros(Ntot,1); DeltaX_Stage2 = zeros(Ntot,1);
            DeltaX_Stage1(1) =(Feed(1))/Nwp; DeltaX_Stage2(1) = (Feed(2))/Nwp;
            % Test simulation with the transition matrix with filter
            Process = G1step; % G single step, G1turn
            InitProfile=[InitProfile(1:end,:);InitProfile(1:Nbuffer,:)]/1000;
            ro=InitProfile;
            % % Initial Position by Stage 1 K factor
            X0(1)= ((ro(Nwp))-Process(1,:)*ro)/kFactors(1);
            % Define DeltaX and X0 vectors before the loop to identify the stages
            Current_DeltaX = [DeltaX_Stage1, DeltaX_Stage2];  % Set DeltaX for each stage, 0 for the second stage
            Current_Steps= [Nsteps, SparkOutsteps];
            % Create a VideoWriter object for the entire process
            writerObj = VideoWriter('LFProfileevolution.mp4', 'MPEG-4');
            writerObj.FrameRate = 10; % Adjust the frame rate as needed
            open(writerObj);
            % %
            PreviousProfile = InitProfile;
            % rlim([minRadius, maxRadius]);
            for Stages= 1: No_Stages
                DeltaX =  Current_DeltaX(:,Stages);
                Steps= Current_Steps(Stages);
                k=kFactors(Stages);
                % Adjust the starting index for stage 2
                if Stages == 1
                    startStep = 1;        % Stage 1 starts at step
                else
                    startStep = Nsteps + 1;  % Stage 2 starts
                    disp(['Estimating SparkOut Transformation matrix for the turns: ', num2str(Steps/Nwp)]);
                    X0_sparkout=kFactors(Stages-1)/kFactors(Stages)*(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
                    % % X0_sparkout=(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
                    X0= X0_sparkout;
                    % G1step=G1step_sparkout;
                end
                ActualProfile = PreviousProfile;
                AveRad = dw/2-mean(PreviousProfile);
                % maxRadius = AveRad * 1.6; %
                HF = figure;
                % set(HF, 'Position', get(0, 'ScreenSize')); % Maximizes the figure window
                % rlim([minRadius, maxRadius]); % Set radial limits
                HP = polarplot(2 * pi / Ntot:2 * pi / Ntot:2 * pi, AveRad - 1e3*ActualProfile);
                hold on;
                % Add support blade and control wheel positions (these are static overlays)
                % Plot the process positions (translating them to angular positions):
                SizeProcess = length(Process);
                theta_Blade = 2 * pi * Col_Blade / SizeProcess; % Convert to angular position
                theta_ContrWh = 2 * pi * Col_ContrWh / SizeProcess; % Convert to angular position
                % theta_Nwp = 2 * pi * Nwp / SizeProcess; % Convert to angular position
                % Plot lines or points representing these positions:
                polarplot([0,0], [AveRad*0.8, AveRad * 1.1],'y', 'LineWidth', 2);
                polarplot([theta_Blade, theta_Blade], [AveRad*0.8, AveRad * 1.1], 'r', 'LineWidth', 2); % Support Blade position
                polarplot([theta_ContrWh, theta_ContrWh], [AveRad*0.8, AveRad * 1.1], 'b', 'LineWidth', 2); % Control Wheel position
                % polarplot([theta_Nwp, theta_Nwp], [AveRad*0.8, AveRad * 1.1], 'g', 'LineWidth', 2); % Contact Filter position
                text(theta_Blade, AveRad * 1.1 + 0.02 * AveRad, 'Support Blade', 'Color', 'r', 'FontSize', 8, 'HorizontalAlignment', 'center');
                text(theta_ContrWh, AveRad * 1.1 + 0.02 * AveRad, 'Control Wheel', 'Color', 'b', 'FontSize', 8, 'HorizontalAlignment', 'center');
                text(0, AveRad * 1.1 + 0.02 * AveRad, 'Grinding Wheel', 'Color', 'y', 'FontSize', 8, 'HorizontalAlignment', 'center');
                % HT=text(3/2*pi*0.9,10,['Rev.=',num2str(0)]);
                % text(theta_Nwp, AveRad * 1.1 + 0.02 * AveRad, 'Contact Filter', 'Color', 'g', 'FontSize', 8, 'HorizontalAlignment', 'center');
                title('Profile Evolution in Polar Plot with Support Blade & Control Wheel Position');
                % maxRadius =  (dw / 2 + max(abs(InitProfile))); % Define a maximum radius limit
                % rlim([0, maxRadius]); % Set fixed limits for the plot
                % rlim([minRadius, maxRadius]); % Ensure the radial limits are maintained
                % HT = text(3 / 2 * pi, 0.9 * dw / 2, ['Rev. = ', num2str(0)]);
                % rlim([0, 1.2]);
                % HT = text(3 / 2 * pi, 0.9 * dw / 2, ['Rev. = ', num2str(0)]);
                % % SizeProcess = length(Process);
                % % % plots the evolving workpiece profile
                % % subplot('Position', [0.1, 0.25, 0.8, 0.7]) ; % [left bottom width height]
                % % HP = plot(PreviousProfile); title(sprintf('profile + %d points filter buffer',Nbuffer)); ...
                % %     xlabel('radial position + buffer'); grid on; legend('actual'); % legend('previous','actual')
                % % % AX = axis; xlim([1,SizeProcess]);
                % % maxProfileValue = max(abs(PreviousProfile)) * 1.2; % Use a multiplier to ensure enough space
                % % ylim([-maxProfileValue, maxProfileValue])
                % % % ylim([-1.1*AX(4) 1.1*AX(4)]); % leaves space for the negative peaks generated by the contact with the control wheel
                % % % plots, as a reference, the first row ofthe transition matrix:
                % % subplot('Position',[0.1, 0.1, 0.8, 0.15]) ; % [left bottom width height]
                % % imagesc(Process(1,:)); title('Support Blade|Control Wheel|Contact filter locations'); xlim([1,SizeProcess]);
                % % MaxProcess1 = max(abs(Process(1,:)));
                % % hold on; image(Col_Blade +(-1:1), [0.5 1 1.5], MaxProcess1*1.1*[1 1 1]);
                % % image(Col_ContrWh +(-1:1), [0.5 1 1.5], MaxProcess1*1.1*[1 1 1]);
                % % image(Nwp +[-Nbuffer 0 Nbuffer], [0.5 1 1.5], MaxProcess1*1.1*[1 1 1]);
                % % hold off
                % % % Create a VideoWriter object
                % % writerObj = VideoWriter('process_animation.mp4', 'MPEG-4');
                % % writerObj.FrameRate = 10; % Adjust the frame rate as needed
                % % open(writerObj);
                % Loop over steps (Nd*rev) Roughing stage
                for II = startStep:(startStep + Steps - 1)
                    ActualProfile = Process*PreviousProfile+k*(X0+(II-1)*DeltaX);
                    IntersFlag=((ActualProfile(1)-PreviousProfile(Nwp))*DeltaW)/DT;
                    % IntersFlag
                    if  IntersFlag(1)>0
                        % rSim1_rem(II)=(rSim(1)-rSimTemp(Nwp))*WP_length*pi*dw*2;
                        rSim1_rem(II)=((ActualProfile(1)-PreviousProfile(Nwp))*DeltaW)/DT;
                        % % cumulated wear
                        rSim1(II)= ActualProfile(1);
                    else
                        % % cumulated wear
                        rSim1(II)=PreviousProfile(1);
                        rSimTemp= circshift(PreviousProfile,1);
                        rSimTemp(1)= PreviousProfile(Nwp);
                        rSim1_rem(II)=0;
                        % FeedTerm(II)=0;
                    end
                    % rSim = rSimTemp;
                    % HP(1).YData = PreviousProfile; HP(2).YData = ActualProfile;
                    % HP(1).YData = ActualProfile;
                    % Update the polar plot data
                    % Clamp ActualProfile values to stay within limits
                    % ActualProfile = max(min(ActualProfile, maxRadius), 0);
                    % if mod(II,Ntot*ResDisplay/360)==0 % Figure update every Nd/ResDisplay revolution
                    AveRad = dw / 2 - mean(ActualProfile);
                    Deviations=ActualProfile-mean(ActualProfile);
                    % Dev=circshift(Dev,-CI(theta));
                    % HP.ThetaData = 2 * pi / Ntot:2 * pi / Ntot:2 * pi;
                    set(HP,'RData' ,AveRad-Deviations*1e3);
                    % set(HT,'String',['Rev.=',num2str())]); % number of revolution
                    % rlim([0, maxRadius]);
                    % Update the text information
                    % HT.String = ['Rev. = ', num2str(II / Nwp)];
                    drawnow;
                    % Write the current frame to the video file
                    frame = getframe(HF);
                    writeVideo(writerObj, frame);
                    PreviousProfile = ActualProfile;
                    % end
                end
                % Store the last profile after the final iteration
                lastProfile = ActualProfile;
            end
            % Close the VideoWriter object
            close(writerObj);
            % Save the last profile to a .mat file
            save('Profilefinalemat', 'lastProfile');
        end
    end

% Return relevant results
Estimations.Det_G = det(G1step); % is equal to +/- 1-k, from Laplace extension to compute the determinant
% % Geo Profiles with just centering actual depends on the requested sample
% points of workpiece
Estimations.AnalyticalProfileCentredData=CentredData;
Estimations.AnalyticalProfileUnCentredData= FinalProfiles;  % micron raw profiles without any type of changes of output
% G1 step for N revolutions
% Estimations.G1Step_N= G1step_n;
% G1 step for SparkOut revolutions
% % Estimations.G1Step__SparkOutStep= G1step_SparkOutStage;
% % Eigen vectors on N step mtrix (sorted by frequency)
% Estimations.Eigvectors_nStep= EigVectnstep;
% % % % RDToT (to comapare with HF Simulation)
InternVar.RDToT= RDToT_wear;
InternVar.MRR= MRRrate;
InternVar.MRRTheor=TheorMRR;
InternVar.FeedVector= FeedTerm;
InternVar.HFDataOut= DataOut;
end


% Compute final profiles computation function
function [FinalProfiles,rInit,CentredData,G1step_n,AngularDomain,G1step_SparkOutStage,FeedTerm,MRR,RDToT,TheorMRR,DataOut] = computeFinalProfile(beta_s,Beta,Alpha,dw,G1step,~,EigVect1step,EigVal1step,turns,Nwp,Ntot,Feed,InitProfile,~,kFactors,Simulation,WP_rps,~,~,~,verbose,SparkOutturns,DetachMent,ProcessParameters)
Nsteps = turns*Nwp; WP_length= 146;
DeltaW= pi*dw/Nwp;
DT= (2*pi)/Nwp/(WP_rps*2*pi);
% % SparkOut Steps
% SparkOutsteps= ceil(SparkoutTime*WP_rps)*Nwp;
SparkOutsteps= SparkOutturns*Nwp;
% InvG1step = inv(G1step);
% % fprintf('G1step condition number: %f',cond(G1step))
DeltaX_Stage1 = zeros(Ntot,1); X0 = zeros(Ntot,1); DeltaX_Stage2 = zeros(Ntot,1);
% % To compute the ToT diameteral removal for the MRR Computation
DiamRemoval = (Feed(1)/(1-sin(Beta)/sin((Beta+Alpha))+sin(Alpha)/sin((Beta+Alpha)))*2)/Nwp;
% % DeltaX
% DeltaX(1) = DiamRemoval;
DeltaX_Stage1(1) =(Feed(1)*cos(beta_s))/Nwp; DeltaX_Stage2(1) = (Feed(2))/Nwp;
% projecting with respect to the simulation feed signal
% horizontal direction material removel
% If Intial Profile is given
if isempty(InitProfile)
    NoiseAmpl = 1e-12; % 5e-2;
    rInit = NoiseAmpl*(rand(Ntot,1)-0.5); % white noise initial profile
else
    rInit= InitProfile;
    Nbuffer= Ntot-Nwp;
    rInit= [rInit(1:end,:);rInit(1:Nbuffer,:)]/1000; % to mm
end

% % % % % % % % % % % % % % % % % % %
% % % Loop over initial Profiles
% Initialize
FinalProfiles= zeros(Ntot,size(rInit,2));
No_Stages = size(nonzeros(kFactors ~= 0),1);
if No_Stages==1
    SparkOutsteps=0; else
end
rSim1_rem= zeros((Nsteps+SparkOutsteps),1);
rSim1= zeros((Nsteps+SparkOutsteps),1);
CumrSim_RDToT= zeros((Nsteps+SparkOutsteps),size(rInit,2));
CumrSim_rem= zeros((Nsteps+SparkOutsteps),size(rInit,2));

% Profile Centering parameters
NPointsPlaneMeas = (Ntot);
OneTurn = 2*pi*(1:NPointsPlaneMeas)'/NPointsPlaneMeas; % angular coordinate
CosOneTurn = cos(OneTurn); SinOneTurn = sin(OneTurn);

% Precompute constants for Analytical Estimations using geometric series!!!!
if isequal(Simulation,'off')
    InvG1step = G1step^(-1);
    eyeNtot = eye(Ntot);
    InvG1step2 = InvG1step^2;
    EigVal1stepNsteps = diag(EigVal1step.^Nsteps);
    EigVal1stepNstepsMinus1 = diag(EigVal1step.^(Nsteps-1));
    invTerm = ((eyeNtot - InvG1step)^2)^(-1);
    eigTerm = EigVect1step * EigVal1stepNstepsMinus1 / EigVect1step;
    constInvG1stepTerm = Nsteps * InvG1step2 - (Nsteps + 1) * InvG1step + eigTerm;
    % % N step using G1STEP Computation
    G1step_n= EigVect1step * EigVal1stepNsteps / EigVect1step;
    % % 12/09/2025
    % G1step_n= diag(diag(G1step_n)); 
    AngularDomain= G1step_n;
    % %
    X0termseries=((eye(Ntot)-G1step)^(-1))*(eye(Ntot)-G1step_n);
else
    G1step_n=[];
    AngularDomain= [];
end

% % Loop Over Initial Profiles
for InitProf= 1:size(rInit,2)
    % as per hf (processing initial profiles)
    FFTInitPro=fft(-rInit(:,InitProf)); %  FFT
    LFlg=Nwp-length(rInit);
    if mod(length(rInit),2)==0
        if LFlg>0
            ZRs=zeros((Nd-length(rInit))/2,1);
            FFTInitProDense=[FFTInitPro(1:end/2);ZRs;FFTInitPro(end/2+1);ZRs;FFTInitPro(end/2+2:end)];
            RD=ifft(FFTInitProDense)*Ntot/length(rInit);
        elseif LFlg==0
            RD=-rInit;
        elseif LFlg<0
            FFTInitProDense=[FFTInitPro(1:Ntot/2);FFTInitPro(end/2+1);FFTInitPro(end-Ntot/2+2:end)];
            RD=ifft(FFTInitProDense)*Ntot/length(rInit);
        end
    else
        % Odd number of points (beacause of buffer points)
        if LFlg > 0
            ZRs = zeros((Nd - length(rInit)) / 2, 1);
            % For odd, we pad differently because the symmetry is handled differently
            FFTInitProDense = [FFTInitPro(1:(end-1)/2); ZRs; FFTInitPro((end+1)/2); ZRs; FFTInitPro((end+3)/2:end)];
            RD = ifft(FFTInitProDense) * Ntot / length(rInit);
        elseif LFlg == 0
            RD = -rInit; % No change needed
        elseif LFlg < 0
            % For truncating the odd-length FFT
            FFTInitProDense = [FFTInitPro(1:(Ntot-1)/2); FFTInitPro((end+1)/2); FFTInitPro(end-(Ntot-1)/2+1:end)];
            RD = ifft(FFTInitProDense) * Ntot / length(rInit);
        end
    end
    % % % % % % % % % % % % % % % % % %
    rSim = RD;
    % Condition to deal to work with simulation or not (time consuming!!)
    ro=rSim;
    % % Initial Position by Stage 1 K factor
    % X0(1)= ((ro(Nwp))-G1step(1,:)*ro)/(kFactors(1));
    % % New X0 estimation (15/05/2025) GB
    % % G1stepMR = G1step(1,:);
    % % G1stepMR((Ntot-1-2*Nbuffer):end) = kFactors(1)/(1-kFactors(1))*G1stepMR((Ntot-1-2*Nbuffer):end); % the coeff infront of the GW filter is changed from (1-k) to k, to remove the filtered value around wp
    % % X0(1)= -G1stepMR*ro/kFactors(1);
    % % % Xo excatly equal to the HF Model (20/05/2025)
    FlagInitProf=1; 
    % % stage1 & stage 2 revolutions
    MachProc.InfeedD=Feed(1)*1000; % % Feed is in mm
    MachProc.DynStruct.OmegaW= (WP_rps*2*pi);
    InitialProfile.DataUF=flip(InitProfile);
    X0(1)= HF_X0Estimation((InitialProfile),MachProc,Nsteps/length(InitProfile),SparkOutsteps/length(InitProfile),dw,FlagInitProf,ProcessParameters);
    % % % % % % % % % % % % % % % 
    DataOut.Xopos=X0(1);
    % Define DeltaX and X0 vectors before the loop to identify the stages
    Current_DeltaX = [DeltaX_Stage1, DeltaX_Stage2];  % Set DeltaX for each stage, 0 for the second stage
    Current_Steps= [Nsteps, SparkOutsteps];
    % transformation matrix for the requested turns and with Final Profile
    % estimation with the effect of feed
    if isequal(Simulation,'on')
        if verbose
            disp(['Staring LF Iterative for the turns: ', num2str(Nsteps/Nwp)]);
        end
        % % Initialize detach counter
        DetachCounterStage1=0; DetachCounterStage2=0;
        for Stages= 1: No_Stages
            DeltaX =  Current_DeltaX(:,Stages);
            Steps= Current_Steps(Stages);
            k=kFactors(Stages);
            % Adjust the starting index for stage 2
            if Stages == 1
                startStep = 1;        % Stage 1 starts at step 1
            else
                startStep = Nsteps + 1;  % Stage 2 starts
                if verbose
                    disp(['Staring LF Iterative for the spark out turns: ', num2str(Steps/Nwp)]);
                end
                % % X0_sparkout=((kFactors(Stages-1)*(kFactors(Stages-1)*kfactorcoeff_cw))/(kFactors(Stages)*kFactors(Stages-1)*kfactorcoeff_cw))*(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
                % % not accounting for control wheel k factor
                X0_sparkout=kFactors(Stages-1)/(kFactors(Stages))*(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
                % % X0_sparkout=(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
                X0= X0_sparkout;
                % G1step=G1step_sparkout;
            end
            % Loop over steps (Nd*rev) Roughing stage
            for II = startStep:(startStep + Steps - 1)
                % rSimTemp = G1step*rSim+k*(X0+(II-1)*DeltaX);
                rSim1Temp = G1step(1,:)*rSim+k*(X0(1)+(II-1)*DeltaX(1));
                IntersFlag=((rSim1Temp-rSim(Nwp))*DeltaW)/DT;
                % IntersFlag
                if strcmp(DetachMent, 'on')
                    Value = IntersFlag(1) > 0; % on condition
                else
                    Value = (IntersFlag(1) == 1); % Off condition, negative removal
                end
                % if IntersFlag(1)>0
                if Value
                    % rSim1_rem(II)=(rSim(1)-rSimTemp(Nwp))*WP_length*pi*dw*2;
                    rSim1_rem(II)=((rSim1Temp-rSim(Nwp))*DeltaW)/DT;
                    % % cumulated wear
                    rSim1(II)= rSim1Temp(1);
                else
                    % % cumulated wear
                    rSim1(II)=rSim(1);
                    rSim1Temp= rSim(Nwp);
                    rSim1_rem(II)=0;
                    % FeedTerm(II)=0;
                    % Increment the respective counter based on the stage
                    if Stages == 1
                        DetachCounterStage1 = DetachCounterStage1 + 1;
                    else
                        DetachCounterStage2 = DetachCounterStage2 + 1;
                    end
                end
                rSim = circshift(rSim,1);
                rSim(1)=rSim1Temp;
            end
            % % phase offfset to deal with non integer steps
            % if mod(Steps, Nwp) ~= 0
            % phaseOffset = (Current_Steps(1)+Steps)/Nwp; % for residual phase offset
            % rSim = circshift(rSim(1:Nwp), -round(phaseOffset*2*pi)); % Align phase
            % end
        end
        % if isequal(Stage,'Roughing')
        %     G1step_SparkOutStage=[];
        %     % FeedTerm=[];
        % else
        %     X0=(X0(1)+(II-1)*DeltaX(1));
        %     % % % Second stage profile given by roughing stage
        %     isequal(Stage,'Roughing+SparkOut');
        %     % Spark out phase
        %     for JJ = 1:SparkOutsteps
        %         % rSimTemp = G1step * rSim + k * x0;
        %         % rSim1(II)= rSimTemp(1);
        %         % rSim1_rem(II)= rSim(1)-rSimTemp(Nwp);
        %         % rSim = rSimTemp;
        %         % % % variables are empty, if iteration approach is used
        %         % % (not computed)
        %         rSimTemp = G1step*rSim+k_sparkout*(X0);
        %         IntersFlag=((rSimTemp(1)-rSim(Nwp))*DeltaW)/DT;
        %         if  IntersFlag(1)>0
        %             FeedTerm(II)= k*(X0(1));
        %             rSim1_rem(II)=((rSimTemp(1)-rSim(Nwp))*DeltaW)/DT;
        %             rSim1(II)= rSimTemp(1);
        %         else
        %             % % cumulated wear
        %             rSim1(II)=rSim(1);
        %             rSimTemp= circshift(rSim,1);
        %             rSimTemp(1)= rSim(Nwp);
        %             rSim1_rem(II)=0; FeedTerm=0;
        %         end
        %         rSim = rSimTemp;
        %     end
        G1step_SparkOutStage=[];
        FeedTerm=[];
        % end
        % % Store r sim Var in micron
        CumrSim= rSim*1000;
        FinalProfiles(:, InitProf)=CumrSim;
        % % Removal & RDToT to compare with High Fidelity model (16/09)
        CumrSim_RDToT(:,InitProf)= rSim1;
        CumrSim_rem(:,InitProf)= (rSim1_rem);
    else
        if verbose
            disp(['Staring LF using Geometeric Series Formulations for the turns: ', num2str(Nsteps/Nwp)]);
        end
        % DeltaX=DeltaX_Stage1;
        % % Analytical Estimation
        % % rAna= EigVect1step*diag(EigVal1step.^Nsteps)/EigVect1step*rInit(:,InitProf)+inv((eye(Ntot)-InvG1step)^2)*(Nsteps*InvG1step^2-(Nsteps+1)*InvG1step+EigVect1step*diag(EigVal1step.^(Nsteps-1))/EigVect1step)*DeltaX;
        % condition to deal with stage rough or spark out for now!!!!!!
        % 14/06/2024
        Current_DeltaX = [DeltaX_Stage1, DeltaX_Stage2];
        DeltaX= Current_DeltaX(:,1);
        % % % closed feed term 12/12
        % % % FeedTerm= invTerm * constInvG1stepTerm * DeltaX * kFactors(1);
        FeedTerm=[]; DetachCounterStage1=0; DetachCounterStage2=0;
        % rAna_rough = G1step_n * rSim ...
        %     +  kFactors(1)*(X0(1)+invTerm * constInvG1stepTerm * DeltaX);
        % % 01/12/2025
        rAna_rough = G1step_n * rSim ...
            +  kFactors(1)*(X0termseries*X0 +invTerm * constInvG1stepTerm * DeltaX);
        % if isequal(Stage,'Roughing')
        if No_Stages==1
            rAna= rAna_rough;
            G1step_SparkOutStage=[];
        else
            % X0_sparkout=kFactors(1)/kFactors(2)*(X0+(Nsteps + 1 -1)*Current_DeltaX(2));
            % % % X0_sparkout=(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
            % X0= X0_sparkout;
            % % isequal(Stage,'Roughing+SparkOut');
            % % % Spark Out Phase
            X0_sparkout=((kFactors(No_Stages-1))/(kFactors(No_Stages)))*(X0+(Nsteps + 1 -1)*Current_DeltaX(:,No_Stages-1));
            % % X0_sparkout=(X0+(Nsteps + 1 -1)*Current_DeltaX(:,Stages-1));
            X0= X0_sparkout;
            SparkOut_EigVal1stepNsteps = diag(EigVal1step.^SparkOutsteps);
            % % SparkOut Using G1Step
            G1step_SparkOutStage= EigVect1step * SparkOut_EigVal1stepNsteps / EigVect1step;
            rAna= G1step_SparkOutStage * rAna_rough + ((eye(Ntot)-G1step)^(-1))*(eye(Ntot)-(G1step_SparkOutStage))*kFactors(2)*X0;
            % % % phase offfset to deal with non integer steps
            % if mod(SparkOutsteps, Nwp) ~= 0
            %     phaseOffset = mod(SparkOutsteps, Nwp); % for residual phase offset
            %     rAna = circshift(real(rAna), -phaseOffset); % Align phase
            % end
        end
        % % % % % % % % % % % % % % % % % end
        CumrAna = rAna * 1000; % Store Var in micron
        FinalProfiles(:, InitProf)= CumrAna;
        % % If, using the formulations (a faster computation)
        % the wear and MRR cannot be computed
        CumrSim_RDToT(:, InitProf)=0;
        CumrSim_rem(:, InitProf)=0;
    end
    % % % Centering 14/06/2024
    % (actual profiles without harmonic reduction) taking real part as imag
    % part is almost zero (round off errorrs)
    Profile = (FinalProfiles(:, InitProf));
    CosOffset = CosOneTurn'*Profile/Nwp*2; % "X dir" offset
    SinOffset = SinOneTurn'*Profile/Nwp*2; % "Y dir" offset
    % centered data:
    % meanAnaProfile = mean(anaProfile);
    CentredData.GeoProfiles_Ana(:,InitProf) = (Profile)-SinOneTurn*SinOffset-CosOneTurn*CosOffset;
    % % % % Cumulated Wear (should be comparable with RDToT Of HF model)
    RDToT=  CumrSim_RDToT(:, InitProf);
    MRR= CumrSim_rem(:, InitProf);
    % % % Angular domain G1
    % % % % if isequal(Simulation,'on')
    % % % %     AngularDomain= G1step_n;
    % % % % else
    % % % %     AngularDomain= G1step_n+FeedTerm;
    % % % % end
end
% % % formula MRR
duration= (turns/WP_rps);
TheorMRR.AvgMRRTheor=DiamRemoval*Nsteps*WP_length*pi*dw/duration/2; % mm3/sec
% %
Diametral_Infeed=DiamRemoval(1)*Nwp;
Infeed=Diametral_Infeed/(1-sin(Beta)/sin((Beta+Alpha))+sin(Alpha)/sin((Beta+Alpha)));  % mm
TheorMRR.MRRReal= Infeed*WP_rps*WP_length*pi*dw; % mm3/sec
% Return the results
DataOut.DetachCountersStage1=DetachCounterStage1;
DataOut.DetachCountersStage2=DetachCounterStage2;
end
% % % % % % % % % % End of analysis