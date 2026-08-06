% % % HF xo estimation

% % Xo estimation exactly equal to HF model first contact with workpiece
% % % % % (Modified: 20/05/2025)

function [X0_Estimation]= HF_X0Estimation(InitialProfile,MachProc,Nrev,n,dw,FlagInitProf,ProcessParameters)

% Temporarily suppress a specific warning
% Suppress the specific colon operand warning line 83
warning('off','all')


NormalSpecificEnergy_J_mm_3=61.881; MachProc.DynStruct.Kc= 2.4802e+04; 

MachProc.ds_rpm=1000; MachProc.wplength=146;
% % 146

delta_X0=0.5/1000;
%  % Stand off over-compensation [0.5 micron]


MachProc.cuttingspeed= ProcessParameters.ds*pi*MachProc.ds_rpm/60/1000;

MachProc.DynStruct.EnergyPar=(NormalSpecificEnergy_J_mm_3*MachProc.wplength)/(MachProc.cuttingspeed);

MachProc.DynStruct.Ev= MachProc.DynStruct.EnergyPar;

OPT.SafCoeff=40; hw= ProcessParameters.hw; dcw= ProcessParameters.dcw; ds= ProcessParameters.ds; gamma= ProcessParameters.gamma;

Nd= length(InitialProfile.DataUF);

fz=fzGen([Nrev  n], [MachProc.InfeedD/1000, 0],Nd);

DataOutRotInit=InitialProfile;

SafCoeff=OPT.SafCoeff;  vf= fz;

CI = @(x) (1 + mod(x-1, Nd));


betaS=asin(2*hw/(ds+dw));

betaCW=asin(2*hw/(dcw+dw));

beta=betaS+betaCW;

alpha=pi/2-gamma*pi/180-(betaS);

k1=sin(beta)/sin((alpha)+(beta));

k2=sin(alpha)/sin((alpha)+(beta));

% Computation of the indeces corresponding to the delayed positions at
% rubbing wheel and work rest
ct1=round(alpha/2/pi*Nd); % Index associted to K1

ct2=round((pi-beta)/2/pi*Nd); % Index associted to K2

% Reading external profile coming from measurement (ROUNDTEST Mitutoyo)
if FlagInitProf
    for nCurv=1:length(DataOutRotInit)
        FFTInitPro=fft(-DataOutRotInit(nCurv).DataUF/1000); %Calcola la FFT del pro
        if mod(length(DataOutRotInit(nCurv).DataUF/1000),2)~=0, error('Odd profile points not supported yet'); end
        LFlg=Nd-length(DataOutRotInit(nCurv).DataUF/1000);
        if LFlg>0
            ZRs=zeros((Nd-length(DataOutRotInit(nCurv).DataUF/1000))/2,1);
            FFTInitProDense=[FFTInitPro(1:end/2);ZRs;FFTInitPro(end/2+1);ZRs;FFTInitPro(end/2+2:end)];
            RD{nCurv}=ifft(FFTInitProDense)*Nd/length(DataOutRotInit(nCurv).DataUF/1000);
        end
        if LFlg==0
            RD{nCurv}=-DataOutRotInit(nCurv).DataUF/1000;
        end
        if LFlg<0
            FFTInitProDense=[FFTInitPro(1:Nd/2);FFTInitPro(end/2+1);FFTInitPro(end-Nd/2+2:end)];
            RD{nCurv}=ifft(FFTInitProDense)*Nd/length(DataOutRotInit(nCurv).DataUF/1000);
        end
    end
else % If the filename is not indicated, a random profile deviation is generated
    RD{1}=NoisAmpl*rand(Nd,1); % Vector of radius deviations [mm]   NoisAmpl*rand(Nd,1)*1000
end

RDbuff= cell2mat(RD);

% Accuracy verification
KFact=10*max(0.5*((MachProc.DynStruct.Kc*Nd)/(MachProc.DynStruct.Ev*dw/2*MachProc.DynStruct.OmegaW+(MachProc.DynStruct.Kc*Nd))),2.1*SafCoeff*max(vf)/dw);

MaxAngle0 = acos((dw/2-SafCoeff*vf/KFact)/(dw/2)); % Angle where to search the engagement

EngagInd0=-floor(MaxAngle0/(2*pi/Nd)):ceil(MaxAngle0/(2*pi/Nd));

X0=min(RDbuff(CI(EngagInd0)))+delta_X0;
delta_GOB=dGOB_Comp(RDbuff, dw, CI, ct1, 0, EngagInd0, Nd);
delta_GOC=dGOC_Comp(RDbuff, dw, dcw, CI, ct2, 0, EngagInd0, Nd);
db=-k1*delta_GOB; dr=k2*delta_GOC;
X=X0+db+dr;
X0_Estimation=X;

% % % % % % X0_Estimation=-4.623724775146132e-04;

end


function delta_GOB=dGOB_Comp(RD, dw, CI, ct1, theta, SearchAng, Nd)

% See article of Cui et al.
RangeSpan=-ceil(SearchAng/Nd):ceil(SearchAng/Nd);
AngSpan=RangeSpan*2*pi/Nd;
Range=CI(theta-ct1+RangeSpan);
delta_GOB=dw/2-max((dw/2-RD(Range)).*cos(AngSpan)');
end

function delta_GOC=dGOC_Comp(RD, dw, dcw, CI, ct2, theta, SearchAng, Nd)

% See article of Cui et al.
RangeSpan=-ceil(SearchAng/Nd):ceil(SearchAng/Nd);
AngSpan=RangeSpan*2*pi/Nd;
Range=CI(theta-ct2+RangeSpan);
delta_GOC=min((dw/2+dcw/2)-((dw/2-RD(Range)).*cos(AngSpan)')-sqrt((dcw/2)^2-((dw/2-RD(Range)).*sin(AngSpan)').^2));
end


