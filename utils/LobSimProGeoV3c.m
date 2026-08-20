function varargout=LobSimProGeoV3c(dw,ds,dcw,vf,gamma,hw,DynStruct,varargin)
%.
% Test EXAMPLE:
%
% [FFTProf,GeoProf,DetachIndex,InternVar]=LobSimProGeoV2(68,500,300,0.9,0.8,0.0005,18,22,50);
%
% FFTProf=LobSimProGeo(dw,ds,dcw,KK,lc,vf,gamma,Hw,NRev);
% It launches centerless grinding dynamic simulation with a z-buffer
% workpiece representation and returns the final workpiece profile in terms
% of harmonic components (other possible output are detailed in the following).
% A process animation is also provided as option. By default, the initial
% workpiece profile is an ideal circle affected by uniform noise.
% In case the dynamics is not defined (see 'Dynamics' OPTION in the following),
% a pseudo-infinite stiffness and damping are considered.
%
% Mandatory INPUT:
%  - dw: workpiece diameter in mm
%  - ds: operating wheel diameter in mm
%  - dcw: regulating wheel diameter
%  - vf: instantaneous feed per revolute in mm/rev (1 value for each
%    angular increment.
%  - gamma: blade angle in deg
%  - Hw: workpiece heigth in mm
%  - NRev: revolutions number
%  - DynStruct: DynStruct describes the operating wheel - workpiece
%    relative dynamics, but it contains also the other parameters needed in
%    such a kind of simulation:
%
%         DynStruct.Mq [Kg] : modal mass matrix
%         DynStruct.Kq [N/mm]: modal stiffness matrix
%         DynStruct.Cq [N*s/mm]: modal dampining matrix
%         DynStruct.T: Transformation singleton matrix (row vector to pass
%              from modal DoFs to the physica relative DoF Operating wheel-
%              Workpiece
%         DynStruct.OmegaW [rad/s]: workpiece angular velocity
%         DynStruct.Ev [N/mm^2*s]: Volumetric Specific Energy/Vs*b
%         DynStruct.Ee [N/mm]: Edge specific force (Specific force for unit length * b)
%         DynStruct.Kc [N/mm]: contact stiffness
%
%         Example:
%         DynExampl.Mq=eye(3);DynExampl.Cq=eye(3);DynExampl.Kq=eye(3);
%         DynExampl.Cq=eye(3); DynExampl.T=ones(1,3);DynExampl.OmegaW=2*pi;
%         DynExampl.Ev=1e4;DynExampl.Ee=1e3; DynExampl.Kc=1e4;
%
% FFTProf=LobSimProGeo(dw,ds,dcw,KK,lc,vf,gamma,Hw,NRev,OPTION1,value1,OPTION2,value2,...);
%
% OPTIONS:
%
% ...'InitProf',filename...: it loads the initial profile (deviation error) from the file
% filename. The file contains Mitutoyo measurements. If filename is '',
% then a GUI appears to browse for the file. Measurements in [micron].
% ...'SecId',integer: it loads the nth profile contained in 'InitProf'.
% Default value=1.
%
% ...'InitProfId',filename...: it loads the ideal initial profile (deviation error) from the file
% filename. The Ideal initial profile does not come from Mitutoyo measurements,
% but it is a simple txt file containing the profile that can be loaded as
% it is. If filename is '', then a GUI appears to browse for the file. Measurements in [micron].
%
% Options 'InitProf' and 'InitProfId' are mutually excluding.
%
% ...'Display','on'/'off'...: it displays a video of the operation
%
% ...'VideoFile',filename...: it grabs the displayed video of the operation in the
% file named filename in mp4 format, adding the file name extension. If
% 'Display' option is not active, the video is void
%
% ...'Dynamics','on'/'off'...: if 'off' the sole static stiffness is
% considered (default is 'on')
%
% % ...'options',OPT...: it sets the advanced options contained in the
% structured variable OPT. Its fields can be:
%
%     .Nd: dimension of workpiece array (default 7200)
%     .NoisAmpl: amplitude of the white noise used to initialize the
%                workpiece profile if a user profile is not available in mm (default
%                0.1)
%     .CoeffAmpl: amplification of profile error for display (default 1000)')
%     .ResDisplay: resolution for workpiece display refresh in deg (default 5)
%     .FFTComputationDeadlock: if Display option is active, number of revolutions triggering the FFT
%                              of workpiece profile to be visulized (default 10)
%     .NLobes: band lobes limit of the FIR filter for contact filtering
%              (defaul 100)
%     .RlimSpan: ratio between axes limits and workpiece radius, to be
%                frozen, allowing a full display of the workpiece even when defect
%                grow (default 1.6)
%     .NewMarkParam: it sets the newmark "gamma" parameter for the
%               integration at NewMarkParam. Default value=0.5
%               it can be always set, but it make sense
%               only if Dynamics is active. Beta parameter
%               is set to NewMarkParam/2.
%     .GeoAccParam: Theoretical number of intersecting radii that define the geometrical accuracy of the z-buffer
%     .SearchAng: span angle for contact point searching [rad]
%     .SafCoeff: safety coefficient for the computation of the engagement
%                zone extent
%
%      EXAMPLE:
%
%      OPT.Nd=7200;
%      OPT.NoisAmpl=0.1e-3;
%      OPT.CoeffAmpl=1000;
%      OPT.ResDisplay=5;
%      OPT.FFTComputationDeadlock=10;
%      OPT.RlimSpan=1.6;
%      OPT.NRevIntern=10;
%      OPT.NewMarkParam=0.5;
%      OPT.GeoAccParam=3;
%      OPT.SearchAng=5*pi/180; [rad]
%      OPT.SafCoeff=3;
%      OPT.OptimalityTolerance=1e-6;
%      OPT.delta_X0=1e-3;
%
% [FFTProf,GeoProf]=LobSimProGeo(dw,ds,dcw,KK,lc,vf,gamma,Hw,NRev);
%
% GeoProf is the geometrical profile of the workpiece Measurements in mm.
%
% [FFTProf,GeoProf,DetachIndex]=(dw,ds,dcw,KK,lc,vf,gamma,Hw,NRev);
%
%  DetachIndex is a structured variable composed a follows:
%  .ang: angular vector of workpiece positions
%  .value: ratio between the detachment and not detachment simulation
%          instant
%
% [FFTProf,GeoProf,DetachIndex,InternVar]=LobSimProGeo(dw,ds,dcw,KK,lc,vf,gamma,Hw,NRev);
%
%  InternVar is a vector of structured variable sampling the evoution of the
%  simulation by providing some internal states and workpiece profile.
%  The sampling rate is governed by the advanced option "NRevIntern"
%  (default value: 20).
%       InternVar.Sam(ctIntern).Ain.Nlbs: Number of lobes vector for
%                                     intersection before clipping
%       InternVar.Sam.(ctIntern).Ain.Rev:  Sampling revolution number
%       InternVar.Sam(ctIntern).Ain.FFT:  FFT of the intersection before
%                                     clipping
%       InternVar.Sam(ctIntern).Aout.Nlbs: [The same than above, but for the
%                                      final intersection, after clipping
%       InternVar.Sam(ctIntern).Aout.Rev:
%       InternVar.Sam(ctIntern).Aout.FFT:
%       InternVar.Sam(ctIntern).RD.Nlbs:  [The same than above, but for the
%                                      profile reduction]
%       InternVar.Sam(ctIntern).RD.Rev:
%       InternVar.Sam(ctIntern).RD.FFT:
%       InternVar.Aval: History of work-wheel intersection after clipping
%       InternVar.AvalNCl: History of work-wheel intersection before clipping
%       InternVar.RDToT: History of work profile reduction
%
% M. Leonesio 04/03/2022 (v1)   GB 6/5/22: adde initial profile as optional output 4
% ML 15/07/2022 : Readrotondimetro2
% ML 26/09/2022 : No eliminazione effetto gradino dalla FFT
% ML, LP: 14/03/2023, Profile Sign Change (RDbuff)
% Variable Initialization
cToT=0;
ctRev=0;
ctDet=0;
ctIntern=0;
NLB=0;
DispFLAGBool=0;
FlagvideoGrab=0;
ctTouch=0; % Variables to take into account Touch/DeTouch occurrence
ctDetach=0; %"                                                     %
NFixedInput=7;
SecId=1;
DataOutRotInit=0; % stand off reduction to avoid any initial detachment [1 micron] 


% Default advanced options value
Nd=7200; % Size of workpiece array
NoisAmpl=1e-3; % Initial Noise on workpiece surface (if profile is not defined) [mm]
CoeffAmpl=1e3; % Deviation amplifiction factor for display
ResDisplay=5; % Degrees of rotation for each animation frame [deg]
FFTComputationDeadlock=2; % Number of revolutions for on-line FFT computation
RlimSpan=1.6; % Coefficients to resize the radius limit of the plot
NRevIntern=20; % Revolutions number lag for internal variable storage, i.e, the variable are sored each 'NRevIntern' angular steps.
NewMarkParam=0.5; % (0--> Forward Eulero; 1--> Backward Eulero)
GeoAccParam=3; % Number of engaged radii that we must intersect to guarantee the accuracy
SearchAng=5*pi/180; % [rad] Angle span for searching the real contact point
SafCoeff=50;
OptimalityTolerance=1e-6;
delta_X0=8e-3; % Stand off over-compensation [0.5 micron]

% Options Flags
CtOpt=1;
FlagInitProf=0;
FlagvideoGrab=0;
FlagInitProfId=0;
FlagDyn=1;

% Managing options
while CtOpt<nargin-NFixedInput
    switch varargin{CtOpt}
        case 'InitProf'
            DataOutRotInit=varargin{CtOpt+1};FlagInitProf=1;
        case 'SecId'
            SecId=varargin{CtOpt+1};
        case 'InitProfId'
            InitProfIdFilename=varargin{CtOpt+1};FlagInitProfId=1;
        case 'Dynamics'
            OptDyn=varargin{CtOpt+1};
            if not(strcmpi(OptDyn,'on') || strcmpi(OptDyn,'off'))
                error('Invalid value for Display parameter');
            elseif strcmpi(OptDyn,'off')
                FlagDyn=0;
            end
        case 'Display'
            DispFLAG=varargin{CtOpt+1};
            if not(strcmpi(DispFLAG,'on') || strcmpi(DispFLAG,'off'))
                error('Invalid value for Display parameter');
            elseif strcmpi(DispFLAG,'on')
                DispFLAGBool=1;
            else
                DispFLAGBool=0;
            end
        case 'VideoFile'
            VideoFileName=varargin{CtOpt+1};FlagvideoGrab=1;
        case 'options'
            if isfield(varargin{CtOpt+1},'Nd'), Nd=varargin{CtOpt+1}.Nd; end % integer
            if isfield(varargin{CtOpt+1},'NoisAmpl'), NoisAmpl=varargin{CtOpt+1}.NoisAmpl; end % Initial Noise on workpiece surface (if profile is not defined)
            if isfield(varargin{CtOpt+1},'CoeffAmpl'), CoeffAmpl=varargin{CtOpt+1}.CoeffAmpl; end % Deviation amplifiction factor for display
            if isfield(varargin{CtOpt+1},'ResDisplay'), ResDisplay=varargin{CtOpt+1}.ResDisplay; end % [deg]
            if isfield(varargin{CtOpt+1},'FFTComputationDeadlock'), FFTComputationDeadlock=varargin{CtOpt+1}.FFTComputationDeadlock; end %Number of revolutions for on-line FFT computation
            if isfield(varargin{CtOpt+1},'RlimSpan'), RlimSpan=varargin{CtOpt+1}.RlimSpan; end
            if isfield(varargin{CtOpt+1},'NRevIntern'), NRevIntern=varargin{CtOpt+1}.NRevIntern; end
            if isfield(varargin{CtOpt+1},'GeoAccParam'), GeoAccParam=varargin{CtOpt+1}.GeoAccParam; end
            if isfield(varargin{CtOpt+1},'NewMarkParam')
                NewMarkParam=varargin{CtOpt+1}.NewMarkParam;
                if NewMarkParam<0 || NewMarkParam>1
                    error('Invalid value for newmark parameter');
                end
            end
            if isfield(varargin{CtOpt+1},'SearchAng'), SearchAng=varargin{CtOpt+1}.SearchAng; end
            if isfield(varargin{CtOpt+1},'SafCoeff'), SafCoeff=varargin{CtOpt+1}.SafCoeff; end
            if isfield(varargin{CtOpt+1},'OptimalityTolerance'), OptimalityTolerance=varargin{CtOpt+1}.OptimalityTolerance; end
            if isfield(varargin{CtOpt+1},'delta_X0'), delta_X0=varargin{CtOpt+1}.delta_X0; end
        otherwise
            error(['Undefined option ', varargin{CtOpt}]);
    end
    CtOpt=CtOpt+2;
    if FlagInitProf && FlagInitProfId
        error('InitProf and InitProfId options are incompatibles');
    end
end

if not(DispFLAGBool)
    HWait=waitbar(0,'Please wait...');
end

% Definition of wrapping function for circular array
CI = @(x) (1 + mod(x-1, Nd));

% Characteristic angles and coefficient computation (definition of the
% case)

Nrev=length(vf)/Nd;

NRevIntern=min([NRevIntern, floor(Nd/2)]);
% Pre-allocation
RDToT=[]; % Storing the cumulated wear in a non-circular variable (to get WP evolution)
FnOut=[];
StOut=[];
ErrOut=[];
ExtOut=[];
AvalNCl=[];
Aval=[];
nRad=[];

betaS=asin(2*hw/(ds+dw));
betaCW=asin(2*hw/(dcw+dw));
beta=betaS+betaCW;
alpha=pi/2-gamma*pi/180-(betaS);
k1=sin(beta)/sin((alpha)+(beta));
k2=sin(alpha)/sin((alpha)+(beta));
hw

% Initial profile definition

% Check options incompatibility
if FlagInitProf && FlagInitProfId
    error('Incompatible options FlagInitProf and FlagInitProfId')
end

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

if FlagvideoGrab
    for nCurv=1:length(DataOutRotInit)
        videoFWriter{nCurv} = vision.VideoFileWriter ([VideoFileName,num2str(nCurv),'.mp4'], ...
            'FrameRate', 20, ...
            'FileFormat', 'MPEG4', ...
            'Quality', 90 ...
            );
    end
end

% If an ideal profile must be loaded
if FlagInitProfId
    if InitProfIdFilename
        DataProf = load(InitProfIdFilename);
    else
        [FILENAME, PATHNAME, ~] = uigetfile('*.txt', 'IdealProfile');
        DataProf=load([PATHNAME,FILENAME]);
    end
    [rw,cl]=size(DataProf);
    if cl>1, DataProf=DataProf'; end
    FFTInitPro=fft(DataProf);
    if mod(length(DataProf),2)~=0, error('Odd profile points not supported yet'); end
    ZRs=zeros((Nd-length(DataProf))/2,1);
    FFTInitProDense=[FFTInitPro(1:end/2); ZRs; FFTInitPro(end/2+1); ZRs; FFTInitPro(end/2+2:end)];
    RD=ifft(FFTInitProDense)*Nd/length(DataProf);
end % Vector of radius deviations [mm]
if nargout>4
    varargout{5}=RD;    % initial material removal
end

% Initial visualization of the WP
if DispFLAGBool
    AveRad=dw/2-mean(RD{nCurv}); % Computation of the average radius
    HPF=figure;
    HP=polarplot(2*pi/Nd:2*pi/Nd:2*pi,AveRad-CoeffAmpl*RD{nCurv});hold on
    %HP=polarplot(2*pi/Nd:2*pi/Nd:2*pi,AveRad/1000-RD);hold on
    polarplot([0,0], [AveRad, AveRad*1.1]);
    polarplot([-alpha,-alpha], [AveRad*0.8, AveRad*1.2]);
    polarplot([pi-beta,pi-beta], [AveRad*0.8, AveRad*1.2]);
    title(['\gamma=',num2str(gamma),' \beta=',num2str(beta*180/pi),...
        ' \alpha=',num2str(alpha*180/pi)]);
    rlim([0 dw/2*RlimSpan]);
    HT=text(3/2*pi*0.9,10,['Rev.=',num2str(0)]);
    text(3/2*pi*0.85,20,['Deviation error X',num2str(CoeffAmpl)]);
    HTL=text(3/2*pi*0.9,25,['N Lobes=',num2str(NLB)]);
    HTD=text(3/2*pi*0.9,30,'Detach Index = NA');
end

% Computation of the indeces corresponding to the delayed positions at
% rubbing wheel and work rest
ct1=round(alpha/2/pi*Nd); % Index associted to K1
ct2=round((pi-beta)/2/pi*Nd); % Index associted to K2

% Definition of matrices for time integration (see document)
DT=(2*pi)/Nd/DynStruct.OmegaW; % definition of time step
SysD=size(DynStruct.Mq); % System dimensions
A=[DynStruct.Mq DynStruct.Cq DynStruct.Kq; -DT*NewMarkParam*eye(SysD) eye(SysD) zeros(SysD);...
    zeros(SysD) -DT*NewMarkParam*eye(SysD) eye(SysD)]; % Matrix A in the for one-shot integration (see notes about Dynamic integration)
ExpMat=[DynStruct.T' zeros(SysD) zeros(SysD);...
    zeros(SysD(1),1) eye(SysD) zeros(SysD);...
    zeros(SysD(1),1) zeros(SysD) eye(SysD)]; % Auxiliary matrix to rearrange Dofs (T is a singleton row matrix: contribution of mode shapes to the contact DoF)
SolFin=[zeros(1,SysD(1)) zeros(1,SysD(1)) DynStruct.T]*inv(A)*ExpMat;
D=inv(A);
Z=[D(1,2)*DT*(1- NewMarkParam) D(1,2)+D(1,3)*DT*(1- NewMarkParam) D(1,3);...
    D(2,2)*DT*(1- NewMarkParam) D(2,2)+D(2,3)*DT*(1- NewMarkParam) D(2,3);...
    D(3,2)*DT*(1- NewMarkParam) D(3,2)+D(3,3)*DT*(1- NewMarkParam) D(3,3)];

%%% Dynamics initialization
% Initial q state
qVect=[zeros(SysD(1),1);zeros(SysD(1),1);zeros(SysD(1),1)]; % States [acc, vel, pos]

% Initial solution
VectDyn=[0;qVect(1:SysD(1))*DT*(1-NewMarkParam)+qVect(SysD(1)+1:2*SysD(1));...
    qVect(SysD(1)+1:2*SysD(1))*DT*(1- NewMarkParam)+qVect(2*SysD(1)+1:3*SysD(1))]; % Vector to compute the displacement
st=SolFin*VectDyn;
% Accuracy verification
KFact=10*max(0.5*((DynStruct.Kc*Nd)/(DynStruct.Ev*dw/2*DynStruct.OmegaW+(DynStruct.Kc*Nd))),2.1*SafCoeff*max(vf)/dw);   % Rough stiffness factor estimation (conservative as only 1 kc is considered)
NumEngRad=floor(acos((mean(dw/2-RD{1})-SafCoeff*max(vf)/KFact)/mean(dw/2-RD{1}))/(2*pi/Nd)); % Theoretical number of engaged radii
if NumEngRad<GeoAccParam % If the theoretical number of engaged radii is under a given threshold, raises an exception
    error('Not enough geometrical accuracy. Increase z-buffer dimension')
end
% Computation of the minimum tolerance
Fmin=dw*pi/Nd*DynStruct.Ee;
VectDyn=[Fmin;qVect(1:SysD(1))*DT*(1-NewMarkParam)+qVect(SysD(1)+1:2*SysD(1));...
    qVect(SysD(1)+1:2*SysD(1))*DT*(1-NewMarkParam)+qVect(2*SysD(1)+1:3*SysD(1))];
st_min=min(SolFin*VectDyn,1e-8)+1e-9; %  Computation of the displacement induced by the force (1 st attempt)
optionsFsolve=optimoptions('fsolve','FunctionTolerance',st_min,'Display','none',...
    'MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'TolX',1e-6,'OptimalityTolerance',OptimalityTolerance);  % GB:     'MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'TolX',1e-12,'OptimalityTolerance',1e-12);

% Start iterations
 MaxAngle0 = acos((dw/2-SafCoeff*vf/KFact)/(dw/2)); % Angle where to search the engagement
 EngagInd0=-floor(MaxAngle0/(2*pi/Nd)):ceil(MaxAngle0/(2*pi/Nd));
 
 for nCurv=1:length(RD),
     RDbuff=RD{nCurv};
     X0=min(RDbuff(CI(EngagInd0)))+delta_X0;
     delta_GOB=dGOB_Comp(RDbuff, dw, CI, ct1, 0, EngagInd0, Nd);
     delta_GOC=dGOC_Comp(RDbuff, dw, dcw, CI, ct2, 0, EngagInd0, Nd);
     db=-k1*delta_GOB; dr=k2*delta_GOC;     
     X=X0+db+dr;
     
     for theta=1:(Nd*Nrev) % theta is the index identifying the angular position of the workpiece
         X=X+vf(theta)/Nd; % Feed position along op wheel-WP center line
         % Quantities derived from the article of Cui et al.
         delta_GOB=dGOB_Comp(RDbuff, dw, CI, ct1, theta, SearchAng, Nd);
         delta_GOC=dGOC_Comp(RDbuff, dw, dcw, CI, ct2, theta, SearchAng, Nd);
         
         db=-k1*delta_GOB; dr=k2*delta_GOC; % Center displacement components
         % along wheel-WP contact line due to defect on  blade, control
         % wheel [mm]
         
         CentDist=ds/2+dw/2-X+db+dr; % Sommation of contribution [mm]
         %         Xscope(theta)=X;
         %         CentDistScope(theta)=CentDist-311;
         %         if theta>20*Nd, keyboard; end
         if FlagDyn
             Test = ForceComp(ds,dw,CentDist+st,DynStruct.Ev,DynStruct.Ee,DynStruct.Kc,KFact,SafCoeff,RDbuff,DT,max(vf),theta,Nd,CI); %GB   x= st
             [st,Err,exitflag]=fsolve(@(x) (SolFin*[ForceComp(ds,dw,CentDist+x,DynStruct.Ev,DynStruct.Ee,...
                 DynStruct.Kc,KFact,SafCoeff,RDbuff,DT,max(vf),theta,Nd,CI);qVect(1:SysD(1))*DT*(1-NewMarkParam)+qVect(SysD(1)+1:2*SysD(1));...
                 qVect(SysD(1)+1:2*SysD(1))*DT*(1-NewMarkParam)+qVect(2*SysD(1)+1:3*SysD(1))]-x),st,optionsFsolve);
             
             if exitflag<1 % anomalous termination?
                 fprintf('Err: %f; ExitFlag: %i; st: %f\n',Err,exitflag, st);
             end
         else
             [st,Err,exitflag]=fsolve(@(x) (ForceComp(ds,dw,CentDist+x,DynStruct.Ev,DynStruct.Ee,...
                 DynStruct.Kc,KFact,SafCoeff,RDbuff,DT,max(vf),theta,Nd,CI)/DynStruct.Kq-x),st,optionsFsolve);
             %st=0; % Only for DEBUG
         end
         [Fn,Red,RedNCl,a_actual_n, a_0_n]=ForceComp(ds,dw,CentDist+st,DynStruct.Ev,DynStruct.Ee,...
             DynStruct.Kc,KFact,SafCoeff,RDbuff,DT,max(vf),theta,Nd,CI); % It computes the grinding force
         if Fn==0 % Detachment is occurring
             ctDetach=ctDetach+1; % Increment the counter of the detachment time steps
         else
             ctTouch=ctTouch+1; % Increment the counter of the touchment time steps
         end
         % Profile update
         AvalNCl(theta)=max(RedNCl.Val); Aval(theta)=max(Red.Val); % They are the wheel-WP max radial reduction
         a_actual(theta)=a_actual_n; a_0(theta)=a_0_n;; % They are the wheel-WP infeed (actual and nominal
         
         % before and after the clipping (AvalNCl and Aval) respectively.
         
         RDbuff(CI(theta+Red.IDX))=RDbuff(CI(theta+Red.IDX))+Red.Val; % Update the reduction values
         
         % Internal states evolution
         qVect=[D(:,1) Z]*[DynStruct.T'*Fn;qVect];
         
         cToT=cToT+1; % Increasing time steps counter
         RDbuffToT{nCurv}(theta)=RDbuff(CI(theta)); % Storing the cumulated wear in a non-circular variable (to get WP evolution)
         FnOut{nCurv}(theta)=Fn;
         MRROut{nCurv}(theta)=sum(Red.Val.*(dw-RDbuff(CI(theta+Red.IDX)))*pi/Nd)/DT;
         StOut{nCurv}(theta)=st;
         ErrOut{nCurv}(theta)=Err;
         ExtOut{nCurv}(theta)=exitflag;
         nRad{nCurv}(theta)=length(Red.IDX);
         AvalOut{nCurv}=Aval;
         AvalNClOut{nCurv}=AvalNCl;
         RDOut{nCurv}=RDbuff;
         a_actualOut{nCurv}=a_actual;
         a_0Out{nCurv}=a_0;
         
         % Internal Variable saving, if requested
         if nargout>3
             if mod(theta,Nd*NRevIntern)==0 % The variable are saved after NRevIntern revolutions
                 ctIntern=ctIntern+1; % Counter increment
                 [NlbsI,FFT_IN]=fftg(1/Nd,AvalNCl(end-Nd*NRevIntern+1:end),0); % FFT of intersection before clipping
                 [NlbsO,FFT_OUT]=fftg(1/Nd,Aval(end-Nd*NRevIntern+1:end),0); % FFT of intersection after clipping
                 [NlbsR,FFT_RED]=fftg(1/Nd,RDbuffToT{nCurv}(end-Nd*NRevIntern+1:end),0); % Total profile storing
                 InternVar{nCurv}.Sam(ctIntern).Ain.Nlbs=NlbsI; % X axis of FFT, i.e. lobes number
                 InternVar{nCurv}.Sam(ctIntern).Ain.Rev=theta/Nd; % Revolutions axis of the STFT (more precisiley,
                 % Short Revolutions Fourier Transform
                 InternVar{nCurv}.Sam(ctIntern).Ain.FFT=FFT_IN; % FFT components before clipping
                 InternVar{nCurv}.Sam(ctIntern).Aout.Nlbs=NlbsO; % X axis of FFT, i.e. lobes number
                 InternVar{nCurv}.Sam(ctIntern).Aout.Rev=theta/Nd; % Revolutions axis of the STFT (more precisiley,
                 % Short Revolutions Fourier Transform
                 InternVar{nCurv}.Sam(ctIntern).Aout.FFT=FFT_OUT; % FFT components after clipping
                 InternVar{nCurv}.Sam(ctIntern).RDbuff.Nlbs=NlbsR; % X axis of FFT, i.e. lobes number
                 InternVar{nCurv}.Sam(ctIntern).RDbuff.Rev=theta/Nd; % Revolutions axis of the STFT (more precisiley,
                 % Short Revolutions Fourier Transform
                 InternVar{nCurv}.Sam(ctIntern).RDbuff.FFT=FFT_RED; % FFT components of the cumulated radius reduction clipping
             end
         end
         if mod(theta,Nd*ResDisplay/360)==0 % Figure update every Nd/ResDisplay revolution
             ctRev=ctRev+1; % it counts the number of revolutions
             if DispFLAGBool % If the display option is active
                 figure(gcf) % Resume the figure
                 AveRad=dw/2-mean(RDbuff); % Computation of the average radius
                 Dev=RDbuff-mean(RDbuff); Dev=circshift(Dev,-CI(theta)); % Computation of the average radius deviation
                 if ctRev>360/ResDisplay*FFTComputationDeadlock % Enter this code every "deadlock" time steps
                     [NLbs, STF] =fftg(1/Nd,RDbuffToT{nCurv}(theta-FFTComputationDeadlock*Nd+1:theta)',0); % it computes the FFT of the clipped intersection
                     val=max(abs(STF(NLbs>2 & NLbs<Nd/2))); NLB=NLbs(abs(STF)==val); % Identify the maximum component
                     set(HTL,'String',['N Lobes=',num2str(NLB)]); % Update N lobes text
                     set(HTD,'String',['Detach Index = ',num2str(ctDetach/ctTouch)]); % Update detach index text
                 end
                 set(HP,'RData',AveRad-Dev*CoeffAmpl); % Update figure
                 set(HT,'String',['Rev.=',num2str(round(ctRev*ResDisplay/360))]); % Update number of revolution text
                 if FlagvideoGrab % If the video must be grabbed it goes into the code
                     F=getframe; % get the frame
                     Vid(:,:,:,1)=F.cdata; % Put the frame in the video array
                     videoFWriter{nCurv}(Vid(:,:,:,1)); % Write the frame in the video
                 end
                 drawnow % Draw
             else
                 waitbar(round(ctRev*ResDisplay/360)/Nrev,HWait);
             end
         end
         if ctRev>360/ResDisplay*FFTComputationDeadlock %It stores the detach index every "deadlock" steps
             ctDet=ctDet+1;
             DetInd{nCurv}.value(ctDet)=ctDetach/ctTouch;
             DetInd{nCurv}.ang(ctDet)=theta;
         end
     end
     cToT=0;
     ctRev=0;
     ctDet=0;
     ctIntern=0;
     Vid=[];
 end
if FlagvideoGrab, release(videoFWriter{nCurv}); end % Close the video file
% % % if FlagInitProf % The number of radii is given if the workpiece is built in the function, otherwise is taken by
% % %     %  the measurement file
% % %     FFTLength=length(DataOutRotInit(1).DataUF);
% % % else
% % %     FFTLength=Nd;
% % % end
FFTLength=Nd;
% Nominal kinematic FFT (useful to get the sole waviness components
% Commented beacuse uncoherent with dynamical simulation...
% Stp=vf(theta)/(1+k2-k1); % Nominal infeed
% Tr=Stp/FFTLength:Stp/FFTLength:Stp; % Triangular shape given by the kinematic removal
% [~, FFTSimProfTr] =fftg(1/FFTLength,...
%     Tr-mean(Tr),0); % FFT of the nomunal components
%
if ~exist('DetInd') % Create detach index variable if not previously created, for
    % simulation shorter than Deadlock (very rarely...)
    DetInd=cell(length(RD));
end
for nCurv=1:length(RD),
    if nargout
        switch nargout
            case 1
                [NLbs, FFTSimProf] =fftg(1/FFTLength,...
                    RDOut{nCurv}-mean(RDOut{nCurv}),0);
                varargout{1}{nCurv}.value=FFTSimProf;%-FFTSimProfTr;
                varargout{1}{nCurv}.Nlobes=NLbs;
            case 2
                [NLbs, FFTSimProf] =fftg(1/FFTLength,...
                    RDOut{nCurv}-mean(RDOut{nCurv}),0);
                varargout{1}{nCurv}.value=FFTSimProf;%-FFTSimProfTr;
                varargout{1}{nCurv}.Nlobes=NLbs;
                varargout{2}{nCurv}=RDOut{nCurv};
            case 3
                [NLbs, FFTSimProf] =fftg(1/FFTLength,...
                    RDOut{nCurv}-mean(RDOut{nCurv}),0);
                varargout{1}{nCurv}.value=FFTSimProf;%-FFTSimProfTr;
                varargout{1}{nCurv}.Nlobes=NLbs;
                varargout{2}{nCurv}=RDOut{nCurv};
                varargout{3}{nCurv}=DetInd{nCurv};
            case {4, 5}  % for case 5, the initial profile already saved
                [NLbs, FFTSimProf] =fftg(1/FFTLength,...
                    RDOut{nCurv}-mean(RDOut{nCurv}),0);
                varargout{1}{nCurv}.value=FFTSimProf;%-FFTSimProfTr;
                varargout{1}{nCurv}.Nlobes=NLbs;
                varargout{2}{nCurv}=RDOut{nCurv};
                varargout{3}{nCurv}=DetInd{nCurv};
                InternVar{nCurv}.Aval=AvalOut{nCurv};
                InternVar{nCurv}.AvalNCl=AvalNClOut{nCurv};
                InternVar{nCurv}.RDToT=RDbuffToT{nCurv};
                InternVar{nCurv}.Fn=FnOut{nCurv};
                InternVar{nCurv}.MRR=MRROut{nCurv};
                InternVar{nCurv}.st=StOut{nCurv};
                InternVar{nCurv}.Err=ErrOut{nCurv};
                InternVar{nCurv}.ExtOut=ExtOut{nCurv};
                InternVar{nCurv}.nRad=nRad{nCurv};
                InternVar{nCurv}.a_actual=a_actualOut{nCurv};
                InternVar{nCurv}.a_0=a_0Out{nCurv};
                varargout{4}{nCurv}=InternVar;
            otherwise
                error('Error: Too many outputs!');
        end
    end
end
if DispFLAGBool
    close(HPF)
else
    close(HWait)
end

function [Force, varargout] = ForceComp(ds, dw, CentDist, Ev, Ee, Kc, KFact, SafCoeff, RD, DT, vf, theta, Nd, CI)
    % See the NOTES for explanations

    % Calculate the maximum engagement angle
    MaxAngle = acos((dw / 2 - SafCoeff * vf / KFact) / (dw / 2));
    if any(imag(MaxAngle))
        error('Imaginary index');
    end

    % Determine the number of engaged radii
    EngagInd = -floor(MaxAngle / (2 * pi / Nd)) : ceil(MaxAngle / (2 * pi / Nd));
    AlphaVect = EngagInd * 2 * pi / Nd;

    % Calculate the auxiliary angle Gamma
    GammaVect = asin((CentDist / (ds / 2)) * sin(AlphaVect));
    IDimag = imag(GammaVect) ~= 0;
    GammaVect(IDimag) = [];
    AlphaVect(IDimag) = [];
    EngagInd(IDimag) = [];

    % Calculate Beta angle
    BetaVect = GammaVect - AlphaVect;

    % Calculate the new radius
    NewRad = CentDist * sin(BetaVect) ./ sin(GammaVect);
    NewRad(isnan(NewRad)) = CentDist - (ds / 2);

    % Compute the theoretical reduction without compliance
    Reduction = (dw / 2) - RD(CI(theta + EngagInd)) - NewRad';
    RedNCl.Val = Reduction;
    RedNCl.IDX = EngagInd;

    % Find valid reductions
    IDX = Reduction >= 0;
    if all(IDX)
        error('Engagement underestimated. Decrease KFact');
    end

    % Calculate the force vector
    Dh = dw * pi / Nd;
    ForceVect = Dh * (Ev * Reduction(IDX) / DT + Ee) ./ (1 + Dh * Ev / (Kc * DT));

    % Calculate the real reduction with compliance
    ReductionReal = ((Reduction(IDX) - Dh * (Ev * Reduction(IDX) + DT * Ee) ./ (Kc * DT + Dh * Ev)) + ...
                     abs(Reduction(IDX) - Dh * (Ev * Reduction(IDX) + DT * Ee) ./ (Kc * DT + Dh * Ev))) / 2;

    % Handle no intersection case
    if isempty(Reduction(IDX))
        Red.Val = 0;
    else
        Red.Val = ReductionReal;
    end
    Red.IDX = EngagInd(IDX);
    RedNCl.Val(IDX) = ReductionReal;

    % Adjust force vector for zero reductions
    ForceVect(ReductionReal == 0) = Dh * Ee;

    % Calculate the total force
    Force = sum(ForceVect);

    % Handle optional output arguments
    if nargout > 1
        varargout{1} = Red;
    end
    if nargout > 2
        varargout{2} = RedNCl;
        varargout{3} = sum(Red.Val);
        varargout{4} = sum(Reduction(IDX));
    end
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

end                 % GB spostato qui da riga 526 per avere nested functions