% Ex:
% % Theor formulations for contact length estimation using High Fidelity
% for C estimation (Wp & wheel constant)
% % (23/09) MRR Estimations- Real, Normal Force, 

% (25/09) MRR across
% DiameteraL remocal, wrong gamma value!!!, Wrong duration formula!!!!

% % (19/11) beta_s as output

% % % corected the conversion micron to mm (feed velocity) 02/06/2025

% % corrected CW Rr (07/04/2025)

function [lg, lf, lcs, TheorEstimations] = contact_length(hw,gamma,Feed_vel,Wp_rps,Nrev,Wp_length,dw,ds,ds_rpm,dcw,Normal_specificenergy,C,C_CW,Rr,Eeq)
% This function calculates the geometric contact length (lg),
% the contact length due to deflection (lf), and the combined contact length (lcs).
% Inputs:
%   hw - height
%   gamma= deg
%   beta= deg
%  feed_vel micron/sec
%  wp_rps, wp_diameter, wp_length, ds (grinding wheel) diameter, wheel rpm
%   dw - Diameter of the workpiece in mm
%   ds - Diameter of the grinding wheel in mm
%  tangential specific energy j/mm3, force ratio
% % % % % % % % % % % % % % % %   E - Young's modulus (Not considered)
% % % % % % % % % % % % % % % %   Rr - roughness factotr (Not considered)
% Outputs:
%   lg - Geometric contact length in mm
%   lf - Contact length due to deflection in mm
%   lcs - Combined contact length in mm
% % % Beta Corrected (18/09)
% Estimating parametrs from the input
% Input Data Processinng
beta_s= asin(2*hw/(dw+ds)); 

TheorEstimations.beta_s= beta_s; 

betaCW=asin(2*hw/(dcw+dw));

Beta= beta_s+betaCW;

Alpha= pi/2-gamma*pi/180-(beta_s);

Diametral_Infeed= (Feed_vel/1000)/(Wp_rps); % mm

% Infeed= Radial_Infeed/(1-sin(Beta*pi/180)/sin((Beta+Alpha)*pi/180)+sin(Alpha*pi/180)/sin((Beta+Alpha)*pi/180)); % mm

Infeed=Diametral_Infeed/(1-sin(Beta)/sin((Beta+Alpha))+sin(Alpha)/sin((Beta+Alpha)));  % mm

TheorEstimations.TheorInfeed= Infeed; 

% Alpha = Alpha/180*pi; Beta = Beta/180*pi; % radians for contact estimation

% Grinding Contacts
% K1 = sin(Beta)/sin(Alpha+Beta);
% 
% K2 = sin(Alpha)/sin(Alpha+Beta);

MRR_REAL= Infeed*Wp_rps*Wp_length*pi*dw; % mm3/sec
% % % % % % % % % % % 
TheorEstimations.RealMRR= MRR_REAL; 

CuttingSpeed= ds*pi*ds_rpm/60/1000; % m/sec

Fn= MRR_REAL*(Normal_specificenergy)/CuttingSpeed; % N
% % % % % % % % % % % 
TheorEstimations.Fn= Fn;

% duration= round(Nrev/Wp_rps); 
% 
% DiamRemoval= (Feed_vel/1000)*duration; 
% 
% Theor_AvgMRR= DiamRemoval*Wp_length*pi*dw/duration/2;
% % % % % % % % % % % % % 
TheorEstimations.Theor_MRR= [];

% Calculate the geometric contact length (lg) in mm
lg = sqrt(Infeed*(1/(1/dw+1/ds)));

% Calculate the contact length due to deflection (lf) in mm
% lf = sqrt(8 * Rr^2 * Fn * (1/(1/dw+1/ds)) / (pi* E));

% % % CW and OW differentiation 
if isempty(C_CW)
    de = (1/(1/dw+1/ds));  % Equivalent Diameter (d_e)
    C= ((8*(Rr)^2)/(pi*Wp_length*Eeq));
    lf = sqrt(C * Fn * de);
else
    de = (1/(1/dw+1/dcw)); 
    C= ((8*(Rr)^2)/(pi*Wp_length*Eeq));
    lf = sqrt(C * Fn * de);
end

% Calculate the combined contact length (lcs) in mm
lcs = sqrt(lg^2 + lf^2);
end