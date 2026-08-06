
% % Secant Stiffness Estimation (03/04/2025)
% % (Chekecked the value with GB estimated value from the word file)

function SecantStiffness = calculateSecantStiffness(P, Eeq, DeltaR,Length_WP)
% This function calculates the Secant Stiffness based on the formula estimated by GB.
% Inputs:
%   P       - Applied force (N)
%   Eeq     - Equivalent Young's Modulus (Pa)
%   DeltaR  - Change in radius (m)
%   params  - Some parameter used in the formula

% Calculate Secant Stiffness
ds = (P / (pi * Eeq)) * (log(4 * pi * Eeq * DeltaR / P)-1);

% Calculate secant compliance example (SecComplEx)
SecComplEx = ds / (P * Length_WP);  % Compliance example

% Calculate Secant Stiffness (N/mm)
SecantStiffness = 1 / SecComplEx;

end