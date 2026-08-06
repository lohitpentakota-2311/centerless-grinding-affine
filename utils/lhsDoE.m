% % % Latin hyper cube design
% % % added spark out revolutions 
% % % (New approach:  26/07/2025)
function [DoE]= lhsDoE(WorkHeight,WpDiameter,RadialInfeed,Controlwheel_rpm,Bladeangle,No_revolutions,SparkOut,num_samples,BladeangleConstraint)

% Generate 100 samples
% num_samples = 100;

% Define parameter ranges
% WorkHeight = [3, 18]; % Fixed discrete values
% WpDiameter = [16, 36]; % Fixed discrete range
% RadialInfeed = [2, 8, 10, 12]; % Discrete integers
% Controlwheel_rpm = [8, 15, 20]; % Discrete integers
% Bladeangle = [30, 20]; % Discrete integers
% No_revolutions = [8, 10, 15]; % Discrete integers

% % lhs design
lhs_samples = lhsdesign(num_samples, 7);

% Map WorkHeight to fixed discrete values (3 or 18)
WorkHeight_samples = WorkHeight(ceil(lhs_samples(:, 1) * numel(WorkHeight)));

% Map WpDiameter to fixed discrete values (16 or 36)
WpDiameter_samples = WpDiameter(ceil(lhs_samples(:, 2) * numel(WpDiameter)));

% Map RadialInfeed to discrete values
RadialInfeed_samples = RadialInfeed(ceil(lhs_samples(:, 3) * numel(RadialInfeed)));

% Map Controlwheel_rpm to discrete values
Controlwheel_rpm_samples = Controlwheel_rpm(ceil(lhs_samples(:, 4) * numel(Controlwheel_rpm)));

% Map Bladeangle to discrete values
Bladeangle_samples = Bladeangle(ceil(lhs_samples(:, 5) * numel(Bladeangle)));

% Map No_revolutions to discrete values
No_revolutions_samples = No_revolutions(ceil(lhs_samples(:, 6) * numel(No_revolutions)));

SparkOut_Revolutions_samples= SparkOut(ceil(lhs_samples(:, 6) * numel(SparkOut)));

% Ensure Bladeangle has no more than 20 samples at 20
if BladeangleConstraint
else
idx_blade20 = find(Bladeangle_samples == 20);
if numel(idx_blade20) > 30
    % Change excess 20s to 30
    excess = numel(idx_blade20) - 30;
    Bladeangle_samples(idx_blade20(1:excess)) = 30;
end
end 

% Combine all samples into a single matrix
final_samples = [WorkHeight_samples(:), ...
    WpDiameter_samples(:), ...
    RadialInfeed_samples(:), ...
    Controlwheel_rpm_samples(:), ...
    Bladeangle_samples(:), ...
    No_revolutions_samples(:),SparkOut_Revolutions_samples(:)];

% Display or save the result
% disp(final_DoE);

labels = {'WorkHeight', 'WpDiameter', 'RadialInfeed', 'Controlwheel_rpm', 'Bladeangle', 'No_revolutions','SparkOut_Revolutions'};

DoE = [labels; num2cell(final_samples)];
end