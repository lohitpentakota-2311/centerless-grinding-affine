% Freq to SpaceDomain Signal recontruction
% (ComplexInput:Symmetric on selected harmonics by Nof Roundness,Harmonicrange)

function [Roundness,Estimations]= freqToSpatialRoundness(Input,Nd)

% Unnormalize FFT Coeff
for i = 1:size(Input,2)
    Input(:, i) = Input(:, i) * Nd / 2;
end

Estimations.positive_half = Input(1:end,1:end);

Estimations.negative_half = conj(Input(end:-1:2,1:end));  % For signals, flipping and conj is sufficien

% adding All other harmonics to zero
Zrs= zeros(Nd-size(Estimations.negative_half,1)*2-1,size(Estimations.positive_half,2));

% Combine both halves to form the conjugate symmetric vector to Ifft (%
% Construct the conjugate symmetric vector)
conjugate_symmetric_vector = [Estimations.positive_half;Zrs;Estimations.negative_half];

% Perform the inverse FFT to space
RealGeoProf = ifft(conjugate_symmetric_vector);

imag_threshold = 1e-10; % for round off errors
% Set imaginary parts that are below the threshold to zero
RealGeoProf(imag(RealGeoProf) < imag_threshold) = real(RealGeoProf(imag(RealGeoProf) < imag_threshold));

% Verify that the result is real
if ~isreal(RealGeoProf)
    error('exceeds the imag_threshold, check for conjugate symmetric for ifft.');
end

% Profile Centering removing one per turm
NPointsPlaneMeas = length(RealGeoProf);
OneTurn = 2*pi*(1:NPointsPlaneMeas)'/NPointsPlaneMeas; % angular coordinate
CosOneTurn = cos(OneTurn); SinOneTurn = sin(OneTurn);
% Initialize
CentredRealGeoProf= zeros(NPointsPlaneMeas,size(RealGeoProf,2));
for Center= 1:size(RealGeoProf,2)
    CosOffset = CosOneTurn'*RealGeoProf(:,Center)/NPointsPlaneMeas*2; % "X dir" offset
    SinOffset = SinOneTurn'*RealGeoProf(:,Center)/NPointsPlaneMeas*2; % "Y dir" offset
    % centered data:
    CentredRealGeoProf(:,Center) = (RealGeoProf(:,Center))-SinOneTurn*SinOffset-CosOneTurn*CosOffset;
end

% Store for plot
Estimations.RealGeoProf= RealGeoProf;

% compute roundness
Roundness= zeros(1,size(CentredRealGeoProf,2));
for Sample= 1:size(CentredRealGeoProf,2)
    Roundness(Sample)= roundness(CentredRealGeoProf(:,Sample));
end

end

