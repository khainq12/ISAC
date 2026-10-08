function [rHat, fitInfo] = estimateDiversitySlope(snrdB, BER, numHighSNRPoints, floorVal)
% [rHat, fitInfo] = estimateDiversitySlope(snrdB, BER, numHighSNRPoints, floorVal)
%
% Empirical diversity order = -(slope of log10(BER) vs log10(SNR_linear)),
% estimated by ordinary least-squares linear regression over the top
% numHighSNRPoints SNR values of an existing BER/SER-vs-SNR curve.
%
% At high SNR, PEP ~ 1/SNR^r where r = rank(Phi(delta)) (see theory doc
% "LyThuyet_AFDM_OTFS_OFDM.docx", Muc 3, Eq 3.3), so this slope is a direct
% numeric estimate of r.
%
% Inputs:
%   snrdB             - 1 x K vector of tested SNR points [dB], any order
%   BER                - 1 x K vector of BER (or SER) values, same curve as snrdB
%   numHighSNRPoints  - number of highest-SNR points used in the fit window
%                       (default 3; must be >= 2)
%   floorVal          - scalar; BER values below floorVal are clamped up to
%                       floorVal before taking log10, to avoid -Inf when a
%                       finite-trial Monte Carlo run observed zero errors.
%                       No default is assumed - caller must supply it (e.g.
%                       1/(N*bps*numTrials) for BER, 1/(N*numTrials) for
%                       SER/ML-BER).
%
% Outputs:
%   rHat    - scalar, estimated diversity order (positive; = -slope)
%   fitInfo - struct with fields:
%       .snrdBUsed, .BERUsed  : the numHighSNRPoints values actually used
%       .clampedMask          : logical, which used points were clamped
%       .logSNR, .logBER      : x, y vectors used in the regression
%       .slope, .intercept    : raw polyfit(x,y,1) coefficients
%       .R2                   : coefficient of determination of the fit

if nargin < 3 || isempty(numHighSNRPoints)
    numHighSNRPoints = 3;
end
if numHighSNRPoints < 2
    error('estimateDiversitySlope:tooFewPoints', 'numHighSNRPoints must be >= 2.');
end
if nargin < 4 || isempty(floorVal)
    error('estimateDiversitySlope:noFloorVal', ...
        'floorVal must be supplied (e.g. 1/(N*bps*numTrials)); no default is assumed.');
end
if length(snrdB) ~= length(BER)
    error('estimateDiversitySlope:sizeMismatch', 'snrdB and BER must have the same length.');
end
if length(snrdB) < numHighSNRPoints
    error('estimateDiversitySlope:notEnoughData', ...
        'Only %d SNR points supplied but numHighSNRPoints=%d requested.', length(snrdB), numHighSNRPoints);
end

[snrdBSorted, order] = sort(snrdB(:).');
BERSorted = BER(:).';
BERSorted = BERSorted(order);

idx = (length(snrdBSorted) - numHighSNRPoints + 1):length(snrdBSorted);
snrdBUsed = snrdBSorted(idx);
BERUsed = BERSorted(idx);

clampedMask = BERUsed < floorVal;
if sum(clampedMask) > 1
    warning('estimateDiversitySlope:manyClamped', ...
        '%d of %d high-SNR points were clamped to floorVal (zero errors observed) - slope estimate may be biased low; consider more trials or a lower SNR window.', ...
        sum(clampedMask), numHighSNRPoints);
end
BERClamped = max(BERUsed, floorVal);

snrLinear = 10.^(snrdBUsed/10);
logSNR = log10(snrLinear);
logBER = log10(BERClamped);

p = polyfit(logSNR, logBER, 1);
slope = p(1);
intercept = p(2);
rHat = -slope;

yFit = polyval(p, logSNR);
SSres = sum((logBER - yFit).^2);
SStot = sum((logBER - mean(logBER)).^2);
if SStot == 0
    R2 = 1;
else
    R2 = 1 - SSres/SStot;
end

fitInfo = struct('snrdBUsed', snrdBUsed, 'BERUsed', BERUsed, 'clampedMask', clampedMask, ...
    'logSNR', logSNR, 'logBER', logBER, 'slope', slope, 'intercept', intercept, 'R2', R2);

end
