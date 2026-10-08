function [fdHz, fdNorm] = velocityToDoppler(velocityKmh, carrierFreq, N, samplingFreq)
% [fdHz, fdNorm] = velocityToDoppler(velocityKmh, carrierFreq, N, samplingFreq)
% Converts a target radial velocity to Doppler shift, using the same monostatic
% (round-trip) convention as sampleCode.m: f_D = 2*v/lambda = 2*v*f_c/c.
%   Inputs:
%      - velocityKmh  : target velocity [km/h] (scalar or vector)
%      - carrierFreq  : carrier frequency [Hz]
%      - N            : channel matrix size (for normalized Doppler)
%      - samplingFreq : sampling frequency [Hz]
%   Outputs:
%      - fdHz   : physical Doppler shift [Hz]
%      - fdNorm : normalised digital Doppler shift (N*fdHz/samplingFreq)

lightspeed = physconst("lightspeed");
carrierWavelength = lightspeed/carrierFreq;

velocityMs = velocityKmh/3.6;
fdHz = 2*velocityMs/carrierWavelength;
fdNorm = fdHz*N/samplingFreq;

end
