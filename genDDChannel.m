function [H, pathDelays, pathDopplers, fadingCoefs] = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, includeCCP, equalStrength, varargin)
% [H, pathDelays, pathDopplers, fadingCoefs] = genDDChannel(channelType, N, samplingFreq, maxDelay, maxDoppler, P, includeCCP, equalStrength, [c1])
% Generates a random doubly-dispersive channel matrix realization for either an AWGN or Rayleigh-fading scenario.
% Path 1 is always the sensing target; paths 2:P (Rayleigh only, P>1) are extra paths.
%   Inputs:
%      - channelType   : 'AWGN' or 'Rayleigh'
%      - N             : channel matrix size
%      - samplingFreq  : sampling frequency [Hz]
%      - maxDelay      : maximum normalised integer delay index of any path
%      - maxDoppler    : maximum normalised digital Doppler shift of any path
%      - P             : number of paths (P == 1 required for 'AWGN')
%      - includeCCP    : include chirp-cyclic prefix effect (true/false)
%      - equalStrength : if true, all P paths are equal-strength independent Rayleigh
%                        fading (genuine multipath diversity, for BER/SER simulation);
%                        if false, paths 2:P are attenuated "clutter" around a dominant
%                        target path 1 (for sensing simulation, so the target stays
%                        identifiable amid weaker multipath)
%      - c1            : first AFDM chirp frequency (only required if includeCCP is true)
%   Outputs:
%      - H            : N x N channel matrix
%      - pathDelays   : 1 x P path delays [s], path 1 is the target
%      - pathDopplers : 1 x P path Doppler shifts [Hz], path 1 is the target
%      - fadingCoefs  : 1 x P path fading coefficients, path 1 is the target

targetDelayIdx = randi([1 maxDelay], 1, 1);
targetDopplerIdx = (2*rand(1,1) - 1)*maxDoppler;

pathDelays = targetDelayIdx/samplingFreq;
pathDopplers = targetDopplerIdx/N*samplingFreq;

switch lower(channelType)
    case 'awgn'
        if P > 1
            error("AWGN channel supports only a single (target) path.");
        end
        fadingCoefs = 1;
    case 'rayleigh'
        fadingCoefs = (randn(1,1) + 1j*randn(1,1))/sqrt(2);
        if P > 1
            if equalStrength
                clutterAtten = 1;      % genuine multipath: all paths equal-strength Rayleigh
            else
                clutterAtten = 0.15;   % clutter paths are weaker than the target on average, so the target stays identifiable
            end
            clutterDelays = randi([0 maxDelay], 1, P-1)/samplingFreq;
            clutterDopplers = (2*rand(1,P-1) - 1)*maxDoppler/N*samplingFreq;
            clutterFading = clutterAtten*(randn(1,P-1) + 1j*randn(1,P-1))/sqrt(2);
            pathDelays = [pathDelays, clutterDelays];
            pathDopplers = [pathDopplers, clutterDopplers];
            fadingCoefs = [fadingCoefs, clutterFading];
        end
    otherwise
        error("Unknown channel type: %s", channelType);
end

if includeCCP
    H = cconvLTVChannel(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, true, varargin{1});
else
    H = cconvLTVChannel(N, pathDelays, pathDopplers, fadingCoefs, samplingFreq, false);
end

end
