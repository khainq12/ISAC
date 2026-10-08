function [c1, c2] = afdmParams(N, maxDoppler, guardWidth)
% [c1, c2] = afdmParams(N, maxDoppler, guardWidth)
% Chirp frequencies satisfying AFDM integer orthogonality, as in sampleCode.m.
%   Inputs:
%      - N          : number of subcarriers
%      - maxDoppler : maximum normalised digital Doppler shift
%      - guardWidth : AFDM guard width
%   Outputs:
%      - c1, c2 : AFDM chirp frequencies

c1 = (2*(maxDoppler + guardWidth) + 1)/(2*N);
c2 = 1/(N^2*pi);

end
