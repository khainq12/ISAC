function [S, cfg] = waveformMod(waveform, X, N, maxDoppler, guardWidth)
% [S, cfg] = waveformMod(waveform, X, N, maxDoppler, guardWidth)
% Common entry point to modulate symbols X with AFDM, OTFS or OFDM.
%   Inputs:
%      - waveform   : 'AFDM', 'OTFS' or 'OFDM'
%      - X          : N x T transmit symbol matrix
%      - N          : number of subcarriers
%      - maxDoppler : maximum normalised digital Doppler shift (AFDM only)
%      - guardWidth : AFDM guard width (AFDM only)
%   Outputs:
%      - S   : modulated sample matrix
%      - cfg : struct with waveform config needed by the channel/demodulator
%              (cfg.includeCCP, and cfg.c1/cfg.c2 or cfg.K/cfg.L as applicable)

switch upper(waveform)
    case 'AFDM'
        [c1, c2] = afdmParams(N, maxDoppler, guardWidth);
        S = AFDMmod(X, c1, c2);
        cfg = struct('includeCCP', true, 'c1', c1, 'c2', c2);
    case 'OTFS'
        [K, L] = otfsGridKL(N);
        Xg = reshape(X, K, L, []);
        S = OTFSmod(Xg);
        cfg = struct('includeCCP', false, 'K', K, 'L', L);
    case 'OFDM'
        S = OFDMmod(X);
        cfg = struct('includeCCP', false);
    otherwise
        error("Unknown waveform: %s", waveform);
end

end
