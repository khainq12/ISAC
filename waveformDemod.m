function [X] = waveformDemod(waveform, S, N, cfg)
% X = waveformDemod(waveform, S, N, cfg)
% Common entry point to demodulate samples S with AFDM, OTFS or OFDM.
%   Inputs:
%      - waveform : 'AFDM', 'OTFS' or 'OFDM'
%      - S        : modulated sample matrix
%      - N        : number of subcarriers
%      - cfg      : struct returned by waveformMod
%   Output:
%      - X : N x T matrix of demodulated symbols

switch upper(waveform)
    case 'AFDM'
        X = AFDMdemod(S, cfg.c1, cfg.c2);
    case 'OTFS'
        Xg = OTFSdemod(S, [cfg.K, cfg.L]);
        X = reshape(Xg, N, []);
    case 'OFDM'
        X = OFDMdemod(S);
    otherwise
        error("Unknown waveform: %s", waveform);
end

end
