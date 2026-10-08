function [Rdamaged, damageInfo] = applyChannelDamage(R, samplingFreq, damageType, damageParams)
% [Rdamaged, damageInfo] = applyChannelDamage(R, samplingFreq, damageType, damageParams)
%
% Applies a deterministic, waveform-agnostic impairment to a received
% sample-domain vector, layered on top of the usual Rayleigh
% doubly-dispersive channel + AWGN (genDDChannel.m). Defined purely in
% terms of the SAMPLE-DOMAIN signal (raw N-sample block, DFT bins / time
% indices), independent of any waveform's own subcarrier/delay-Doppler-grid
% mapping, so it is fair to apply identically across OFDM/OTFS/AFDM.
%
% Inputs:
%   R            - N x 1 received samples (already includes true channel H
%                  and noise)
%   samplingFreq - sampling frequency [Hz] (needed only for damageType
%                  'freqNotch')
%   damageType   - 'none' | 'freqNotch' | 'timeBurst'
%   damageParams - struct, fields depend on damageType:
%       'freqNotch' : .fLowHz, .fHighHz (band edges in Hz, wrapped to
%                     (-fs/2, fs/2]), .atten (linear amplitude attenuation
%                     factor, e.g. 0.05)
%       'timeBurst' : .startIdx (1-based sample index), .burstLen
%                     (contiguous sample count, must satisfy
%                     startIdx+burstLen-1 <= N), .atten
%
% Outputs:
%   Rdamaged   - N x 1 damaged received vector
%   damageInfo - struct recording which bins/samples were affected
%       'freqNotch' -> .binMask (1 x N logical), .fWrapped (1 x N, Hz)
%       'timeBurst' -> .sampleIdx (1 x burstLen)
%       'none'      -> empty struct

N = length(R);
R = R(:);

switch lower(damageType)
    case 'none'
        Rdamaged = R;
        damageInfo = struct();

    case 'freqnotch'
        f = (0:N-1)*(samplingFreq/N);
        fWrapped = f;
        fWrapped(f > samplingFreq/2) = f(f > samplingFreq/2) - samplingFreq;
        mask = (fWrapped >= damageParams.fLowHz) & (fWrapped <= damageParams.fHighHz);
        Rf = fft(R);
        Rf(mask) = Rf(mask) * damageParams.atten;
        Rdamaged = ifft(Rf);
        damageInfo = struct('binMask', mask, 'fWrapped', fWrapped);

    case 'timeburst'
        endIdx = damageParams.startIdx + damageParams.burstLen - 1;
        if endIdx > N
            error('applyChannelDamage:burstOutOfRange', ...
                'timeBurst: startIdx+burstLen-1 (%d) exceeds N (%d); no wraparound supported.', endIdx, N);
        end
        idx = damageParams.startIdx:endIdx;
        Rdamaged = R;
        Rdamaged(idx) = Rdamaged(idx) * damageParams.atten;
        damageInfo = struct('sampleIdx', idx);

    otherwise
        error('applyChannelDamage:unknownType', 'Unknown damageType: %s', damageType);
end

end
