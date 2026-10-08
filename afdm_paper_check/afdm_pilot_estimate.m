function [gainsHat, delaysHat, dopplersHat] = afdm_pilot_estimate(yE, xpilot, N, P, lmax, alpha_max, xi_nu, c1, A, fractional)
%AFDM_PILOT_ESTIMATE Embedded-pilot ML delay/Doppler/gain estimator,
% following Bemani, Ksairi, Kountouris (arXiv:2204.12798), Section VI.
%
% Uses SUCCESSIVE CANCELLATION (matching pursuit): find the single
% strongest path, subtract its contribution from the residual, repeat.
% For the integer-Doppler case the paper's own eq. 73 shows the P
% candidate templates are exactly orthogonal, so successive cancellation
% and simultaneous top-P peak-picking coincide. For the fractional case
% the paper states (73) is only an ASSUMED approximation ("we exploit
% this approximation and assume that (73) holds") - diagnostic testing
% (noiseless, known ground truth) showed simultaneous top-P peak-picking
% can pick a spurious near-duplicate of the strongest path and miss a
% weaker one when templates are not perfectly orthogonal; successive
% cancellation is the standard fix for this kind of correlated-dictionary
% search and is used here as our own extension on top of the paper's
% description.
%
% yE      : received samples (Nx1) - using the full received vector for
%           simplicity; the guard-null frame design means data symbols
%           contribute ~0 here, so this is equivalent to using only the
%           formal estimation window.
% xpilot  : known pilot symbol value (scalar)
% N,P     : frame length, number of paths to estimate
% lmax    : maximum delay
% alpha_max: maximum integer Doppler
% xi_nu   : fractional Doppler guard (0 for integer case)
% c1      : AFDM chirp parameter used by the transmitter
% A       : NxN unitary AFDM analysis (DAFT) matrix (A' is used as tx)
% fractional: true to also refine fractional Doppler offset

delayCands = 0:lmax;
if fractional
    kv = alpha_max + xi_nu;
else
    kv = alpha_max;
end
alphaCands = -kv:kv;

v1 = A(1,:)'; % = A'*e1, constant across candidates

nCand = numel(delayCands)*numel(alphaCands);
templates = zeros(N, nCand);
candList = zeros(nCand,2);
idx = 1;
for li = delayCands
    for ai = alphaCands
        Hs = Gen_channel_mtx(N, 1, 1, li, ai/N, c1);
        templates(:,idx) = A * (Hs * v1); % col 1 of A*Hs*A' = A*(Hs*v1)
        candList(idx,:) = [li, ai];
        idx = idx + 1;
    end
end

aFine = -0.5:0.05:0.5;

residual = yE;
delaysHat = zeros(1,P); dopplersHat = zeros(1,P); gainsHat = zeros(1,P);

for k = 1:P
    corrVal = (templates' * residual) / xpilot;
    peakEnergy = abs(corrVal).^2 ./ max(sum(abs(templates).^2,1)', 1e-12);
    [~, best] = max(peakEnergy);
    li = candList(best,1);
    aiInt = candList(best,2);

    if fractional
        bestVal = -inf; bestA = 0; bestTmpl = templates(:,best);
        for af = aFine
            nu = aiInt + af;
            Hs = Gen_channel_mtx(N, 1, 1, li, nu/N, c1);
            tmpl = A * (Hs * v1);
            val = abs(tmpl'*residual/xpilot)^2 / max(sum(abs(tmpl).^2), 1e-12);
            if val > bestVal
                bestVal = val; bestA = af; bestTmpl = tmpl;
            end
        end
        nuHat = aiInt + bestA;
        tmplUsed = bestTmpl;
    else
        nuHat = aiInt;
        tmplUsed = templates(:,best);
    end

    hHat = (tmplUsed' * residual) / (xpilot * max(sum(abs(tmplUsed).^2), 1e-12));
    delaysHat(k) = li; dopplersHat(k) = nuHat; gainsHat(k) = hHat;

    residual = residual - xpilot * hHat * tmplUsed; % cancel and continue
end

end
