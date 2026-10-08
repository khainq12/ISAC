function H = build_delay_doppler_channel_frac(L, delays, dopplers, gains, Lh)
%BUILD_DELAY_DOPPLER_CHANNEL_FRAC Same structure as build_delay_doppler_channel.m
% (H = sum_p alpha_p * D(nu_p) * S(tau_p)), but S(tau) supports FRACTIONAL
% tau via a windowed-sinc circular fractional-delay operator
% (build_fractional_shift_matrix.m) instead of circshift. Reduces to the
% original for integer delays (verified in smoketest_fractional_shift.m).
if nargin < 5 || isempty(Lh), Lh = 6; end

n = (0:L-1).';
H = zeros(L, L);

for p = 1:numel(gains)
    S_tau = build_fractional_shift_matrix(L, delays(p), Lh);
    D_nu  = diag(exp(1j * 2*pi * dopplers(p) * n / L));
    H = H + gains(p) * D_nu * S_tau;
end

end
