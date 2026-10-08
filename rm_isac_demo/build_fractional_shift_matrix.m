function S = build_fractional_shift_matrix(L, tau, Lh)
%BUILD_FRACTIONAL_SHIFT_MATRIX Circular fractional-delay operator.
% Approximates (S*x)[n] ~= x[(n-tau) mod L] for a possibly fractional tau,
% via a compact windowed-sinc kernel evaluated on the circular (wrapped)
% distance - same windowed-sinc idea as
% https://github.com/eric-hs-rou/fractional_AF_generator (AF_zerodopplercut.m),
% adapted to a circular model to stay consistent with the existing
% circshift-based (circular) delay model in build_delay_doppler_channel.m.
% For integer tau this reduces (up to the window taper) to a pure shift,
% matching circshift's behavior at the peak.
if nargin < 3 || isempty(Lh), Lh = 6; end

S = zeros(L, L);
for col = 1:L
    m = col - 1;
    for row = 1:L
        nn = row - 1;
        d = nn - m - tau;
        d = mod(d + L/2, L) - L/2;  % wrap to nearest circular representative
        if abs(d) <= Lh
            if d == 0
                sincv = 1.0;
            else
                sincv = sin(pi*d) / (pi*d);
            end
            u = abs(d) / Lh;
            w = (1 - u*u)^2;
            S(row, col) = sincv * w;
        end
    end
end
end
