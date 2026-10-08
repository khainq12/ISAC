function [delayRefined, dopplerRefined] = refine_caf_peak_2d_log(CAFabs2, rowMax, colMax, dopplerGrid, delayGrid)
%REFINE_CAF_PEAK_2D_LOG Sub-grid refinement using the standard 3-point
% parabolic interpolation, but applied to log(|CAF|^2) instead of the
% raw power. Near a sinc^2-shaped mainlobe, log(sinc^2(x)) is quadratic
% to leading order (Taylor expansion), which is a much better local match
% for a 3-point quadratic fit than the raw (non-parabolic) power surface -
% a standard trick in spectral peak / DOA sub-bin interpolation.

dopplerRefined = dopplerGrid(rowMax);
delayRefined   = delayGrid(colMax);

nRow = size(CAFabs2, 1);
nCol = size(CAFabs2, 2);
logCAF = log(CAFabs2 + eps);

if rowMax > 1 && rowMax < nRow
    yL = logCAF(rowMax-1, colMax); yC = logCAF(rowMax, colMax); yR = logCAF(rowMax+1, colMax);
    denom = (yL - 2*yC + yR);
    if abs(denom) > 1e-12
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
        dopplerStep = dopplerGrid(2) - dopplerGrid(1);
        dopplerRefined = dopplerRefined + offset*dopplerStep;
    end
end

if colMax > 1 && colMax < nCol
    yL = logCAF(rowMax, colMax-1); yC = logCAF(rowMax, colMax); yR = logCAF(rowMax, colMax+1);
    denom = (yL - 2*yC + yR);
    if abs(denom) > 1e-12
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
        delayStep = delayGrid(2) - delayGrid(1);
        delayRefined = delayRefined + offset*delayStep;
    end
end

end
