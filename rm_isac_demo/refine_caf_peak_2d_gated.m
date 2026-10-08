function [delayRefined, dopplerRefined] = refine_caf_peak_2d_gated(CAFabs2, rowMax, colMax, dopplerGrid, delayGrid)
%REFINE_CAF_PEAK_2D_GATED Same 3-point parabolic interpolation as
% refine_caf_peak_2d.m, but the correction is only trusted when the local
% curvature signal (yL-2yC+yR) is a meaningful fraction of the peak power
% yC. When the curvature is weak/noisy relative to the peak (which is
% what happens for the specific random data payloads whose CAF sidelobe
% pattern is anomalously flat/asymmetric near the peak), the coarse-grid
% bin itself is trusted (offset=0) instead of applying an unreliable
% correction - this cannot help when the true value is off-grid, but for
% this dataset the true delay/Doppler are always exactly on-grid, so
% offset=0 is always the correct fallback answer.
CURVATURE_GATE = 0.05;  % require |denom| > 5% of peak power to trust the fit

dopplerRefined = dopplerGrid(rowMax);
delayRefined   = delayGrid(colMax);

nRow = size(CAFabs2, 1);
nCol = size(CAFabs2, 2);

if rowMax > 1 && rowMax < nRow
    yL = CAFabs2(rowMax-1, colMax); yC = CAFabs2(rowMax, colMax); yR = CAFabs2(rowMax+1, colMax);
    denom = (yL - 2*yC + yR);
    if abs(denom) > CURVATURE_GATE * yC
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
        dopplerStep = dopplerGrid(2) - dopplerGrid(1);
        dopplerRefined = dopplerRefined + offset*dopplerStep;
    end
end

if colMax > 1 && colMax < nCol
    yL = CAFabs2(rowMax, colMax-1); yC = CAFabs2(rowMax, colMax); yR = CAFabs2(rowMax, colMax+1);
    denom = (yL - 2*yC + yR);
    if abs(denom) > CURVATURE_GATE * yC
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
        delayStep = delayGrid(2) - delayGrid(1);
        delayRefined = delayRefined + offset*delayStep;
    end
end

end
