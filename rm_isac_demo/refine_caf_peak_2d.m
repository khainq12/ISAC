function [delayRefined, dopplerRefined] = refine_caf_peak_2d(CAFabs2, rowMax, colMax, dopplerGrid, delayGrid)
%REFINE_CAF_PEAK_2D Sub-grid (continuous-valued) refinement of a CAF peak
% found by coarse discrete search, via standard 3-point parabolic
% interpolation applied independently along the delay and Doppler axes.
% This is a post-processing addition on top of the unmodified compute_caf.m
% output - it does not alter the CAF computation itself, only refines how
% the peak location is read off the already-computed surface. Same
% textbook technique as estimateDelayDopplerRefined.m used elsewhere in
% this project for the Doppler axis; here applied to both axes of a 2-D
% surface.

dopplerRefined = dopplerGrid(rowMax);
delayRefined   = delayGrid(colMax);

nRow = size(CAFabs2, 1);
nCol = size(CAFabs2, 2);

% Doppler axis (rows)
if rowMax > 1 && rowMax < nRow
    yL = CAFabs2(rowMax-1, colMax); yC = CAFabs2(rowMax, colMax); yR = CAFabs2(rowMax+1, colMax);
    denom = (yL - 2*yC + yR);
    if abs(denom) > 1e-12
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
        dopplerStep = dopplerGrid(2) - dopplerGrid(1);
        dopplerRefined = dopplerRefined + offset*dopplerStep;
    end
end

% Delay axis (columns)
if colMax > 1 && colMax < nCol
    yL = CAFabs2(rowMax, colMax-1); yC = CAFabs2(rowMax, colMax); yR = CAFabs2(rowMax, colMax+1);
    denom = (yL - 2*yC + yR);
    if abs(denom) > 1e-12
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
        delayStep = delayGrid(2) - delayGrid(1);
        delayRefined = delayRefined + offset*delayStep;
    end
end

end
