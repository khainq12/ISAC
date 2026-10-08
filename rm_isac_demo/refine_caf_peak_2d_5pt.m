function [delayRefined, dopplerRefined] = refine_caf_peak_2d_5pt(CAFabs2, rowMax, colMax, dopplerGrid, delayGrid)
%REFINE_CAF_PEAK_2D_5PT Sub-grid refinement of a CAF peak using a 5-point
% least-squares quadratic fit (vs. the 3-point exact-interpolation used
% in refine_caf_peak_2d.m). Averaging over 5 points instead of 3 reduces
% sensitivity to a single asymmetric sidelobe sample (which is what was
% found to bias the 3-point version for a fraction of random data
% payloads - the CAF sidelobe pattern near the peak depends on the
% random QPSK payload used as the correlation reference). Falls back to
% the 3-point formula near grid edges where a full 5-point window is not
% available.

dopplerRefined = dopplerGrid(rowMax);
delayRefined   = delayGrid(colMax);

nRow = size(CAFabs2, 1);
nCol = size(CAFabs2, 2);

% Doppler axis (rows)
if rowMax > 2 && rowMax < nRow - 1
    y = CAFabs2(rowMax-2:rowMax+2, colMax);
    offset = lsq5pt(y);
    dopplerStep = dopplerGrid(2) - dopplerGrid(1);
    dopplerRefined = dopplerRefined + offset*dopplerStep;
elseif rowMax > 1 && rowMax < nRow
    yL = CAFabs2(rowMax-1, colMax); yC = CAFabs2(rowMax, colMax); yR = CAFabs2(rowMax+1, colMax);
    offset = parabola3pt(yL, yC, yR);
    dopplerStep = dopplerGrid(2) - dopplerGrid(1);
    dopplerRefined = dopplerRefined + offset*dopplerStep;
end

% Delay axis (columns)
if colMax > 2 && colMax < nCol - 1
    y = CAFabs2(rowMax, colMax-2:colMax+2);
    offset = lsq5pt(y);
    delayStep = delayGrid(2) - delayGrid(1);
    delayRefined = delayRefined + offset*delayStep;
elseif colMax > 1 && colMax < nCol
    yL = CAFabs2(rowMax, colMax-1); yC = CAFabs2(rowMax, colMax); yR = CAFabs2(rowMax, colMax+1);
    offset = parabola3pt(yL, yC, yR);
    delayStep = delayGrid(2) - delayGrid(1);
    delayRefined = delayRefined + offset*delayStep;
end

end

function offset = parabola3pt(yL, yC, yR)
    denom = (yL - 2*yC + yR);
    offset = 0;
    if abs(denom) > 1e-12
        offset = 0.5*(yL - yR)/denom;
        offset = max(min(offset, 0.5), -0.5);
    end
end

function offset = lsq5pt(y)
    % y = [y(-2) y(-1) y(0) y(1) y(2)], x = [-2 -1 0 1 2]
    y = y(:).';
    ym2=y(1); ym1=y(2); y0=y(3); y1=y(4); y2=y(5);
    b = (2*(y2-ym2) + (y1-ym1)) / 10;
    a = (2*(ym2+y2) - (ym1+y1) - 2*y0) / 14;
    offset = 0;
    if abs(a) > 1e-12
        offset = -b/(2*a);
        offset = max(min(offset, 1.0), -1.0);
    end
end
