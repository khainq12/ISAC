clc; clear;
L = 32;
x = zeros(L,1); x(1) = 1;  % impulse at n=0

fprintf('--- Integer tau=5: compare circshift vs fractional matrix ---\n');
S_int = build_fractional_shift_matrix(L, 5, 6);
y_int = S_int * x;
y_ref = circshift(x, 5);
fprintf('max abs diff = %.4e (should be small)\n', max(abs(y_int - y_ref)));
[~, pk] = max(abs(y_int));
fprintf('peak of fractional-matrix output at index %d (0-based), expected 5\n', pk-1);

fprintf('\n--- Fractional tau=5.4: check peak lands between bin 5 and 6 ---\n');
S_frac = build_fractional_shift_matrix(L, 5.4, 6);
y_frac = S_frac * x;
[~, idxSorted] = sort(abs(y_frac), 'descend');
fprintf('top 3 bins (0-based) and magnitudes:\n');
for i=1:3
    fprintf('  bin %d : %.4f\n', idxSorted(i)-1, abs(y_frac(idxSorted(i))));
end

fprintf('\n--- Energy preservation check (tau=5.4) ---\n');
fprintf('sum|x|^2 = %.4f, sum|y_frac|^2 = %.4f (should be close to 1)\n', sum(abs(x).^2), sum(abs(y_frac).^2));

fprintf('=== DONE ===\n');
