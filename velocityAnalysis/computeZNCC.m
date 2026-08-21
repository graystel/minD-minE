%computes the zero mean normalized cross-correlation
function zncc = computeZNCC(I1, I2, idx, half_win, shifts) %change to integer shiftRange?

% Extract template from I1
shiftRange = length(shifts);
template   = I1(idx - half_win : idx + half_win); % this will be compared with the second intensity trace
t_zm = normalizeSignal(template);
zncc = zeros(1, shiftRange);


% Test every possible shift within the search range
for s = 1:shiftRange
    shift = shifts(s);

    % Extract corresponding snippet from I2
    i2_start = idx - half_win + shift;
    i2_end   = idx + half_win + shift;
    snippet2 = I2(i2_start : i2_end); % what is compared to the template

    s2_zm     = normalizeSignal(snippet2); 

    % Calculate Normalized Cross-Correlation
    zncc(s) = sum(t_zm .* s2_zm) / length(template);
end
end

%normalizes a signal to have zero mean and unit variance
function normSig = normalizeSignal(sig) 

sig_mean = mean(sig);
sig_std = std(sig);

if sig_std == 0, sig_std = 1; end % Prevent division by zero
normSig = (sig - sig_mean) / sig_std;
end