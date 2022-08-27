function [emission,excitation] = blueSpectra(emission,excitation)
emission = blueArrayStep(emission);
excitation = blueArrayStep(excitation);
end

function X_new = blueArrayStep(X)
X_new = X(:);
X_new = X_new(2:end);
X_new = [X_new; X_new(end)];
end