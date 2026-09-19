function result = metrics(estimate, reference)
%METRICS Generic error statistics for two aligned vectors.
%   RESULT = METRICS(ESTIMATE, REFERENCE) computes error = estimate -
%   reference. Both vectors must be nonempty, finite, real, and the same
%   length. Nothing gets dropped, so validSampleCount is just the full
%   vector length.
%
%   RESULT contains:
%     error            column vector, estimate - reference
%     RMSE             root-mean-square error
%     MAE              mean absolute error
%     bias             signed mean error
%     maxAbsError      maximum absolute error
%     validSampleCount number of finite aligned comparison samples
%
%   Generic on purpose, no knowledge of datasets, SOC, or any particular
%   model baked in. Just numbers in, numbers out.

arguments
    estimate (:, :) double
    reference (:, :) double
end

assert(isvector(estimate) && isvector(reference), "metrics:NonVectorInput", ...
    "Estimate and reference must both be vectors.");
assert(~isempty(estimate) && ~isempty(reference), "metrics:EmptyInput", ...
    "Estimate and reference must both be nonempty.");
assert(numel(estimate) == numel(reference), "metrics:SizeMismatch", ...
    "Estimate and reference must have equal lengths.");
assert(isreal(estimate) && isreal(reference) && all(isfinite(estimate)) && ...
    all(isfinite(reference)), "metrics:NonFiniteInput", ...
    "Estimate and reference must contain only finite real values.");

errorValues = estimate(:) - reference(:);

result = struct;
result.error = errorValues;
result.RMSE = sqrt(mean(errorValues .^ 2));
result.MAE = mean(abs(errorValues));
result.bias = mean(errorValues);
result.maxAbsError = max(abs(errorValues));
result.validSampleCount = numel(errorValues);
end
