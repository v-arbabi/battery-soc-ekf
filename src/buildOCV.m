function model = buildOCV(socTarget, voltageTarget)
%BUILDOCV Fit the approved OCV-SOC target with a degree-9 polynomial.
%   MODEL = BUILDOCV(SOCTARGET, VOLTAGETARGET) checks that the SOC grid
%   is ordered, unique, finite, and spans exactly [0, 1], then fits a
%   single degree-9 polynomial and gets dOCV/dSOC analytically from it
%   (polyder, not a numerical derivative). Doesn't touch or rebuild the
%   target itself, just fits it.

arguments
    socTarget (:, 1) double
    voltageTarget (:, 1) double
end

polynomialDegree = 9;

assert(numel(socTarget) == numel(voltageTarget), ...
    "buildOCV:SizeMismatch", ...
    "SOC and voltage targets must contain the same number of samples.");
assert(numel(socTarget) >= polynomialDegree + 1, ...
    "buildOCV:InsufficientSamples", ...
    "At least %d target samples are required.", polynomialDegree + 1);
assert(all(isfinite(socTarget)) && all(isfinite(voltageTarget)), ...
    "buildOCV:NonFiniteInput", ...
    "SOC and voltage targets must contain only finite values.");
assert(socTarget(1) == 0 && socTarget(end) == 1, ...
    "buildOCV:InvalidSOCRange", ...
    "The SOC target must span exactly from 0 to 1.");
assert(all(diff(socTarget) > 0), ...
    "buildOCV:InvalidSOCOrder", ...
    "The SOC target must be strictly increasing and unique.");

[polynomialCoefficients, fitStructure] = polyfit( ...
    socTarget, voltageTarget, polynomialDegree);
derivativeCoefficients = polyder(polynomialCoefficients);

fittedTargetVoltage = polyval(polynomialCoefficients, socTarget);
targetResidual = fittedTargetVoltage - voltageTarget;

model = struct;
model.degree = polynomialDegree;
model.coefficients = polynomialCoefficients;
model.derivativeCoefficients = derivativeCoefficients;
model.evaluate = @(socQuery) polyval(polynomialCoefficients, socQuery);
model.evaluateDerivative = @(socQuery) polyval(derivativeCoefficients, socQuery);
model.fitQuality = struct;
model.fitQuality.rmse_V = sqrt(mean(targetResidual .^ 2));
model.fitQuality.maxAbsError_V = max(abs(targetResidual));
model.fitQuality.endpointError_V = [targetResidual(1); targetResidual(end)];
model.fitQuality.residualNorm_V = fitStructure.normr;
model.fitQuality.degreesOfFreedom = fitStructure.df;
end
