function params = ecmParamsHPPC()
%ECMPARAMSHPPC Return the locked HPPC-identified 1-RC parameter set.
%   PARAMS = ECMPARAMSHPPC() returns the fixed values produced by the
%   single identification pass (src/estimateHPPC1RC.m on the 64-window
%   strict-time fit view; see docs/DECISIONS.md, "Two ECM parameter sets").
%   Hardcoded here, not recomputed - no further HPPC re-analysis,
%   refitting, or optimiser reruns, so treat these as frozen.
%
%   This is the nominal process-model parameter set for ekfEstimate.m and
%   for Experiments A, B, C and E. ecmParams.m keeps the literature-derived
%   set untouched, used only as the mismatched model for Experiment D.
%
%   Local US06 forward-voltage fit (just descriptive, not a claim this is
%   "good enough"):
%     bias        = +0.071692468 V
%     RMSE        =  0.083772801 V
%     maxAbsError =  0.575274276 V
%   (see tests/checkUS06HPPC1RCVoltage.m). Smaller than the literature set's
%   fit (+0.089247646 V bias, 0.114294767 V RMSE), but still not small.
%   This carries into the EKF results as a known process-model bias -
%   report it, don't tune it away.

params = struct;
params.R0 = 0.0312244650416;
params.R1 = 0.0113287933401;
params.tau = 5.23962539484;
params.C1 = params.tau / params.R1;
params.Q_Ah = 2.997393193;
params.Q_As = params.Q_Ah * 3600;
end
