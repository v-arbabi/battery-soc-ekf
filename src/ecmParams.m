function params = ecmParams()
%ECMPARAMS Return the provisional constant 1-RC parameter set.
%   PARAMS = ECMPARAMS() returns the fixed Panasonic NCR18650PF parameters
%   approved for the initial implementation at 25 degC. R0/R1 in ohms, C1
%   in farads, tau in seconds, Q_Ah in amp-hours, Q_As in amp-seconds.
%
%   R0, R1, C1 come from the literature at SOC = 0.5, pulled from the
%   same-cell 25 degC FFRLS/UDDS table in docs/ECM_PARAMETER_PROPOSAL.md.
%   These are constant placeholder values, no SOC dependence and no local
%   validation behind them yet.
%
%   Q_Ah is the measured discharge capacity from the local C/20 run.
%   tau and Q_As are just derived from the rest:
%     tau = R1 * C1
%     Q_As = Q_Ah * 3600

params = struct;
params.R0 = 0.004500052;
params.R1 = 0.028981830;
params.C1 = 218.8730953;
params.tau = params.R1 * params.C1;
params.Q_Ah = 2.997393193;
params.Q_As = params.Q_Ah * 3600;
end
