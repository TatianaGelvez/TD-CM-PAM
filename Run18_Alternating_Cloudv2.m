%%% =====================================================================
%%% Run file to solve the passive cavitation problem in temporal signals.
%%% Methods used: FD-DAS, FD-RCB, FD-enet, FD-spTV, TD-MB-PAM
%%% =====================================================================
clear variables; close all; clc

%%% ---------------------------------------------------------------------
%%% Paths
%%% ---------------------------------------------------------------------
addpath("Support/", "Support/FD_Methods", "Support/BM4D", "Support/BM4D/bm4d", "Figures/", "Data/")

%%% ---------------------------------------------------------------------
%%% Load data
%%% ---------------------------------------------------------------------
load('Data/Alternating_Cloud_Easy2_Conv1period=9.mat')
Sgnl1 = Sgnl(:,41:end);
load('Data/Alternating_Cloud_Easy2_Conv1period=10.mat')
Sgnl = [Sgnl1  Sgnl(:,1:40)];
load('Data/Parameters_Alternating_Cloud_Easy2_Conv.mat')
load('Data/Alternating_Cloud_Easy2_Conv1period=9.mat')
%%% ---------------------------------------------------------------------
%%% Experiment setup
%%% ---------------------------------------------------------------------
par.name       = "Alternating_v2_Cloud_Easy2_Convperiod=9";   % Title of the experiment
par.doSave     = true;   % Save figures
par.doPlot     = true;   % Plot figures

%%% ---------------------------------------------------------------------
%%% State-of-the-art beamforming parameters for CSM and regularization
%%% parameters
%%% ---------------------------------------------------------------------
par.FD_K       = 130;    % Number of frequency groups
par.FD_overlap = 0.9;    % Overlap ratio
par.FD_f_recon = par.fcentrale;  % Reconstruction frequency

par.FD_lambdaL1     = 0.5;
par.FD_lambdaL2     = 0.3;
par.FD_lambdaTVlat  = 0.25;
par.FD_lambdaTVax   = 0.03;
par.FD_lambdaL1TV   = 0.3;
par.FD_rho          = 4;
par.FD_inv_opTV     = [];
par.FD_inv_op_enet  = [];

%%% ---------------------------------------------------------------------
%%% Time domain beamforming parameters
%%% ---------------------------------------------------------------------
par.thr_taus        = [7.5e-1 1.0];   % Step-size limits
par.thr_sparse      = [0.95 0.95];  % Sparsification thresholds los demas 0.9925
par.thr_smoothness  = [0.1 0.1];
par.thr_gamma       = [5e-4 2.0];     % TV regularization limits
par.thr_mu          = [4e-3 0.5];     % Denoising parameter limits
par.thr_reso        = 7e-3 / ((par.xLim(2)-par.xLim(1))/(par.NPixlsX-1));

% Parameter search ranges (can be customized)
par.RangeTau     = [1]; %%maximum 3
par.RangeLmbd    = [5]; %%maximum 5
par.RangeGamma   = [1]; %%maximum 5
par.RangeMu      = [1]; %%maximum 5
par.RangeSubImg  = [1];

%%% ---------------------------------------------------------------------
%%% Fixed experiment parameters (overriding search ranges)
%%% ---------------------------------------------------------------------
par.tau_v     = 0.95;%[0.75 ; 0.75 ; 1]; 0.75 0.95
%par.lambd_v   = [0.0278 0.14 0.150];
%par.gamma_v   = 0.001;   % Alternative value tested: 0.005 or 0.04
par.mu_v      = 5e-4;%0.003;%0.03;   % Alternative value tested: 0.005  0.01 %%% todo los demas 0.1
par.myMu      = [5e-16; 5e-9; 5e-8; 5e-7; 5e-10; 5e-10; 5e1];
par.replica   = 1;
par.NIter     = 4;%20;              % Number of iterations
par.NSampls   = 250;            % Minimum number of samples
par.Sgnl      = Sgnl;

hc.map_GT      = map_gt;
hc.MapsGt      = MapsGt;



%%% ---------------------------------------------------------------------
%%% Main execution
%%% ---------------------------------------------------------------------
Main_InversePassiveCloudsAlternating(par, hc);
