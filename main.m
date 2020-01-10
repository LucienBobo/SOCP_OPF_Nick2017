
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Lucien Bobo
%
% When publishing results based on this code, please cite:
%
% L. Bobo, A. Venzke, S. Chatzivasileiadis, "Second-Order Cone Relaxations
% of the Optimal Power Flow for Active Distribution Grids", 2020. Available
% online: https://arxiv.org/abs/2001.00898
%
% M.  Nick,  R.  Cherkaoui,  J.-Y.  LeBoudec,  and  M.  Paolone,  "An  exact
% convex  formulation  of  the  optimal  power  flow  in  radial distribution
% networks including transverse components", IEEE Transactions on Automatic
% Control, 2017.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
clc;

% Choice of a grid model (roughly following the matpower data format, extended with regulator and capacitor bank data)
mpc=grid_IEEE34;
%mpc=grid_IEEE123;


% Builds an optimizer object

Opt_Nick=SOCP_Nick_build(mpc);


% Some of the key simulation parameters

maxtotpv2maxtotload=2.5; % total installed PV capacity in proportion to total peak load in the network
loadfactor=1; % load at each bus, in proportion of peak load at respective bus (0 = no load, 1 = peak load)
pvfactor=1; % PV output at each bus, in proportion of peak production capacity (0 = no PV, 1 = full PV)
voltmin=0.9; % upper bound for voltage magnitudes, in p.u.
voltmax=1.1; % lower bound for voltage magnitudes, in p.u.
ampmax=4; % upper bound for current magnitudes, in p.u.


% Prepares parameters and runs the optimisation problem
% (see inside file for obtaining solution results)

prepare_and_run_optimisation;
