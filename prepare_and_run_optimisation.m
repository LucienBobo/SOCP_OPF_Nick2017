
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


%% Preparing parameters

% Computing the individual PV capacity at each bus so that total installed PV capacity is equal to peak load
totLoad=sum(mpc.bus(:,3));
isLoad=mpc.bus(:,3)>0;
numPVs=size(mpc.gen(2:end,:),1);
indivPV=maxtotpv2maxtotload*totLoad/numPVs;

Smax=(indivPV*1.1)^2; % square apparent power limit at each PV station (only used if PV_limits=APPARENT)
qlim2pmax=0.5; % ratio between reactive power capacity and available active power production at each PV station (only used if PV_limits=Q_BOUNDS)
powerfactor=0.9; % max power factor at each PV station (used only if PV_limits=POWER_FACTOR)

regratio = ones(size(mpc.regs,1),1); % setting voltage regulator ratios to 1 (equiv. to omitting them)
Qtarget=0; % setting a reactive power target to 0 (only used i obj_func=Q_TARGET or setFixedQ0=true)

% Loads
p_load=loadfactor*mpc.bus(:,3);
q_load=loadfactor*mpc.bus(:,4);

% Generators
n_gen=size(mpc.gen,1);
gen_Pup=[200;ones(n_gen-1,1)*pvfactor*indivPV]; % Pmax
gen_Plo=[-200;zeros(n_gen-1,1)]; % Pmin
gen_Qup=[200;ones(n_gen-1,1)*pvfactor*indivPV*qlim2pmax]; % Qmax
gen_Qlo=[-200;-ones(n_gen-1,1)*pvfactor*indivPV*qlim2pmax]; % Qmin

vmin=voltmin^2;
vmax=voltmax^2;
lmax=ampmax^2;

% Computation of Pmax and Qmax is omitted here, as constraint (12g) is omitted in SOCP_Nick_build.m (see lines 333-337)
Pmax=zeros(size(mpc.bus,1)-1,1);
Qmax=zeros(size(mpc.bus,1)-1,1);


%% Solving problem

fprintf('Computing problem solution...\n');

[solution,infeasible]=Opt_Nick({p_load',q_load',gen_Plo',gen_Pup',gen_Qlo',gen_Qup',regratio,vmin,vmax,lmax,Pmax',Qmax',Qtarget,powerfactor,Smax});
% solution cell array returns: {Objective,s0,s,St,Sb,v,f,Sthat,Stbar,Sbhat,Sbbar,vbar,fbar}

if infeasible~=0; fprintf('Problem could not be solved\n');
else; fprintf('Problem solved\n');


    %% Obtaining solution results
    
    obj=solution{1}; % objective value
    f=solution{7}; % current magnitudes
    v=solution{6}; % voltage magnitudes
    St=solution{4};
    
    root=find(mpc.bus(:,2)==3);
    s0=St(mpc.branch(:,1)==root); % injection at root bus
	rootNodeP=real(s0); % active power injection at root bus
    rootNodeQ=imag(s0); % reactive power injection at root bus
	
    s=-solution{3}; % injection at each non-root bus
    PV_prod=sum(real(s(:))+p_load(1:end-1)); % total PV output in the grid
    losses=real(s0)+sum(real(s(:))); % total losses in the grid

	% computes the largest slack on relaxed constraints
    n=size(mpc.bus,1)-1;
    slacks=zeros(n,1);
    v0=mpc.bus(root,end)^2;
    b=mpc.branch(:,5);
    for l=1:n    
        up=mpc.branch(l,1);
        if up==root
            vup=v0;
        else
            vup=v(up);
        end
        left=value(f(l))*value(vup);
        right=abs(value(St(l))+sqrt(-1)*value(vup)*value(b(l)))^2;
        slacks(l,1)=left-right;

    end
    maxSlack=max(slacks(:,1)); % largest slack on relaxed constraints
	
	exact=(maxSlack<1e-2); % estimating exactness of the problem, using a threshold of 1e-2 on largest slack on relaxed constraints

end