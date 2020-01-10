
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


function [Opt_Nick]=SOCP_Nick_build(mpc)

% This function builds a YALMIP optimizer
% - some parameters are defined here, and cannot be modified once the optimizer is built (see section 2)
% - some parameters are defined when calling the optimizer, hence they appear as variables here (see section 3)

% Parameter "mpc" roughly follows the matpower data format, though extended with regulator and capacitor bank data


%% 1. Settings

constants;

RELAX=true; % true for augmented relaxed problem (SOCP) // false for augmented problem (NLP)

if RELAX
    fprintf('Building optimizer for augmented relaxed OPF from (Nick et al, 2017)...\n');
else
    fprintf('Building optimizer for augmented OPF from (Nick et al, 2017)...\n');
end

setFixedQ0=false; % true for adding an equality constraint on reactive power injection at root bus (enforced to value 'Qtarget')

PV_limits=APPARENT; % different options for bus withdrawal feasibility space:
					% - Q_BOUNDS for independent active and reactive power bounds
					% - POWER_FACTOR for power factor constraint
					% - APPARENT for apparent power constraint

obj_func=P_COST; % different options for objective function:
				 % - P_COST for minimisation of total cost of active power injection/withdrawal
				 % - Q_TARGET for minimisation of deviation from reactive power target
				 % - EMPTY for simply searching for a feasible solution
lPenalty=0; % penalty on current magnitudes in the objective function

amp_lim=true; % true for including limit on current magnitudes // false for omitting it
shunts=true; % true for including line shunts in the model // false for omitting them
% /!\ if shunts are omitted here, they should also be set to 0 in prepare_and_run_optimisation.m, line 62, to ensure consistent computation of the relaxed constraints slack


%% 2. Model parameters encoded within the optimizer object

n=(size(mpc.bus,1)-1); % number of non-root buses in model
root=find(mpc.bus(:,2)==3); % array index of root node
busid=mpc.bus(:,1); % list of bus ids, in case they are different than array indexes
z=mpc.branch(:,3)+sqrt(-1)*mpc.branch(:,4); % line impedances
if shunts
	b=mpc.branch(:,5); % line shunts
else
	b=zeros(n,1);
end
upstream=mpc.branch(:,1); % id of bus upstream of each line
downstream=mpc.branch(:,2); % id of bus downstream of each line
regbus=mpc.regs(:,1); % buses where voltage regulators are located
caps=mpc.caps; % matrix with capacitor bank data (first column: bus where bank is located // second column: single phase reactance [MVAr])
v0 = mpc.bus(root,end)^2; % square root bus voltage magnitude
n_gen=size(mpc.gen,1); % number of buses with generators (e.g. PV stations + slack bus)
genbus=mpc.gen(:,1); % bus ids of buses with generators
gencostlin=mpc.gencost(:,6); % linear cost factors of generators' active power production


%% 3. Model parameters defined as variables in the optimizer object

p_load = sdpvar(1,n+1);
q_load = sdpvar(1,n+1);
gen_Plo = sdpvar(1,n_gen);
gen_Pup = sdpvar(1,n_gen);
gen_Qlo = sdpvar(1,n_gen);
gen_Qup = sdpvar(1,n_gen);
vmin = sdpvar(1,1);
vmax = sdpvar(1,1);
lmax = sdpvar(1,1);
Pmax = sdpvar(1,n);
Qmax = sdpvar(1,n);
regratio = sdpvar(size(regbus,1),1);

Qtarget=sdpvar(1,1); % target for net Q absorbtion by the distribution network
Smax=sdpvar(1,1); % max apparent power at PV station (used only if PV_limits=APPARENT)
powerfactor=sdpvar(1,1); % max power factor at each PV station (used only if PV_limits=POWER_FACTOR)

%% 4. Problem Variables
% Notation is as in (Nick et al, 2017) and (Bobo et al, 2020)

s0 = sdpvar(1,1,'full','complex'); % complex power withdrawal at root bus
s = sdpvar(1,n,'full','complex'); % complex power withdrawal at each non-root bus
v = sdpvar(1,n); % square voltage magnitude at each non-root bus
f = sdpvar(1,n); % square current magnitude in central element of each line
St = sdpvar(1,n,'full','complex');
Sb = sdpvar(1,n,'full','complex');
Sthat = sdpvar(1,n,'full','complex');
Stbar = sdpvar(1,n,'full','complex');
Sbhat = sdpvar(1,n,'full','complex');
Sbbar = sdpvar(1,n,'full','complex');
vbar = sdpvar(1,n);
fbar = sdpvar(1,n);

if PV_limits==POWER_FACTOR || PV_limits==APPARENT
    p_gen=sdpvar(1,n_gen);
    q_gen=sdpvar(1,n_gen);
end


%% 5. Problem Constraints
% Constraint numbers are as in (Nick et al, 2017)

Constraints = [];

for l = 1:n % for each line
    
    up=upstream(l); % id of upstream bus
	upindex=find(busid==up); % array index of upstream bus
    dn=downstream(l); % id of downstream bus
	dnindex=find(busid==dn); % array index of downstream bus
    if isempty(find(up==regbus,1))
        k=1;
    else
        k=regratio(find(up==regbus,1))^2;
    end
    if upindex==root
        vup=k*v0;
        vbarup=k*v0;
    else
        vup=k*v(upindex);
        vbarup=k*vbar(upindex);
    end
    
    % 8a
    sm=s(dnindex);
    for  j=find(upstream==dn)'
        sm=sm+St(j);
    end
    sm=sm+z(l)*f(l)-sqrt(-1)*(vup+v(dnindex))*b(l);
    Constraints = [Constraints St(l)==sm];
    
    % 8b
    left=v(dnindex)-vup;
    right=-2*real(conj(z(l))*(St(l)+sqrt(-1)*vup*b(l)))+f(l)*abs(z(l))^2;
    Constraints = [Constraints left==right];
    
    if RELAX
        % 10
        Constraints = [Constraints rcone(St(l)+sqrt(-1)*vup*b(l),f(l)/2,vup)];
    else
        % 8c
        left=f(l)*vup;
        right=real(St(l))*real(St(l))...
            + (imag(St(l))+vup*b(l))*(imag(St(l))+vup*b(l));
        Constraints = [Constraints left==right];
    end
    
    % 8d
    sm=s(dnindex);
    for  j=find(upstream==dn)'
        sm=sm+St(j);
    end
    Constraints = [Constraints Sb(l)==sm];

    % 9g/9h
    if isempty(find(genbus==dn,1)) % no generator at this bus
        Constraints = [Constraints s(dnindex)==p_load(dnindex)+q_load(dnindex)*sqrt(-1)];
        if ~isempty(find(caps(:,1)==dn,1))
            warning('Capacitor banks on non-generator buses ignored');
        end
    else
        g=find(genbus==dn);
        
        c=find(caps(:,1)==dn,1);
        if isempty(c)
            c_Cap=0;
        else
            c_Cap=caps(find(caps(:,1)==dn,1),2);
        end        
        
        switch PV_limits
            case Q_BOUNDS
                pmin=p_load(dnindex)-gen_Plo(g);
                pmax=p_load(dnindex)-gen_Pup(g);
                qmin=q_load(dnindex)-gen_Qlo(g);
                qmax=q_load(dnindex)-gen_Qup(g);
                Constraints = [Constraints s(dnindex)<=pmax+qmax*sqrt(-1)];
                Constraints = [Constraints s(dnindex)>=pmin+qmin*sqrt(-1)];
            case POWER_FACTOR % With flexible capacitor banks
                Constraints = [Constraints gen_Plo(g)<=p_gen(g)<=gen_Pup(g)];
                Constraints = [Constraints -p_gen(g)*tan(acos(powerfactor))<=q_gen(g)<=p_gen(g)*tan(acos(powerfactor))];
                Constraints = [Constraints s(dnindex)>=(p_load(dnindex)-p_gen(g))+sqrt(-1)*(q_load(dnindex)-q_gen(g)-c_Cap)];
                Constraints = [Constraints s(dnindex)<=(p_load(dnindex)-p_gen(g))+sqrt(-1)*(q_load(dnindex)-q_gen(g))];
            case APPARENT % With fixed capacitor banks
                Constraints = [Constraints gen_Plo(g)<=p_gen(g)<=gen_Pup(g)];
                rot_matrix=[1/sqrt(2) 1/sqrt(2) 0 0;1/sqrt(2) -1/sqrt(2) 0 0;0 0 1 0; 0 0 0 1];
                rot_cone_vars=[Smax/2;1;p_gen(g);q_gen(g)];
                if RELAX
                    Constraints = [Constraints cone(rot_matrix*rot_cone_vars)];
                else
                    Constraints = [Constraints p_gen(g)*p_gen(g)+q_gen(g)*q_gen(g)<=Smax]; % Smax is in square value already
                end
                Constraints = [Constraints real(s(dnindex))==p_load(dnindex)-p_gen(g)];
                Constraints = [Constraints imag(s(dnindex))==q_load(dnindex)-q_gen(g)-c_Cap];
%                 Constraints = [Constraints imag(s(dnindex))>=q_load(dnindex)-q_gen(g)-c_Cap];
%                 Constraints = [Constraints imag(s(dnindex))<=q_load(dnindex)-q_gen(g)];
            case APPARENT_cCaps % With continuous capacitor banks
                warning('Apparent cCaps also implemented in Nick');
                Constraints = [Constraints gen_Plo(g)<=p_gen(g)<=gen_Pup(g)];
                rot_matrix=[1/sqrt(2) 1/sqrt(2) 0 0;1/sqrt(2) -1/sqrt(2) 0 0;0 0 1 0; 0 0 0 1];
                rot_cone_vars=[Smax/2;1;p_gen(g);q_gen(g)];
                if RELAX
                    Constraints = [Constraints cone(rot_matrix*rot_cone_vars)];
                else
                    Constraints = [Constraints p_gen(g)*p_gen(g)+q_gen(g)*q_gen(g)<=Smax]; % Smax is in square value already
                end
                Constraints = [Constraints real(s(dnindex))==p_load(dnindex)-p_gen(g)];
                Constraints = [Constraints imag(s(dnindex))<=q_load(dnindex)-q_gen(g)];
                Constraints = [Constraints imag(s(dnindex))>=q_load(dnindex)-q_gen(g)-c_Cap];
        end
                
    end
    
    % 11a
    sm=s(dnindex);
    for  j=find(upstream==dn)'
        sm=sm+Sthat(j);
    end
    sm=sm-sqrt(-1)*(vbarup+vbar(dnindex))*b(l);
    Constraints = [Constraints Sthat(l)==sm];
    
    % 11b
    left=vbar(dnindex)-vbarup;
    right=-2*real(conj(z(l))*(Sthat(l)+sqrt(-1)*vbarup*b(l)));
    Constraints = [Constraints left==right];
    
	% 11c
	sm=s(dnindex);
	for  j=find(upstream==dn)'
		sm=sm+Stbar(j);
	end
	sm=sm+z(l)*fbar(l)-sqrt(-1)*(vup+v(dnindex))*b(l);
	Constraints = [Constraints Stbar(l)==sm];

	% 11d
	comp1=real(Sbhat(l))+sqrt(-1)*(imag(Sbhat(l))-vbar(dnindex)*b(l));
	comp2=real(Sbbar(l))+sqrt(-1)*(imag(Sbhat(l))-vbar(dnindex)*b(l));
	comp3=real(Sbhat(l))+sqrt(-1)*(imag(Sbbar(l))-v(dnindex)*b(l));
	comp4=real(Sbbar(l))+sqrt(-1)*(imag(Sbbar(l))-v(dnindex)*b(l));
	if RELAX
		Constraints = [Constraints rcone(comp1,fbar(l)/2,v(dnindex))];
		Constraints = [Constraints rcone(comp2,fbar(l)/2,v(dnindex))];
		Constraints = [Constraints rcone(comp3,fbar(l)/2,v(dnindex))];
		Constraints = [Constraints rcone(comp4,fbar(l)/2,v(dnindex))];
	else
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp1)^2+imag(comp1)^2];
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp2)^2+imag(comp2)^2];
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp3)^2+imag(comp3)^2];
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp4)^2+imag(comp4)^2];
	end

	% 11e
	comp1=real(Sthat(l))+sqrt(-1)*(imag(Sthat(l))+vbarup*b(l));
	comp2=real(Stbar(l))+sqrt(-1)*(imag(Sthat(l))+vbarup*b(l));
	comp3=real(Sthat(l))+sqrt(-1)*(imag(Stbar(l))+vup*b(l));
	comp4=real(Stbar(l))+sqrt(-1)*(imag(Stbar(l))+vup*b(l));
	if RELAX
		Constraints = [Constraints rcone(comp1,fbar(l)/2,vup)];
		Constraints = [Constraints rcone(comp2,fbar(l)/2,vup)];
		Constraints = [Constraints rcone(comp3,fbar(l)/2,vup)];
		Constraints = [Constraints rcone(comp4,fbar(l)/2,vup)];
	else
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp1)^2+imag(comp1)^2];
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp2)^2+imag(comp2)^2];
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp3)^2+imag(comp3)^2];
		Constraints = [Constraints fbar(l)*v(dnindex)>=real(comp4)^2+imag(comp4)^2];
	end
		

	% 11f
	sm=s(dnindex);
	for  j=find(upstream==dn)'
		sm=sm+Stbar(j);
	end
	Constraints = [Constraints Sbbar(l)==sm];
    
    % 11g
    sm=s(dnindex);
    for  j=find(upstream==dn)'
        sm=sm+Sthat(j);
    end
    Constraints = [Constraints Sbhat(l)==sm];
    
    % 12c
    Constraints = [Constraints v(dnindex)>=vmin];    
    
    % 12d
    Constraints = [Constraints vbar(dnindex)<=vmax];
    
    % added: voltage constraints for downstream side of voltage regulator
    if ~isempty(find(dn==regbus,1))
        k=regratio(find(dn==regbus,1))^2;
        Constraints = [Constraints k*v(dnindex)>=vmin];
        Constraints = [Constraints k*vbar(dnindex)<=vmax];
    end
    
    if amp_lim
        % 12e - ampacity limits
        comp0=real(Sb(l))+sqrt(-1)*(imag(Sb(l)));
        comp1=real(Sbhat(l))+sqrt(-1)*(imag(Sbhat(l)));
        comp2=real(Sbbar(l))+sqrt(-1)*(imag(Sbhat(l)));
        comp3=real(Sbhat(l))+sqrt(-1)*(imag(Sbbar(l)));
        comp4=real(Sbbar(l))+sqrt(-1)*(imag(Sbbar(l)));
        if RELAX
			Constraints = [Constraints rcone(comp1,lmax/2,v(dnindex))];
			Constraints = [Constraints rcone(comp2,lmax/2,v(dnindex))];
			Constraints = [Constraints rcone(comp3,lmax/2,v(dnindex))];
			Constraints = [Constraints rcone(comp4,lmax/2,v(dnindex))];
        else
			Constraints = [Constraints lmax*v(dnindex)>=real(comp1)^2+imag(comp1)^2];
			Constraints = [Constraints lmax*v(dnindex)>=real(comp2)^2+imag(comp2)^2];
			Constraints = [Constraints lmax*v(dnindex)>=real(comp3)^2+imag(comp3)^2];
			Constraints = [Constraints lmax*v(dnindex)>=real(comp4)^2+imag(comp4)^2];
        end

        % 12f - ampacity limits
        comp0=real(St(l))+sqrt(-1)*(imag(St(l)));
        comp1=real(Sthat(l))+sqrt(-1)*(imag(Sthat(l)));
        comp2=real(Stbar(l))+sqrt(-1)*(imag(Sthat(l)));
        comp3=real(Sthat(l))+sqrt(-1)*(imag(Stbar(l)));
        comp4=real(Stbar(l))+sqrt(-1)*(imag(Stbar(l)));
        if RELAX
			Constraints = [Constraints rcone(comp1,lmax/2,vup)];
			Constraints = [Constraints rcone(comp2,lmax/2,vup)];
			Constraints = [Constraints rcone(comp3,lmax/2,vup)];
			Constraints = [Constraints rcone(comp4,lmax/2,vup)];
		else
			Constraints = [Constraints lmax*vup>=real(comp1)^2+imag(comp1)^2];
			Constraints = [Constraints lmax*vup>=real(comp2)^2+imag(comp2)^2];
			Constraints = [Constraints lmax*vup>=real(comp3)^2+imag(comp3)^2];
			Constraints = [Constraints lmax*vup>=real(comp4)^2+imag(comp4)^2];
        end 
    end
    
%     % 12g - omitted for simplicty as does not shrink the feasible space (see Nick's paper, p.5)
%     Constraints = [Constraints real(St(l))<=real(Stbar(l))];
%     Constraints = [Constraints real(Stbar(l))<=Pmax(l)];
%     Constraints = [Constraints imag(St(l))<=imag(Stbar(l))];
%     Constraints = [Constraints imag(Stbar(l))<=Qmax(l)];

if setFixedQ0
	Constraints = [Constraints imag(St(upstream==root))==Qtarget];
end

end

%%


p=sdpvar(n_gen,1);
for i=1:n_gen
    bus=find(busid==genbus(i));
    if bus==root
        p(i)=real(St(upstream==root));
    else
        p(i)=-real(s(bus))+p_load(bus);
    end
end
P_Objective=sum(gencostlin.*p);

Q_Objective=(imag(St(upstream==root))-Qtarget)^2;

switch obj_func
    case P_COST
        Objective=P_Objective;
    case Q_TARGET
        Objective=Q_Objective;
    case EMPTY
        Objective=[];
end

Objective = Objective + lPenalty*sum(f);

%%

if RELAX
    ops = sdpsettings('solver','mosek','verbose',0,'showprogress',0);
else
    ops = sdpsettings('solver','ipopt', 'showprogress',0,'verbose',0);
    ops.ipopt.max_iter=200;
end
    
if obj_func==EMPTY
    Opt_Nick=optimizer(Constraints,Objective,ops,...
        {p_load,q_load,gen_Plo,gen_Pup,gen_Qlo,gen_Qup,regratio,vmin,vmax,lmax,Pmax,Qmax,Qtarget,powerfactor,Smax},...
        {P_Objective,s0,s,St,Sb,v,f,Sthat,Stbar,Sbhat,Sbbar,vbar,fbar});
else
    Opt_Nick=optimizer(Constraints,Objective,ops,...
        {p_load,q_load,gen_Plo,gen_Pup,gen_Qlo,gen_Qup,regratio,vmin,vmax,lmax,Pmax,Qmax,Qtarget,powerfactor,Smax},...
        {Objective,s0,s,St,Sb,v,f,Sthat,Stbar,Sbhat,Sbbar,vbar,fbar});
end

fprintf('Optimizer object built\n');

end