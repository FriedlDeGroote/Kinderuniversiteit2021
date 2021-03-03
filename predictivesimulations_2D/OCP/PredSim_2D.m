%% OCP for muscle-driven 2D predictive simulations of walking
%
% Author: Antoine Falisse
% Date: 10/11/2018
%
clear all;
clc
close all;

import casadi.*

%% User inputs
num_set = [1,1,1,1,1,1,0];% 1 or 0 depending on what you want to do
% num_set = [0,1,1,1,1,1,0,];

idx_ww = 1;                 % Index row in matrix settings
subject = 'subject1';       % subject
setup.derivatives =  'AD';  % Derivative supplier

nametrial_walk.id = 'gait_14';
nametrial_walk.IK     = ['IK_',nametrial_walk.id];
time_IC = [3.5,4.2];

%% Settings
solveProblem            = num_set(1);
analyseResults          = num_set(2);
loadResults             = num_set(3);
saveSensitivityResults  = num_set(4);
reconstructGaitCycle    = num_set(5);
writeIKmotion           = num_set(6);
checkBoundsIG           = num_set(7);
% Variable settings
settings = [1.33, 4, 50, 1, 1, 1, 1, 3];

%% Loop over settings
for www = 1:length(idx_ww)
ww = idx_ww(www);

%% Settings
% Variable
v_tgt       = settings(ww,1);    % target velocity
tol_ipopt   = settings(ww,2);    % tolerance (means 1e-(tol_ipopt))
N           = settings(ww,3);    % number of mesh intervals
W.act       = settings(ww,4);    % weight muscle activations
W.back      = settings(ww,5);    % weight back torque excitations
W.GRF       = settings(ww,6);    % weight ground reaction forces
W.acc       = settings(ww,7);    % weight ground reaction forces
exp_act     = settings(ww,8);    % power muscle activations
% Fixed
W.u = 0.001;
% Filename for saving data
v_tgt_id = round(v_tgt,2);
savename = ['_case',num2str(ww),'_v',num2str(v_tgt_id*100),...
    '_T',num2str(tol_ipopt),'_N',num2str(N),'_act',num2str(W.act*10),...
    '_back',num2str(W.back*10),'_GRF',num2str(W.GRF*10),...
    '_acc',num2str(W.acc*10),'_expAct',num2str(exp_act),...
    '_',setup.derivatives];

%% External function
pathmain = pwd;
% Path to external function
if strcmp(getenv('COMPUTERNAME'),'GBW-L-W2003')
    pathexternal = ['C:\Users\u0101727\Documents\',...
        'Visual Studio 2015\Projects\external-CasADi\',...
        'PredSim_s1_2D\','PredSim_s1_2D-install\bin'];              
end
% Load external function
cd(pathexternal);
switch setup.derivatives
    case {'AD'}  
        F = external('F','PredSim_s1_2D.dll');                
    case 'FD_F'
        F = external('F','PredSim_s1_2D.dll',struct('enable_fd',true,...
            'enable_forward',false,'enable_reverse',false,...
            'enable_jacobian',false,'fd_method','forward'));
end
cd(pathmain);
% Indices in external function output
% Joint torques
jointi.pelvis.tilt  = 1; 
jointi.pelvis.tx    = 2;
jointi.pelvis.ty    = 3;
jointi.hip.l        = 4;
jointi.hip.r        = 5;
jointi.knee.l       = 6;
jointi.knee.r       = 7;
jointi.ankle.l      = 8;
jointi.ankle.r      = 9;
jointi.trunk.ext    = 10;
% Ground reaction forces
GRFi.r              = 11:12;
GRFi.l              = 13:14;
jointi.all          = jointi.pelvis.tilt:jointi.trunk.ext;
nq.all              = length(jointi.all); % #dofs
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.ty;
nq.abs              = length(ground_pelvisi); % #absolute dofs
nq.leg              = 3;
nq.trunk            = 1; 
GRFi.all            = [GRFi.r,GRFi.l];
nGRF                = length(GRFi.all); % #GRFs

%% Model info
body_mass = 62;
body_weight = body_mass*9.81;

%% Collocation scheme
[pathRepo,~,~] = fileparts(pathmain);
pathCollocationScheme = [pathRepo,'\CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[tau_root,C,D,B] = CollocationScheme(d,method);

%% MT-parameters 
muscleNames = {'hamstrings_r','bifemsh_r','glut_max_r','iliopsoas_r',...
    'rect_fem_r','vasti_r','gastroc_r','soleus_r','tib_ant_r'};  
jointNames = {'pelvis_tilt','pelvis_tx','pelvis_ty','hip_flexion_l',...
    'hip_flexion_r','knee_angle_l','knee_angle_r','ankle_angle_l',...
    'ankle_angle_r','lumbar_extension'};
pathmusclemodel = [pathRepo,'\MuscleModel'];
addpath(genpath(pathmusclemodel));    
musi = MuscleIndices_gait1018(muscleNames);
NMuscle = length(muscleNames)*2;
load([pathmusclemodel,'\MTparameters_',subject,'.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];
pathpolynomial = [pathRepo,'\Polynomials'];
addpath(genpath(pathpolynomial));
tempload = ...
    load([pathpolynomial,'\muscle_spanning_joint_INFO_',subject,'.mat']);
[Indmusi,mai] = MomentArmIndices_gait1018(muscleNames,...
    tempload.muscle_spanning_joint_INFO);

%% CasADi functions
pathCasADiFunctions = [pathRepo,'\CasADiFunctions'];
addpath(genpath(pathCasADiFunctions));
load([pathpolynomial,'\muscle_spanning_joint_INFO_',subject,'.mat']);
load([pathpolynomial,'\MuscleInfo_',subject,'.mat']);
musi_pol = musi;
NMuscle_pol = NMuscle/2;
CasADiFunctions_gait1018

%% Get experimental data
[pathRepo,~,~] = fileparts(pathmain);
pathData = [pathRepo,'\OpenSimModel\',subject];
pathVariousFunctions = [pathRepo,'\VariousFunctions'];
addpath(genpath(pathVariousFunctions));
% IK (in radian)
pathIK = [pathData,'\IK\',nametrial_walk.IK,'.mot'];
Qs = getIK(pathIK,jointNames);

%% Bounds
pathBounds = [pathRepo,'\Bounds'];
addpath(genpath(pathBounds));
[bounds,scaling] = getBounds_gait1018(NMuscle,nq,jointi);

%% Initial guess
pathIG = [pathRepo,'\Guess'];
addpath(genpath(pathIG));
% guess = getGuess_gait1018(N,nq,NMuscle,body_weight,scaling,v_tgt,jointi,nGRF);
guess = getGuess_gait1018_IK(Qs,N,nq,NMuscle,body_weight,scaling,time_IC,jointi,nGRF);

%% Observe bounds and initial guess
if checkBoundsIG
    pathPlots = [pathRepo,'\Plots'];
    addpath(genpath(pathPlots));
    plotBoundsIG_s1
end

%% Formulate the NLP
if solveProblem
    % Empty NLP
    w   = {};
    w0  = [];
    lbw = [];
    ubw = [];
    J   = 0;
    g   = {};
    lbg = [];
    ubg = [];
    % Static parameters
    % Final time
    tf              = MX.sym('tf',1);
    w               = [w {tf}];
    lbw             = [lbw; 0.1];
    ubw             = [ubw; 1];
    w0              = [w0;  0.6];
    % States
    % Muscle activations
    a0              = MX.sym('a0',NMuscle);
    w               = [w {a0}];
    lbw             = [lbw; bounds.a.lower'];
    ubw             = [ubw; bounds.a.upper'];
    w0              = [w0;  guess.a(1,:)'];
    % Tendon forces
    FTtilde0        = MX.sym('FTtilde0',NMuscle);
    w               = [w {FTtilde0}];
    lbw             = [lbw; bounds.FTtilde.lower'];
    ubw             = [ubw; bounds.FTtilde.upper'];
    w0              = [w0;  guess.FTtilde(1,:)'];
    % Joint positions and velocities
    X0              = MX.sym('X0',2*nq.all);
    w               = [w {X0}];    
    lbw             = [lbw; bounds.QsQdots_0.lower'];
    ubw             = [ubw; bounds.QsQdots_0.upper'];    
    w0              = [w0;  guess.QsQdots(1,:)'];
    % Back activations
    a_b0            = MX.sym('a_b0',nq.trunk);
    w               = [w {a_b0}];
    lbw             = [lbw; bounds.a_b.lower'];
    ubw             = [ubw; bounds.a_b.upper'];
    w0              = [w0;  guess.a_b(1,:)'];
       
    % Pre-allocation states    
    for k=0:N
        Xk{k+1,1} = MX.sym(['X_' num2str(k+1)], 2*nq.all);
    end 
    
    % Initial point
    ak          = a0;
    FTtildek    = FTtilde0;
    Xk{1,1}     = X0;
    a_bk        = a_b0; 
    
    % Initial position pelvis_tx (should be 0)
    pelvis_tx0 = Xk{1,1}(2*jointi.pelvis.tx-1,1).*...
        scaling.QsQdots(2*jointi.pelvis.tx-1);  
    % Final position pelvis_tx
    pelvis_txf = Xk{N+1,1}(2*jointi.pelvis.tx-1,1).*...
        scaling.QsQdots(2*jointi.pelvis.tx-1);
    % Distance traveled
    dist_trav_tot = pelvis_txf-pelvis_tx0;
    
    h = tf/N;
    % loop over mesh points
    for k=0:N-1
        % Controls at mesh points (piecewise-constant in mesh intervals): 
        % Time derivative of muscle activations (states)
        vAk                 = MX.sym(['vA_' num2str(k)], NMuscle);
        w                   = [w {vAk}];
        lbw                 = [lbw; bounds.vA.lower'];
        ubw                 = [ubw; bounds.vA.upper'];
        w0                  = [w0; guess.vA(k+1,:)'];
        % Time derivative of muscle-tendon forces (states)
        dFTtildek           = MX.sym(['dFTtilde_' num2str(k)], NMuscle);
        w                   = [w {dFTtildek}];
        lbw                 = [lbw; bounds.dFTtilde.lower'];
        ubw                 = [ubw; bounds.dFTtilde.upper'];
        w0                  = [w0; guess.dFTtilde(k+1,:)'];  
        % Time derivative of joint velocities (states) 
        Ak                  = MX.sym(['A_' num2str(k)], nq.all);
        w                   = [w {Ak}];
        lbw                 = [lbw; bounds.Qdotdots.lower'];
        ubw                 = [ubw; bounds.Qdotdots.upper'];
        w0                  = [w0; guess.Qdotdots(k+1,:)'];    
        % Ground reaction forces
        GRFck               = MX.sym(['GRFc_' num2str(k)], nGRF);
        w                   = [w {GRFck}];
        lbw                 = [lbw; bounds.GRF.lower'];
        ubw                 = [ubw; bounds.GRF.upper'];
        w0                  = [w0; guess.GRF(k+1,:)'];
        % Back excitations
        e_bk                = MX.sym(['e_b_' num2str(k)], nq.trunk);
        w                   = [w {e_bk}];
        lbw                 = [lbw; bounds.e_b.lower'];
        ubw                 = [ubw; bounds.e_b.upper'];
        w0                  = [w0; guess.e_b(k+1,:)'];
      
        % States at collocation points:     
        % Muscle activations
        akj = {};
        for j=1:d
            akj{j}  = MX.sym(['	a_' num2str(k) '_' num2str(j)], NMuscle);
            w       = {w{:}, akj{j}};
            lbw     = [lbw; bounds.a.lower'];
            ubw     = [ubw; bounds.a.upper'];
            w0      = [w0;  guess.a(k+1,:)'];
        end   
        % Tendon forces
        FTtildekj = {};
        for j=1:d
            FTtildekj{j} = ...
                MX.sym(['FTtilde_' num2str(k) '_' num2str(j)], NMuscle);
            w            = {w{:}, FTtildekj{j}};
            lbw          = [lbw; bounds.FTtilde.lower'];
            ubw          = [ubw; bounds.FTtilde.upper'];
            w0           = [w0;  guess.FTtilde(k+1,:)'];
        end
        % Joint positions and velocities        
        Xkj = {};
        for j=1:d
            Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], 2*nq.all);
            w      = {w{:}, Xkj{j}};
            lbw    = [lbw; bounds.QsQdots.lower'];
            ubw    = [ubw; bounds.QsQdots.upper'];
            w0     = [w0;  guess.QsQdots(k+1,:)'];
        end   
        % Back activations
        a_bkj = {};
        for j=1:d
            a_bkj{j}= MX.sym(['	a_b_' num2str(k) '_' num2str(j)], nq.trunk);
            w       = {w{:}, a_bkj{j}};
            lbw     = [lbw; bounds.a_b.lower'];
            ubw     = [ubw; bounds.a_b.upper'];
            w0      = [w0;  guess.a_b(k+1,:)'];
        end   
        
        % Unscale values
        Xk_nsc          = Xk{k+1,1}.*scaling.QsQdots';
        FTtildek_nsc    = FTtildek.*(scaling.FTtilde');
        Ak_nsc          = Ak.*scaling.Qdotdots';
        for j=1:d
            Xkj_nsc{j} = Xkj{j}.*scaling.QsQdots';
            FTtildekj_nsc{j} = FTtildekj{j}.*scaling.FTtilde';
        end   
        
        % Call external function
        [Tk] = F([Xk_nsc;Ak_nsc]);   
        
        % Get muscle-tendon lengths and lengthening velocities
        % Left leg
        qin_l       = [Xk_nsc(jointi.hip.l*2-1,1),...
            Xk_nsc(jointi.knee.l*2-1,1),Xk_nsc(jointi.ankle.l*2-1,1)];  
        qdotin_l    = [Xk_nsc(jointi.hip.l*2,1),...
            Xk_nsc(jointi.knee.l*2,1),Xk_nsc(jointi.ankle.l*2,1)];  
        [lMTk_l,vMTk_l,MA_l] = f_lMT_vMT_dM(qin_l,qdotin_l);    
        MA_hip_l    =  MA_l(mai(1).mus.l',1);
        MA_knee_l   =  MA_l(mai(2).mus.l',2);
        MA_ankle_l  =  MA_l(mai(3).mus.l',3);    
        % Right leg
        qin_r      = [Xk_nsc(jointi.hip.r*2-1,1),...
            Xk_nsc(jointi.knee.r*2-1,1),Xk_nsc(jointi.ankle.r*2-1,1)];  
        qdotin_r    = [Xk_nsc(jointi.hip.r*2,1),...
            Xk_nsc(jointi.knee.r*2,1), Xk_nsc(jointi.ankle.r*2,1)];      
        [lMTk_r,vMTk_r,MA_r] = f_lMT_vMT_dM(qin_r,qdotin_r);  
        % Indices from left since vector 1:NMuscle
        MA_hip_r    = MA_r(mai(1).mus.l',1);
        MA_knee_r   = MA_r(mai(2).mus.l',2);
        MA_ankle_r  = MA_r(mai(3).mus.l',3);
        % Combine both legs
        lMTk_lr     = [lMTk_l;lMTk_r];
        vMTk_lr     = [vMTk_l;vMTk_r];   
        % Get muscle-tendon forces and derive Hill-equilibrium
        [Hilldiffk,FTk] = f_ForceEquilibrium_FtildeState(...
                ak,FTtildek.*scaling.FTtilde',...
                dFTtildek.*scaling.dFTtilde,lMTk_lr,vMTk_lr);            
        % Loop over collocation points
        Xk_nsc_end          = D(1)*Xk_nsc;
        FTtildek_nsc_end    = D(1)*FTtildek_nsc;
        ak_end              = D(1)*ak;
        a_bk_end            = D(1)*a_bk;
        for j=1:d
            % Expression for the state derivatives at the collocation point
            xp_nsc          = C(1,j+1)*Xk_nsc;
            FTtildep_nsc    = C(1,j+1)*FTtildek_nsc;
            ap              = C(1,j+1)*ak;
            a_bp            = C(1,j+1)*a_bk;
            for r=1:d
                xp_nsc       = xp_nsc + C(r+1,j+1)*Xkj_nsc{r};
                FTtildep_nsc = FTtildep_nsc + C(r+1,j+1)*FTtildekj_nsc{r};
                ap           = ap + C(r+1,j+1)*akj{r};
                a_bp         = a_bp + C(r+1,j+1)*a_bkj{r};
            end 
            % Append collocation equations (implicit formulation)
            % Defect/dynamic constraints are scaled using the same scale
            % factors as was used to scale the states (Patterson&Rao2013).
            % See also Betts2010, section 4.8 Scaling (p166-168)
            % Activation dynamics (implicit formulation)  
            g       = {g{:}, (h*vAk.*scaling.vA - ap)./scaling.a};
            lbg     = [lbg; zeros(NMuscle,1)];
            ubg     = [ubg; zeros(NMuscle,1)]; 
            % Contraction dynamics (implicit formulation)            
            g       = {g{:}, (h*dFTtildek.*scaling.dFTtilde - ...
                FTtildep_nsc)./(scaling.FTtilde')};
            lbg     = [lbg; zeros(NMuscle,1)];
            ubg     = [ubg; zeros(NMuscle,1)];
            % Skeletal dynamics (implicit formulation)   
            xj_nsc  = [...
                Xkj_nsc{j}(2); Ak_nsc(1); Xkj_nsc{j}(4); Ak_nsc(2);...
                Xkj_nsc{j}(6); Ak_nsc(3); Xkj_nsc{j}(8); Ak_nsc(4);...
                Xkj_nsc{j}(10); Ak_nsc(5); Xkj_nsc{j}(12); Ak_nsc(6);...
                Xkj_nsc{j}(14); Ak_nsc(7); Xkj_nsc{j}(16); Ak_nsc(8);...
                Xkj_nsc{j}(18); Ak_nsc(9); Xkj_nsc{j}(20); Ak_nsc(10)];
            g       = {g{:}, (h*xj_nsc - xp_nsc)./(scaling.QsQdots')};
            lbg     = [lbg; zeros(2*nq.all,1)];
            ubg     = [ubg; zeros(2*nq.all,1)];   
            % Back activation dynamics (explicit formulation)  
            dadt    = f_BackActivationDynamics(e_bk,a_bkj{j});
            g       = {g{:}, (h*dadt - a_bp)./scaling.a_b};
            lbg     = [lbg; zeros(nq.trunk,1)];
            ubg     = [ubg; zeros(nq.trunk,1)]; 
            % Add contribution to the end state
            Xk_nsc_end = Xk_nsc_end + D(j+1)*Xkj_nsc{j};
            FTtildek_nsc_end = FTtildek_nsc_end + D(j+1)*FTtildekj_nsc{j};
            ak_end = ak_end + D(j+1)*akj{j};  
            a_bk_end = a_bk_end + D(j+1)*a_bkj{j};    
            % Add contribution to quadrature function
            J = J + 1/(dist_trav_tot)*(...
                W.act*B(j+1)    *(f_sunsqr_exp(akj{j},exp_act))*h + ...
                W.back*B(j+1)   *(sumsqr(e_bk))*h +... 
                W.GRF*B(j+1)    *(sumsqr(GRFck))*h + ...   
                W.acc*B(j+1)    *(sumsqr(Ak))*h + ...                          
                W.u*B(j+1)      *(sumsqr(vAk))*h + ...
                W.u*B(j+1)      *(sumsqr(dFTtildek))*h);             
        end            
                  
        % Add path constraints
        % Residuals of absolute dofs are null
        g           = {g{:},Tk(ground_pelvisi,1)};
        lbg         = [lbg; zeros(nq.abs,1)];
        ubg         = [ubg; zeros(nq.abs,1)];    
        % Hip constraint (left), joint actuated by muscles
        Ft_hip_l    = FTk(mai(1).mus.l',1);
        T_hip_l     = f_T4(MA_hip_l,Ft_hip_l);
        g           = {g{:},(Tk(jointi.hip.l,1)-(T_hip_l))};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];    
        % Hip constraint (right), joint actuated by muscles
        Ft_hip_r    = FTk(mai(1).mus.r',1);
        T_hip_r     = f_T4(MA_hip_r,Ft_hip_r);
        g           = {g{:},(Tk(jointi.hip.r,1)-(T_hip_r))};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];    
        % Knee constraint (left), joint actuated by muscles 
        Ft_knee_l   = FTk(mai(2).mus.l',1);
        T_knee_l    = f_T5(MA_knee_l,Ft_knee_l);
        g           = {g{:},(Tk(jointi.knee.l,1)-(T_knee_l))};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];    
        % Knee constraint (right), joint actuated by muscles
        Ft_knee_r   = FTk(mai(2).mus.r',1);
        T_knee_r    = f_T5(MA_knee_r,Ft_knee_r);
        g           = {g{:},(Tk(jointi.knee.r,1)-(T_knee_r))};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];    
        % Ankle constraint (left), joint actuated by muscles
        Ft_ankle_l  = FTk(mai(3).mus.l',1);
        T_ankle_l   = f_T3(MA_ankle_l,Ft_ankle_l);
        g           = {g{:},(Tk(jointi.ankle.l,1)-(T_ankle_l))};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];    
        % Ankle constraint (right), joint actuated by muscles
        Ft_ankle_r  = FTk(mai(3).mus.r',1);
        T_ankle_r   = f_T3(MA_ankle_r,Ft_ankle_r);
        g           = {g{:},(Tk(jointi.ankle.r,1)-(T_ankle_r))};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];      
        % Trunk constraint, joint actuated by torque actuator
        g           = {g{:},Tk(jointi.trunk.ext,1)./scaling.BackTau-a_bk};
        lbg         = [lbg; 0];
        ubg         = [ubg; 0];
        % Activation dynamics constraint
        tact = 0.015;
        tdeact = 0.06;
        act1 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tdeact);
        act2 = vAk*scaling.vA + ak./(ones(size(ak,1),1)*tact);
        % act1
        g               = {g{:},act1};
        lbg             = [lbg; zeros(NMuscle,1)];
        ubg             = [ubg; inf*ones(NMuscle,1)]; 
        % act2
        g               = {g{:},act2};
        lbg             = [lbg; -inf*ones(NMuscle,1)];
        ubg             = [ubg; ones(NMuscle,1)./(ones(NMuscle,1)*tact)];        
        % hill-equilibrium constraint
        g               = {g{:},Hilldiffk};
        lbg             = [lbg; zeros(NMuscle,1)];
        ubg             = [ubg; zeros(NMuscle,1)];  
        % GRF constraint
        g               = {g{:},GRFck - Tk(GRFi.all,1)./scaling.GRF'};
        lbg             = [lbg; zeros(nGRF,1)];
        ubg             = [ubg; zeros(nGRF,1)];         
        % New NLP variables for states at end of interval
        % Inversion left/right legs
        orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2]; 
        % Qs and Qdots are inverted left/right after a half gait cycle 
        % For pelvis: pelvis tilt and pelvis ty should be equal
        % For trunk: lumbar ext should be equal
        orderQsInv = [
            2*jointi.pelvis.tilt-1:2*jointi.pelvis.ty,...
            2*jointi.hip.r-1:2*jointi.hip.r,...
            2*jointi.hip.l-1:2*jointi.hip.l,...
            2*jointi.knee.r-1:2*jointi.knee.r,...
            2*jointi.knee.l-1:2*jointi.knee.l,...
            2*jointi.ankle.r-1:2*jointi.ankle.r,...
            2*jointi.ankle.l-1:2*jointi.ankle.l,...
            2*jointi.trunk.ext-1:2*jointi.trunk.ext];   
        if k ~= N-1
            % Muscle activations
            ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
            w               = {w{:}, ak};
            lbw             = [lbw; bounds.a.lower'];
            ubw             = [ubw; bounds.a.upper'];
            w0              = [w0;  guess.a(k+2,:)'];
            % Tendon forces
            FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
            w               = {w{:}, FTtildek};
            lbw             = [lbw; bounds.FTtilde.lower'];
            ubw             = [ubw; bounds.FTtilde.upper'];
            w0              = [w0;  guess.FTtilde(k+2,:)'];    
            % Joint positions and velocities   
            w               = {w{:}, Xk{k+2,1}};
            lbw             = [lbw; bounds.QsQdots.lower'];
            ubw             = [ubw; bounds.QsQdots.upper']; 
            w0              = [w0;  guess.QsQdots(k+2,:)'];
            % Back activations
            a_bk            = MX.sym(['a_b_' num2str(k+1)], nq.trunk);
            w               = {w{:}, a_bk};
            lbw             = [lbw; bounds.a_b.lower'];
            ubw             = [ubw; bounds.a_b.upper'];
            w0              = [w0;  guess.a_b(k+2,:)'];
        else
            % Muscle activations
            ak              = MX.sym(['a_' num2str(k+1)], NMuscle);
            w               = {w{:}, ak};
            lbw             = [lbw; bounds.a.lower'];
            ubw             = [ubw; bounds.a.upper'];            
            w0              = [w0;  guess.a(1,orderMusInv)'];
            % Tendon forces
            FTtildek        = MX.sym(['FTtilde_' num2str(k+1)], NMuscle);
            w               = {w{:}, FTtildek};
            lbw             = [lbw; bounds.FTtilde.lower'];
            ubw             = [ubw; bounds.FTtilde.upper'];
            w0              = [w0;  guess.FTtilde(1,orderMusInv)'];    
            % Joint positions and velocities   
            w               = {w{:}, Xk{k+2,1}};
            lbw             = [lbw; bounds.QsQdots.lower'];
            ubw             = [ubw; bounds.QsQdots.upper'];
            inv_X           = guess.QsQdots(1,orderQsInv);     
            dx = guess.QsQdots(end,2*jointi.pelvis.tx-1) - ...
                guess.QsQdots(end-1,2*jointi.pelvis.tx-1);
            inv_X(2*jointi.pelvis.tx-1) = ...
                guess.QsQdots(end,2*jointi.pelvis.tx-1) + dx;            
            w0                = [w0;  inv_X'];   
            % Back activations
            a_bk            = MX.sym(['a_b_' num2str(k+1)], nq.trunk);
            w               = {w{:}, a_bk};
            lbw             = [lbw; bounds.a_b.lower'];
            ubw             = [ubw; bounds.a_b.upper'];            
            w0              = [w0;  guess.a_b(1,:)'];
        end
        % Rescale for equality constraint
        Xk_end = (Xk_nsc_end)./scaling.QsQdots';
        FTtildek_end = (FTtildek_nsc_end)./scaling.FTtilde';
        % Add equality constraints (next interval starts with end values of 
        % states from previous interval).
        g   = {g{:}, Xk_end-Xk{k+2,1}, FTtildek_end-FTtildek, ...
            ak_end-ak, a_bk_end-a_bk};
        lbg = [lbg; zeros(2*nq.all + NMuscle + NMuscle + nq.trunk,1)];
        ubg = [ubg; zeros(2*nq.all + NMuscle + NMuscle + nq.trunk,1)];    
    end
    
    % Impose periodicity at the end: Qs, Qdots, a, FTtilde, a_b          
    QsInvA = [jointi.pelvis.tilt:2*jointi.pelvis.tilt,...
        2*jointi.pelvis.tx,2*jointi.pelvis.ty-1:2*jointi.pelvis.ty,...
        2*jointi.hip.l-1:2*jointi.trunk.ext]';
    QsInvB = [jointi.pelvis.tilt:2*jointi.pelvis.tilt,...
        2*jointi.pelvis.tx,2*jointi.pelvis.ty-1:2*jointi.pelvis.ty,...
        2*jointi.hip.r-1:2*jointi.hip.r,...
        2*jointi.hip.l-1:2*jointi.hip.l,...
        2*jointi.knee.r-1:2*jointi.knee.r,...
        2*jointi.knee.l-1:2*jointi.knee.l,...
        2*jointi.ankle.r-1:2*jointi.ankle.r,...
        2*jointi.ankle.l-1:2*jointi.ankle.l,...
        2*jointi.trunk.ext-1,2*jointi.trunk.ext]';             
    g   = {g{:}, Xk_end(QsInvA)-X0(QsInvB,1)};
    lbg = [lbg; zeros(length(QsInvB),1)];
    ubg = [ubg; zeros(length(QsInvB),1)];         
    % Impose periodicity at the end: muscle activations
    g   = {g{:}, ak_end-a0(orderMusInv,1)};
    lbg = [lbg; zeros(NMuscle,1)];
    ubg = [ubg; zeros(NMuscle,1)];
    % Impose periodicity at the end: tendon forces
    g   = {g{:}, FTtildek_end-FTtilde0(orderMusInv,1)};
    lbg = [lbg; zeros(NMuscle,1)];
    ubg = [ubg; zeros(NMuscle,1)];
    % Impose periodicity at the end: back activations
    g   = {g{:}, a_bk_end-a_b0(1,1)};
    lbg = [lbg; zeros(nq.trunk,1)];
    ubg = [ubg; zeros(nq.trunk,1)];
    
    % Impose average speed
    vel_aver_tot = dist_trav_tot/tf; 
    g   = {g{:}, vel_aver_tot - v_tgt};
    lbg = [lbg; 0];
    ubg = [ubg; 0];  
   
    %% Solve problem
    % Create an NLP solver
    prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
    options.ipopt.hessian_approximation = 'limited-memory';
    options.ipopt.mu_strategy      = 'adaptive';
    options.ipopt.max_iter = 10000;
    options.ipopt.tol = 1*10^(-tol_ipopt);
    switch setup.derivatives
        case 'FD_F' 
            options.common_options.helper_options = ...
                struct('enable_fd',true,'enable_forward',false,...
                'enable_reverse',false,'print_in',false,...
                'fd_method','forward');
    end
    solver = nlpsol('solver', 'ipopt', prob,options);
    % Create and save diary
    p = mfilename('fullpath');
    [~,namescript,~] = fileparts(p);
    [pathRepo,~,~] = fileparts(pathmain);
    pathresults = [pathRepo,'\Results'];
    if ~(exist([pathresults,'\',namescript],'dir')==7)
        mkdir(pathresults,namescript);
    end
    if (exist([pathresults,'\',namescript,'\Diary',savename],'file')==2)
        delete ([pathresults,'\',namescript,'\Diary',savename])
    end 
    diary([pathresults,'\',namescript,'\Diary',savename]);  
    
    sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
        'lbg', lbg, 'ubg', ubg);    
    diary off
    w_opt = full(sol.x);
    g_opt = full(sol.g);  
    % Save results
    setup.tolerance.ipopt = tol_ipopt;
    setup.bounds = bounds;
    setup.scaling = scaling;
    setup.guess = guess;
    setup.lbw = lbw;
    setup.ubw = ubw;
    % Save results
    save([pathresults,'\',namescript,'\w_opt',savename],'w_opt');
    save([pathresults,'\',namescript,'\g_opt',savename],'g_opt');
    save([pathresults,'\',namescript,'\setup',savename],'setup');
end

if analyseResults
    %% Reorganize results
    if loadResults
        p = mfilename('fullpath');
        [~,namescript,~] = fileparts(p);
        pathresults = [pathRepo,'\Results'];
        load([pathresults,'\',namescript,'\w_opt',savename]);
        load([pathresults,'\',namescript,'\g_opt',savename]);
        load([pathresults,'\',namescript,'\setup',savename]);
    end
    
    % Number of design variables
    NControls = NMuscle+NMuscle+nq.all+nGRF+nq.trunk;
    NStates = NMuscle+NMuscle+2*nq.all+nq.trunk;
    NParameters = 1;

    % Static parameters
    tf_opt  = w_opt(1:NParameters);
    % In the loop
    Nwl = NControls+d*(NStates)+NStates;
    % In total
    Nw = NParameters+NStates+N*Nwl;
    % Before the variable corresponding to the first collocation point
    Nwm = NParameters+NStates+NControls;
    
    % Control points
    % Muscle activations and tendon forces
    a_opt = zeros(N+1,NMuscle);
    FTtilde_opt = zeros(N+1,NMuscle);
    for i = 1:NMuscle
        a_opt(:,i) = w_opt(NParameters+i:Nwl:Nw);
        FTtilde_opt(:,i) = w_opt(NParameters+NMuscle+i:Nwl:Nw);
    end
    % Joint positions and velocities (Qs and Qdots)
    q_opt = zeros(N+1,nq.all);
    qdot_opt = zeros(N+1,nq.all);
    count = 0;
    for i = 1:2:2*nq.all
        count = count +1;
        q_opt(:,count) = w_opt(NParameters+NMuscle+NMuscle+i:Nwl:Nw);
        qdot_opt(:,count) = w_opt(NParameters+NMuscle+NMuscle+i+1:Nwl:Nw);
    end
    % Back activations
    a_b_opt = zeros(N+1,nq.trunk);
    for i = 1:nq.trunk
        a_b_opt(:,i) =w_opt(NParameters+NMuscle+NMuscle+2*nq.all+i:Nwl:Nw);
    end    
    % Time derivatives of muscle activations and tendon forces
    vA_opt = zeros(N,NMuscle);
    dFTtilde_opt = zeros(N,NMuscle);
    for i = 1:NMuscle
        vA_opt(:,i) = w_opt(NParameters+NStates+i:Nwl:Nw);
        dFTtilde_opt(:,i) = w_opt(NParameters+NStates+NMuscle+i:Nwl:Nw);
    end
    % Time derivatives of joint velocities
    qdotdot_opt = zeros(N,nq.all);
    for i = 1:nq.all
        qdotdot_opt(:,i) = ...
            w_opt(NParameters+NStates+NMuscle+NMuscle+i:Nwl:Nw);
    end
    % GRF
    GRF_opt = zeros(N,nGRF);
    for i = 1:nGRF
        GRF_opt(:,i) = ...
            w_opt(NParameters+NStates+NMuscle+NMuscle+nq.all+i:Nwl:Nw);
    end
    % Back excitations
    e_b_opt = zeros(N,nq.trunk);
    for i = 1:nq.trunk
        e_b_opt(:,i) = ...
           w_opt(NParameters+NStates+NMuscle+NMuscle+nq.all+nGRF+i:Nwl:Nw); 
    end    
    % Collocation points
    % Muscle activations
    a_opt_ext=zeros(N*(d+1)+1,NMuscle);
    a_opt_ext(1:(d+1):end,:)= a_opt;
    for nmusi=1:NMuscle
        a_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+nmusi:Nwl:Nw);
        a_opt_ext(3:(d+1):end,nmusi) = ...
            w_opt(Nwm+NMuscle+nmusi:Nwl:Nw);
        a_opt_ext(4:(d+1):end,nmusi) = ...
            w_opt(Nwm+NMuscle+NMuscle+nmusi:Nwl:Nw);
    end  
    % Muscle-tendon forces
    FTtilde_opt_ext=zeros(N*(d+1)+1,NMuscle);
    FTtilde_opt_ext(1:(d+1):end,:)= FTtilde_opt;
    for nmusi=1:NMuscle
        FTtilde_opt_ext(2:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle+nmusi:Nwl:Nw);
        FTtilde_opt_ext(3:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle+NMuscle+nmusi:Nwl:Nw);
        FTtilde_opt_ext(4:(d+1):end,nmusi) = ...
            w_opt(Nwm+d*NMuscle+NMuscle+NMuscle+nmusi:Nwl:Nw);
    end
    % Qs & Qdots
    q_opt_ext=zeros(N*(d+1)+1,nq.all);
    q_opt_ext(1:(d+1):end,:)= q_opt;
    q_dot_opt_ext=zeros(N*(d+1)+1,nq.all);
    q_dot_opt_ext(1:(d+1):end,:)= qdot_opt;
    nqi_col = 1:2:2*nq.all;
    for nqi=1:nq.all
        nqi_q = nqi_col(nqi);
        q_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+nqi_q:Nwl:Nw);   
        q_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+nqi_q:Nwl:Nw);  
        q_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+2*nq.all+nqi_q:Nwl:Nw);  
        q_dot_opt_ext(2:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+nqi_q+1:Nwl:Nw);   
        q_dot_opt_ext(3:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+nqi_q+1:Nwl:Nw);  
        q_dot_opt_ext(4:(d+1):end,nqi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+2*nq.all+2*nq.all+nqi_q+1:Nwl:Nw);
    end
    % Back activations
    a_b_opt_ext=zeros(N*(d+1)+1,nq.trunk);
    a_b_opt_ext(1:(d+1):end,:)= a_b_opt;
    for nmusi=1:nq.trunk
        a_b_opt_ext(2:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+d*2*nq.all+nmusi:Nwl:Nw);
        a_b_opt_ext(3:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+d*2*nq.all+nq.trunk+nmusi:Nwl:Nw);
        a_b_opt_ext(4:(d+1):end,nmusi) = w_opt(Nwm+d*NMuscle+...
            d*NMuscle+d*2*nq.all+nq.trunk+nq.trunk+nmusi:Nwl:Nw);
    end   
    
    %% Unscale variables
    % Control points
    % States
    % Qs
    q_opt_unsc.rad = q_opt(1:end-1,:).*repmat(...
    scaling.Qs,size(q_opt(1:end-1,:),1),1); 
    q_opt_unsc.deg = q_opt_unsc.rad;
    dof_roti = [jointi.pelvis.tilt,jointi.hip.l:jointi.trunk.ext];
    q_opt_unsc.deg(:,dof_roti) = q_opt_unsc.deg(:,dof_roti).*180/pi;  
    q_opt_unsc_all.rad = q_opt(1:end,:).*repmat(...
        scaling.Qs,size(q_opt(1:end,:),1),1); 
    q_opt_unsc_all.deg = q_opt_unsc_all.rad;
    q_opt_unsc_all.deg(:,dof_roti) = ...
        q_opt_unsc_all.deg(:,dof_roti).*180/pi;    
    % Qdots
    qdot_opt_unsc.rad = qdot_opt(1:end-1,:).*repmat(...
        scaling.Qdots,size(qdot_opt(1:end-1,:),1),1);
    qdot_opt_unsc.deg = qdot_opt_unsc.rad;
    qdot_opt_unsc.deg(:,dof_roti) = ...
        qdot_opt_unsc.deg(:,dof_roti).*180/pi;    
    qdot_opt_unsc_all.rad = qdot_opt(1:end,:).*repmat(...
        scaling.Qdots,size(qdot_opt(1:end,:),1),1); 
    % a
    a_opt_unsc = a_opt(1:end-1,:).*repmat(...
        scaling.a,size(a_opt(1:end-1,:),1),size(a_opt,2));
    % FTtilde
    FTtilde_opt_unsc = FTtilde_opt(1:end-1,:).*repmat(...
        scaling.FTtilde,size(FTtilde_opt(1:end-1,:),1),1);
    % a_b
    a_b_opt_unsc = a_b_opt(1:end-1,:).*repmat(...
        scaling.a_b,size(a_b_opt(1:end-1,:),1),size(a_b_opt,2));
    % Controls
    % Qdotdots
    qdotdot_opt_unsc.rad = ...
        qdotdot_opt.*repmat(scaling.Qdotdots,size(qdotdot_opt,1),1);
    qdotdot_opt_unsc.deg = qdotdot_opt_unsc.rad;
    qdotdot_opt_unsc.deg(:,dof_roti) = ...
        qdotdot_opt_unsc.deg(:,dof_roti).*180/pi;
    % vA
    vA_opt_unsc = vA_opt.*repmat(scaling.vA,size(vA_opt,1),size(vA_opt,2));
    tact = 0.015;
    tdeact = 0.06;
    e_opt_unsc = computeExcitationRaasch(a_opt_unsc,vA_opt_unsc,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);
    % dFTtilde
    dFTtilde_opt_unsc = dFTtilde_opt.*repmat(...
        scaling.dFTtilde,size(dFTtilde_opt,1),size(dFTtilde_opt,2));
    % GRF
    GRF_opt_unsc = GRF_opt.*repmat(scaling.GRF,size(GRF_opt,1),1);
    % e_b
    e_b_opt_unsc = e_b_opt.*repmat(scaling.e_b,size(e_b_opt,1),...
        size(e_b_opt,2));    
    
    %% Grid    
    % control points
    tgrid = linspace(0,tf_opt,N+1); % mesh points 
    dtime = zeros(1,d+1);
    for i=1:4
        dtime(i)=tau_root(i)*(tf_opt/N);
    end
    % control points and collocation points
    tgrid_ext = zeros(1,(d+1)*N+1);
    for i=1:N
        tgrid_ext(((i-1)*4+1):1:i*4)=tgrid(i)+dtime;
    end
    tgrid_ext(end)=tf_opt;     
 
    %% Residuals and GRF at optimal solution
    Xk_Qs_Qdots_opt             = zeros(N,2*nq.all);
    Xk_Qs_Qdots_opt(:,1:2:end)  = q_opt_unsc.rad;
    Xk_Qs_Qdots_opt(:,2:2:end)  = qdot_opt_unsc.rad;
    Xk_Qdotdots_opt             = qdotdot_opt_unsc.rad;
    out_res_opt = zeros(N,nq.all+nGRF);
    for i = 1:N
        [res] = F([Xk_Qs_Qdots_opt(i,:)';Xk_Qdotdots_opt(i,:)']);
        out_res_opt(i,:) = full(res);    
    end
    assertGRFmax = max(max(abs(out_res_opt(:,GRFi.all)-GRF_opt_unsc)));
    assertBackTmax = max(max(abs(out_res_opt(:,jointi.trunk.ext)-...
        (a_b_opt_unsc)*scaling.BackTau)));
       
    %% Average speed
    % distance traveled
    dist_trav_opt = q_opt_ext(end,jointi.pelvis.tx)*...
        scaling.Qs(jointi.pelvis.tx) - q_opt_ext(1,jointi.pelvis.tx)*...
        scaling.Qs(jointi.pelvis.tx);
    % time elapsed
    time_elaps_opt = tf_opt;
    % Averaged velocity
    vel_aver_opt = dist_trav_opt/time_elaps_opt;     
   
    %% Reconstruct full gait cycle
if reconstructGaitCycle                
    % Different phases
    threshold = 20;
    if exist('HS1','var')
        clear HS1
    end
    % Right heel strike first    
    phase_tran_tgridi = find(GRF_opt_unsc(:,2)<threshold,1,'last');
    if ~isempty(phase_tran_tgridi)        
        if phase_tran_tgridi == N
            temp_idx = find(GRF_opt_unsc(:,2)>threshold,1,'first');
            if ~isempty(temp_idx)
                if temp_idx-1 ~= 0 && ...
                        find(GRF_opt_unsc(temp_idx-1,2)<threshold)
                    phase_tran_tgridi_t = temp_idx;             
                    IC1i = phase_tran_tgridi_t;
                    HS1 = 'r';
                end 
            else            
                IC1i = phase_tran_tgridi + 1; 
                HS1 = 'r';
            end
        else            
            IC1i = phase_tran_tgridi + 1; 
            HS1 = 'r';
        end        
    end
    if ~exist('HS1','var')
        % Left heel strike first
        phase_tran_tgridi = find(GRF_opt_unsc(:,4)<threshold,1,'last');       
        if phase_tran_tgridi == N
            temp_idx = find(GRF_opt_unsc(:,4)>threshold,1,'first');
            if ~isempty(temp_idx)  
                if temp_idx-1 ~= 0 && ...
                        find(GRF_opt_unsc(temp_idx-1,4)<threshold)
                    phase_tran_tgridi_t = temp_idx;             
                    IC1i = phase_tran_tgridi_t;
                    HS1 = 'l';
                else
                    IC1i = phase_tran_tgridi + 1; 
                    HS1 = 'l';
                end 
            else
                IC1i = phase_tran_tgridi + 1; 
                HS1 = 'l';
            end
        else            
            IC1i = phase_tran_tgridi + 1; 
            HS1 = 'l';
        end        
    end
    if isempty(phase_tran_tgridi)
        continue;
    end

    % Joint kinematics
    QsSymA = [jointi.pelvis.tilt,jointi.pelvis.ty,...
        jointi.hip.l:jointi.trunk.ext];
    QsSymB = [jointi.pelvis.tilt,jointi.pelvis.ty,...    
        jointi.hip.r,...
        jointi.hip.l,...
        jointi.knee.r,jointi.knee.l,...
        jointi.ankle.r,jointi.ankle.l,...
        jointi.trunk.ext]; 
    QsSymA_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
        jointi.pelvis.ty,...
        jointi.hip.l:jointi.trunk.ext];
    QsSymB_ptx = [jointi.pelvis.tilt,jointi.pelvis.tx,...
        jointi.pelvis.ty,...    
        jointi.hip.r,...
        jointi.hip.l,...
        jointi.knee.r,jointi.knee.l,...
        jointi.ankle.r,jointi.ankle.l,...
        jointi.trunk.ext]; 
    % Joint positions
    q_opt_GC = zeros(N*2,size(q_opt_unsc.deg,2));
    q_opt_GC(1:N-IC1i+1,:) = q_opt_unsc.deg(IC1i:end,:);   
    q_opt_GC(N-IC1i+2:N-IC1i+1+N,QsSymA) = q_opt_unsc.deg(1:end,QsSymB);
    q_opt_GC(N-IC1i+2:N-IC1i+1+N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:end,jointi.pelvis.tx) + ...
        q_opt_unsc_all.deg(end,jointi.pelvis.tx);        
    q_opt_GC(N-IC1i+2+N:2*N,:) = q_opt_unsc.deg(1:IC1i-1,:);    
    q_opt_GC(N-IC1i+2+N:2*N,jointi.pelvis.tx) = ...
        q_opt_unsc.deg(1:IC1i-1,jointi.pelvis.tx) + ...
        2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')         
        q_opt_GC(:,QsSymA_ptx)  = q_opt_GC(:,QsSymB_ptx);
    end  
    temp_q_opt_GC_pelvis_tx = q_opt_GC(1,jointi.pelvis.tx);
    q_opt_GC(:,jointi.pelvis.tx) = q_opt_GC(:,jointi.pelvis.tx)-...
        temp_q_opt_GC_pelvis_tx;
    q_opt_GUI_GC = zeros(2*N,1+nq.all);
    q_opt_GUI_GC(1:N-IC1i+1,1) = tgrid(:,IC1i:end-1)';
    q_opt_GUI_GC(N-IC1i+2:N-IC1i+1+N,1)  = tgrid(:,1:end-1)' + tgrid(end);
    q_opt_GUI_GC(N-IC1i+2+N:2*N,1) = tgrid(:,1:IC1i-1)' + 2*tgrid(end);    
    q_opt_GUI_GC(:,2:end) = q_opt_GC;
    q_opt_GUI_GC(:,1) = q_opt_GUI_GC(:,1)-q_opt_GUI_GC(1,1);
    pathOpenSim = [pathRepo,'\OpenSim'];
    addpath(genpath(pathOpenSim));
    if writeIKmotion
        JointAngle.labels = {'time','pelvis_tilt','pelvis_tx','pelvis_ty',...
        'hip_flexion_l','hip_flexion_r',...
        'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
        'lumbar_extension'};
    q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
    q_opt_GUI_GC_2(2*N+1:4*N,1) = q_opt_GUI_GC_2(2*N+1:4*N,1) + ...
        q_opt_GUI_GC_2(end,1) + ...
        q_opt_GUI_GC_2(end,1)-q_opt_GUI_GC_2(end-1,1);
    q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) = ...
        q_opt_GUI_GC_2(2*N+1:4*N,jointi.pelvis.tx+1) + ...
        2*q_opt_unsc_all.deg(end,jointi.pelvis.tx);
    JointAngle.data = q_opt_GUI_GC_2;
    filenameJointAngles = [pathRepo,'\Results\',namescript,...
            '\IK',savename,'.mot'];
    write_motionFile(JointAngle, filenameJointAngles)
    end
    % Joint velocities
    qdot_opt_GC = zeros(N*2,size(q_opt,2));
    qdot_opt_GC(1:N-IC1i+1,:) = qdot_opt_unsc.deg(IC1i:end,:);
    qdot_opt_GC(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
        qdot_opt_unsc.deg(1:end,QsSymB_ptx);
    qdot_opt_GC(N-IC1i+2+N:2*N,:) = ...
        qdot_opt_unsc.deg(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        qdot_opt_GC(:,QsSymA_ptx) = qdot_opt_GC(:,QsSymB_ptx);
    end
    % Joint accelerations
    qdotdot_opt_GC = zeros(N*2,size(q_opt,2));
    qdotdot_opt_GC(1:N-IC1i+1,:) = qdotdot_opt_unsc.deg(IC1i:end,:);
    qdotdot_opt_GC(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
        qdotdot_opt_unsc.deg(1:end,QsSymB_ptx);
    qdotdot_opt_GC(N-IC1i+2+N:2*N,:) = qdotdot_opt_unsc.deg(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        qdotdot_opt_GC(:,QsSymA_ptx) = qdotdot_opt_GC(:,QsSymB_ptx);
    end

    % Ground reaction forces
    GRF_opt_GC = zeros(N*2,nGRF);
    GRF_opt_GC(1:N-IC1i+1,:) = GRF_opt_unsc(IC1i:end,:);
    GRF_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = GRF_opt_unsc(1:end,[3:4,1:2]);
    GRF_opt_GC(N-IC1i+2+N:2*N,:) = GRF_opt_unsc(1:IC1i-1,:);
    GRF_opt_GC = GRF_opt_GC./(body_weight/100);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        GRF_opt_GC(:,[3:4,1:2]) = GRF_opt_GC(:,:);
    end  
    
    % Joint kinetics
    tau_opt_GC = zeros(N*2,size(q_opt,2));
    tau_opt_GC(1:N-IC1i+1,1:nq.all) = ...
        out_res_opt(IC1i:end,1:nq.all)./body_mass;
    tau_opt_GC(N-IC1i+2:N-IC1i+1+N,QsSymA_ptx) = ...
        out_res_opt(1:end,QsSymB_ptx)./body_mass;
    tau_opt_GC(N-IC1i+2+N:2*N,1:nq.all) = ...
        out_res_opt(1:IC1i-1,1:nq.all)./body_mass;
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        tau_opt_GC(:,QsSymA_ptx) = tau_opt_GC(:,QsSymB_ptx);
    end

    % Tendon Forces
    orderMusInv = [NMuscle/2+1:NMuscle,1:NMuscle/2];
    FTtilde_opt_GC = zeros(N*2,NMuscle);
    FTtilde_opt_GC(1:N-IC1i+1,:) = FTtilde_opt_unsc(IC1i:end,:);
    FTtilde_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
        FTtilde_opt_unsc(1:end,orderMusInv);
    FTtilde_opt_GC(N-IC1i+2+N:2*N,:) = FTtilde_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        FTtilde_opt_GC(:,:) = FTtilde_opt_GC(:,orderMusInv);
    end

    % Muscle activations
    a_opt_GC = zeros(N*2,NMuscle);
    a_opt_GC(1:N-IC1i+1,:) = a_opt_unsc(IC1i:end,:);
    a_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = a_opt_unsc(1:end,orderMusInv);
    a_opt_GC(N-IC1i+2+N:2*N,:) = a_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        a_opt_GC(:,:) = a_opt_GC(:,orderMusInv);
    end

    % Derivative of tendon forces
    dFTtilde_opt_GC = zeros(N*2,NMuscle);
    dFTtilde_opt_GC(1:N-IC1i+1,:) = dFTtilde_opt_unsc(IC1i:end,:);
    dFTtilde_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = ...
        dFTtilde_opt_unsc(1:end,orderMusInv);
    dFTtilde_opt_GC(N-IC1i+2+N:2*N,:) = dFTtilde_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        dFTtilde_opt_GC(:,:) = dFTtilde_opt_GC(:,orderMusInv);
    end

    % Muscle excitations
    vA_opt_GC = zeros(N*2,NMuscle);
    vA_opt_GC(1:N-IC1i+1,:) = vA_opt_unsc(IC1i:end,:);
    vA_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = vA_opt_unsc(1:end,orderMusInv);
    vA_opt_GC(N-IC1i+2+N:2*N,:) = vA_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        vA_opt_GC(:,:) = vA_opt_GC(:,orderMusInv);
    end
    tact = 0.015;
    tdeact = 0.06;
    e_opt_GC = computeExcitationRaasch(a_opt_GC,vA_opt_GC,...
        ones(1,NMuscle)*tdeact,ones(1,NMuscle)*tact);   
    
    % Back activations
    a_b_opt_GC = zeros(N*2,nq.trunk);
    a_b_opt_GC(1:N-IC1i+1,:) = a_b_opt_unsc(IC1i:end,:);
    a_b_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = a_b_opt_unsc(1:end,:);
    a_b_opt_GC(N-IC1i+2+N:2*N,:) = a_b_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        a_b_opt_GC(:,:) = a_b_opt_GC(:,:);
    end

    % Back excitations
    e_b_opt_GC = zeros(N*2,nq.trunk);
    e_b_opt_GC(1:N-IC1i+1,:) = e_b_opt_unsc(IC1i:end,:);
    e_b_opt_GC(N-IC1i+2:N-IC1i+1+N,:) = e_b_opt_unsc(1:end,:);
    e_b_opt_GC(N-IC1i+2+N:2*N,:) = e_b_opt_unsc(1:IC1i-1,:);
    % If the first heel strike was on left foot then we invert so that we
    % always start with the right foot, easier for comparison afterwards
    if strcmp(HS1,'l')
        e_b_opt_GC(:,:) = e_b_opt_GC(:,:);
    end
end  
%% Save sensitivity results       
    if saveSensitivityResults
        if (exist([pathresults,'\',namescript,...
                '\SensitivityResults.mat'],'file')==2) 
            load([pathresults,'\',namescript,'\SensitivityResults.mat']);
        else
            SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
                (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
                (['NMesh_',num2str(N)]). ...
                (['Wact_',num2str(W.act*10)]). ...
                (['Wback_',num2str(W.back*10)]). ...
                (['WGRF_',num2str(W.GRF*10)]). ...            
                (['Wacc_',num2str(W.acc*10)]). ...                
                (['expAct_',num2str(exp_act)]). ...                
                (['Derivatives_',setup.derivatives]) = struct('q_opt',[]);
        end    
        % Put data in structure
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).q_opt_GC = q_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).qdot_opt_GC = ...
            qdot_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).qdotdot_opt_GC = ...
            qdotdot_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).GRF_opt_GC = ...
            GRF_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).tau_opt_GC = ...
            tau_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).a_opt_GC = ...
            a_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).FTtilde_opt_GC = ...
            FTtilde_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).e_opt_GC = e_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).dFTtilde_opt_GC = ...
            dFTtilde_opt_GC;
        SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
            (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
            (['NMesh_',num2str(N)]). ...
            (['Wact_',num2str(W.act*10)]). ...
            (['Wback_',num2str(W.back*10)]). ...
            (['WGRF_',num2str(W.GRF*10)]). ...            
            (['Wacc_',num2str(W.acc*10)]). ...                
            (['expAct_',num2str(exp_act)]). ...                
            (['Derivatives_',setup.derivatives]).a_b_opt_GC = ...
            a_b_opt_GC; 
        % Save data
        save([pathresults,'\',namescript,'\SensitivityResults.mat'],...
        'SensitivityResults');
    end
end
end