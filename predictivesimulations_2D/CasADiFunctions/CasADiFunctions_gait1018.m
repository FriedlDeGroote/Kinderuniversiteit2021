% This script contains several CasADi-based functions that are
% used when solving the OCPs
import casadi.*

%% Function for polynomial approximation
pathpolynomial = [pathRepo,'\Polynomials'];
addpath(genpath(pathpolynomial));
muscle_spanning_info_m = muscle_spanning_joint_INFO(musi_pol,:);
MuscleInfo_m.muscle    = MuscleInfo.muscle(musi_pol);                  
qin     = SX.sym('qin',1,nq.leg);
qdotin  = SX.sym('qdotin',1,nq.leg);
lMT     = SX(NMuscle_pol,1);
vMT     = SX(NMuscle_pol,1);
dM      = SX(NMuscle_pol,nq.leg);
for i=1:NMuscle_pol      
    index_dof_crossing  = find(muscle_spanning_info_m(i,:)==1);
    order               = MuscleInfo_m.muscle(i).order;
    [mat,diff_mat_q]    = n_art_mat_3_cas_SX(qin(1,index_dof_crossing),...
                                             order);
    lMT(i,1)            = mat*MuscleInfo_m.muscle(i).coeff;
    vMT(i,1)            = 0;
    dM(i,1:nq.leg)      = 0;
    nr_dof_crossing     = length(index_dof_crossing); 
    for dof_nr = 1:nr_dof_crossing
        dM(i,index_dof_crossing(dof_nr)) = ...
            (-(diff_mat_q(:,dof_nr)))'*MuscleInfo_m.muscle(i).coeff;
        vMT(i,1) = vMT(i,1) + (-dM(i,index_dof_crossing(dof_nr))*...
            qdotin(1,index_dof_crossing(dof_nr)));
    end 
end
f_lMT_vMT_dM = Function('f_lMT_vMT_dM',{qin,qdotin},{lMT,vMT,dM});

%% Functions for sum of values to a certain exponent
etempNMuscle = SX.sym('etempNMuscle',NMuscle);
exp          = SX.sym('exp',1);
JtempNMuscle = 0;
for i=1:length(etempNMuscle)
    JtempNMuscle = JtempNMuscle + etempNMuscle(i).^exp;
end
f_sunsqr_exp=Function('f_sunsqr_exp',{etempNMuscle,exp},{JtempNMuscle});

%% Function for sum of products 
% 3 elements (ankle)
ma_temp3 = SX.sym('ma_temp3',3);
ft_temp3 = SX.sym('ft_temp3',3);
J_sptemp3 = 0;
for i=1:length(ma_temp3)
    J_sptemp3 = J_sptemp3 + ma_temp3(i,1)*ft_temp3(i,1);    
end
f_T3 = Function('f_T3',{ma_temp3,ft_temp3},{J_sptemp3});
% 4 elements (hip)
ma_temp4 = SX.sym('ma_temp4',4);
ft_temp4 = SX.sym('ft_temp4',4);
J_sptemp4 = 0;
for i=1:length(ma_temp4)
    J_sptemp4 = J_sptemp4 + ma_temp4(i,1)*ft_temp4(i,1);    
end
f_T4 = Function('f_T4',{ma_temp4,ft_temp4},{J_sptemp4});
% 5 elements (knee)
ma_temp5 = SX.sym('ma_temp5',5);
ft_temp5 = SX.sym('ft_temp5',5);
J_sptemp5 = 0;
for i=1:length(ma_temp5)
    J_sptemp5 = J_sptemp5 + ma_temp5(i,1)*ft_temp5(i,1);    
end
f_T5 = Function('f_T5',{ma_temp5,ft_temp5},{J_sptemp5});

%% Muscle contraction dynamics
pathmusclemodel = [pathRepo,'\MuscleModel'];
addpath(genpath(pathmusclemodel));
% Controls 
FTtilde = SX.sym('FTtilde',NMuscle); % Normalized tendon force
% States 
dFTtilde = SX.sym('dFTtilde',NMuscle); % Normalized derivative tendon force
a = SX.sym('a',NMuscle); % Muscle activations
% Parameters of force-length-velocity curves
load Fvparam
load Fpparam
load Faparam
Hilldiff    = SX(NMuscle,1);
FT          = SX(NMuscle,1);
lMT         = SX.sym('lMT',NMuscle); % Muscle-tendon length
vMT         = SX.sym('vMT',NMuscle); % Muscle-tendon velocity
atendon_SX  = SX.sym('atendon',NMuscle); % Tendon stiffness
for m = 1:NMuscle
    [Hilldiff(m),FT(m)] = ForceEquilibrium_FtildeState(a(m),FTtilde(m),...
        dFTtilde(m),lMT(m),vMT(m),MTparameters_m(:,m),Fvparam,Fpparam,...
        Faparam);
end
f_ForceEquilibrium_FtildeState = Function(...
    'f_ForceEquilibrium_FtildeState',{a,FTtilde,dFTtilde,lMT,vMT},...
    {Hilldiff,FT});
% Linear tendon - variable tendon stiffness
for m = 1:NMuscle
    [Hilldiff(m),FT(m)] = ForceEquilibrium_FtildeState_Tlin(a(m),FTtilde(m),...
        dFTtilde(m),lMT(m),vMT(m),MTparameters_m(:,m),Fvparam,Fpparam,...
        Faparam,atendon_SX(m));
end
f_ForceEquilibrium_FtildeState_Tlin = Function(...
    'f_ForceEquilibrium_FtildeState_Tlin',{a,FTtilde,dFTtilde,lMT,vMT,atendon_SX},...
    {Hilldiff,FT});
% Function to obtain lM and lMtilde from FTtilde and lMT
lM      = SX(NMuscle,1);
lMtilde = SX(NMuscle,1);
for m = 1:NMuscle
    [lM(m),lMtilde(m)] = FiberLength_TendonForce(FTtilde(m),...
        MTparameters_m(:,m),lMT(m));
end
f_FiberLength_TendonForce = Function('f_FiberLength_Ftilde',...
    {FTtilde,lMT},{lM,lMtilde});
% Function to obtain vM and vMtilde from FTtilde, dFTtilde, lMT, and vMT
vM      = SX(NMuscle,1);
vMtilde = SX(NMuscle,1);
for m = 1:NMuscle
    [vM(m),vMtilde(m)] = FiberVelocity_TendonForce(FTtilde(m),...
        dFTtilde(m),MTparameters_m(:,m),lMT(m),vMT(m));
end
f_FiberVelocity_TendonForce = Function('f_FiberVelocity_Ftilde',...
    {FTtilde,dFTtilde,lMT,vMT},{vM,vMtilde});

%% Back activation dynamics
e_b = SX.sym('e_b',nq.trunk); % back excitations
a_b = SX.sym('a_b',nq.trunk); % back activations
dadt = BackActivationDynamics(e_b,a_b);
f_BackActivationDynamics = ...
    Function('f_BackActivationDynamics',{e_b,a_b},{dadt});
