% This script provides bounds and scaling factors for the design variables

function [bounds,scaling] = getBounds_gait1018(NMuscle,nq,jointi)

%% Bounds
% Qs
bounds.Qs.upper = [30,3,1.4,70,70,10,10,45,45,30];
bounds.Qs.lower = -bounds.Qs.upper;
bounds.Qs.lower(jointi.pelvis.tx) = 0;
bounds.Qs.lower(jointi.pelvis.ty) = 0.75;
bounds.Qs.lower([jointi.knee.l,jointi.knee.r]) = -120;
dof_roti = [jointi.pelvis.tilt,jointi.hip.l:jointi.trunk.ext];
bounds.Qs.upper(dof_roti) = bounds.Qs.upper(dof_roti)*pi/180;
bounds.Qs.lower(dof_roti) = bounds.Qs.lower(dof_roti)*pi/180;
% Qdots
bounds.Qdots.upper = [5,12,2,10,10,20,20,10,10,5];
bounds.Qdots.lower = -bounds.Qdots.upper;
bounds.Qdots.lower(jointi.pelvis.tx) = 0;
% Qdotdots
bounds.Qdotdots.upper = [50,30,30,150,150,250,250,250,250,50];
bounds.Qdotdots.lower = -bounds.Qdotdots.upper;
% Muscle activations
bounds.a.lower = zeros(1,NMuscle);
bounds.a.upper = ones(1,NMuscle);
% Tendon forces
bounds.FTtilde.lower = zeros(1,NMuscle);
bounds.FTtilde.upper = 5*ones(1,NMuscle);
% Derivatives of muscle activations
tact = 0.015;
tdeact = 0.06;
bounds.vA.lower = (-1/100*ones(1,NMuscle))./(ones(1,NMuscle)*tdeact);
bounds.vA.upper = (1/100*ones(1,NMuscle))./(ones(1,NMuscle)*tact);
% Derivatives of tendon forces
bounds.dFTtilde.lower = -1*ones(1,NMuscle);
bounds.dFTtilde.upper = 1*ones(1,NMuscle);
% GRF
bounds.GRF.upper = [500,2000,500,2000];
bounds.GRF.lower = [-500,-200,-500,-200];
% Back torque activations
bounds.a_b.lower = -ones(1,nq.trunk);
bounds.a_b.upper = ones(1,nq.trunk);
% Back torque excitations
bounds.e_b.lower = -ones(1,nq.trunk);
bounds.e_b.upper = ones(1,nq.trunk);

%% Scaling
% Qs
scaling.Qs      = max(abs(bounds.Qs.lower),abs(bounds.Qs.upper));
bounds.Qs.lower = (bounds.Qs.lower)./scaling.Qs;
bounds.Qs.upper = (bounds.Qs.upper)./scaling.Qs;
% Qdots
scaling.Qdots      = max(abs(bounds.Qdots.lower),abs(bounds.Qdots.upper));
bounds.Qdots.lower = (bounds.Qdots.lower)./scaling.Qdots;
bounds.Qdots.upper = (bounds.Qdots.upper)./scaling.Qdots;
% Qs and Qdots are intertwined
bounds.QsQdots.lower = zeros(1,2*nq.all);
bounds.QsQdots.upper = zeros(1,2*nq.all);
bounds.QsQdots.lower(1,1:2:end) = bounds.Qs.lower;
bounds.QsQdots.upper(1,1:2:end) = bounds.Qs.upper;
bounds.QsQdots.lower(1,2:2:end) = bounds.Qdots.lower;
bounds.QsQdots.upper(1,2:2:end) = bounds.Qdots.upper;
scaling.QsQdots                 = zeros(1,2*nq.all);
scaling.QsQdots(1,1:2:end)      = scaling.Qs ;
scaling.QsQdots(1,2:2:end)      = scaling.Qdots ;
% Qdotdots
scaling.Qdotdots = max(abs(bounds.Qdotdots.lower),...
    abs(bounds.Qdotdots.upper));
bounds.Qdotdots.lower = (bounds.Qdotdots.lower)./scaling.Qdotdots;
bounds.Qdotdots.upper = (bounds.Qdotdots.upper)./scaling.Qdotdots;
bounds.Qdotdots.lower(isnan(bounds.Qdotdots.lower)) = 0;
bounds.Qdotdots.upper(isnan(bounds.Qdotdots.upper)) = 0;
% GRF
scaling.GRF = max(abs(bounds.GRF.lower),abs(bounds.GRF.upper));
bounds.GRF.lower = (bounds.GRF.lower)./scaling.GRF;
bounds.GRF.upper = (bounds.GRF.upper)./scaling.GRF;
% Back torques
scaling.BackTau = 150;
% Derivatives of muscle activations
scaling.vA = 100;
% Muscle activations
scaling.a = 1;
% Back torque activations
scaling.a_b = 1;
% Back torque excitations
scaling.e_b = 1;
% Derivatives of tendon forces
scaling.dFTtilde = 100;
% Tendon forces
scaling.FTtilde         = max(...
    abs(bounds.FTtilde.lower),abs(bounds.FTtilde.upper)); 
bounds.FTtilde.lower    = (bounds.FTtilde.lower)./scaling.FTtilde;
bounds.FTtilde.upper    = (bounds.FTtilde.upper)./scaling.FTtilde;

%% Hard constraint
% Initial position of pelvis_tx is in O
bounds.QsQdots_0.lower = bounds.QsQdots.lower;
bounds.QsQdots_0.upper = bounds.QsQdots.upper;
bounds.QsQdots_0.lower(2*jointi.pelvis.tx-1) = 0;
bounds.QsQdots_0.upper(2*jointi.pelvis.tx-1) = 0;

end
