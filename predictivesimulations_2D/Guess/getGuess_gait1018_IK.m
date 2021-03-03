% This script provides an inital guess for the predictive simulations

function guess = getGuess_gait1018_IK(Qs,N,nq,NMuscle,body_weight,scaling,time_IC,jointi,nGRF)

% Spline approximation
Qs_spline.data = zeros(size(Qs.allfilt));
Qs_spline.data(:,1) = Qs.allfilt(:,1);
Qdots_spline.data = zeros(size(Qs.allfilt));
Qdots_spline.data(:,1) = Qs.allfilt(:,1);
Qdotdots_spline.data = zeros(size(Qs.allfilt));
Qdotdots_spline.data(:,1) = Qs.allfilt(:,1);
for i = 2:size(Qs.allfilt,2)
    Qs.datafiltspline(i) = spline(Qs.allfilt(:,1),Qs.allfilt(:,i));
    [Qs_spline.data(:,i),Qdots_spline.data(:,i),...
        Qdotdots_spline.data(:,i)] = ...
        SplineEval_ppuval(Qs.datafiltspline(i),Qs.allfilt(:,1),1);
end

%% Qs: based on experimental data (inverse kinematics)
% Pelvis tilt
guess.Qs_all.data(:,jointi.pelvis.tilt) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis_tx
guess.Qs_all.data(:,jointi.pelvis.tx) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qs_all.data(:,jointi.pelvis.ty) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Hip flexion
guess.Qs_all.data(:,jointi.hip.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qs_all.data(:,jointi.hip.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Knee angle
guess.Qs_all.data(:,jointi.knee.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qs_all.data(:,jointi.knee.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle angle
guess.Qs_all.data(:,jointi.ankle.l) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qs_all.data(:,jointi.ankle.r) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Trunk extension
guess.Qs_all.data(:,jointi.trunk.ext) = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Interpolation
Qs_time = Qs_spline.data(:,strcmp(Qs.colheaders(1,:),'time'));
time_expi.Qs(1) = find(round(Qs_time,3) == round(time_IC(1),3));
time_expi.Qs(2) = find(round(Qs_time,3) == round(time_IC(2),3));
step = (Qs_time(time_expi.Qs(2))-Qs_time(time_expi.Qs(1)))/(N-1);
interval = Qs_time(time_expi.Qs(1)):step:Qs_time(time_expi.Qs(2));
guess.Qs = interp1(round(Qs_time,4),guess.Qs_all.data,round(interval,4));
guess.Qs(:,jointi.pelvis.tx) = guess.Qs(:,jointi.pelvis.tx) - ....
    guess.Qs(1,jointi.pelvis.tx);

%% Qdots: based on spline approximation of filtered Qs
% Pelvis tilt
guess.Qdots_all.data(:,jointi.pelvis.tilt) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis_tx
guess.Qdots_all.data(:,jointi.pelvis.tx) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdots_all.data(:,jointi.pelvis.ty) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Hip flexion
guess.Qdots_all.data(:,jointi.hip.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdots_all.data(:,jointi.hip.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Knee angle
guess.Qdots_all.data(:,jointi.knee.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdots_all.data(:,jointi.knee.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle angle
guess.Qdots_all.data(:,jointi.ankle.l) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdots_all.data(:,jointi.ankle.r) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Trunk extension
guess.Qdots_all.data(:,jointi.trunk.ext) = Qdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Interpolation
guess.Qdots = interp1(round(Qs_time,4),guess.Qdots_all.data,...
    round(interval,4));

%% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots
% Pelvis tilt
guess.Qdotdots_all.data(:,jointi.pelvis.tilt) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tilt'));
% Pelvis_tx
guess.Qdotdots_all.data(:,jointi.pelvis.tx) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_tx'));
% Pelvis_ty
guess.Qdotdots_all.data(:,jointi.pelvis.ty) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'pelvis_ty'));
% Hip flexion
guess.Qdotdots_all.data(:,jointi.hip.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_l'));
guess.Qdotdots_all.data(:,jointi.hip.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'hip_flexion_r'));
% Knee angle
guess.Qdotdots_all.data(:,jointi.knee.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_l'));
guess.Qdotdots_all.data(:,jointi.knee.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'knee_angle_r'));
% Ankle angle
guess.Qdotdots_all.data(:,jointi.ankle.l) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_l'));
guess.Qdotdots_all.data(:,jointi.ankle.r) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'ankle_angle_r'));
% Trunk extension
guess.Qdotdots_all.data(:,jointi.trunk.ext) = Qdotdots_spline.data(:,strcmp(Qs.colheaders(1,:),'lumbar_extension'));
% Interpolation
guess.Qdotdots = interp1(round(Qs_time,4),guess.Qdotdots_all.data,...
    round(interval,4));

%% Muscle variables
guess.a = 0.1*ones(N,NMuscle);
guess.vA = 0.01*ones(N,NMuscle);
guess.FTtilde = 0.1*ones(N,NMuscle);
guess.dFTtilde = 0.01*ones(N,NMuscle);

%% GRF
guess.GRF = zeros(N,nGRF);
guess.GRF(:,[2,4]) = body_weight/2;

%% Back actuators
guess.a_b = 0.1*ones(N,nq.trunk);
guess.e_b = 0.1*ones(N,nq.trunk);

%% Scaling
guess.QsQdots   = guess.QsQdots./repmat(scaling.QsQdots,N,1);
guess.Qdotdots  = guess.Qdotdots./repmat(scaling.Qdotdots,N,1);
guess.GRF       = guess.GRF./repmat(scaling.GRF,N,1);
guess.a         = (guess.a)./repmat(scaling.a,N,size(guess.a,2));
guess.FTtilde   = (guess.FTtilde)./repmat(scaling.FTtilde,N,1);
guess.vA        = (guess.vA)./repmat(scaling.vA,N,size(guess.vA,2));
guess.dFTtilde  = (guess.dFTtilde)./repmat(scaling.dFTtilde,N,...
    size(guess.dFTtilde,2));

end
