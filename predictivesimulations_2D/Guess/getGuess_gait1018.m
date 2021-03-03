% This script provides an inital guess for the predictive simulations

function guess = getGuess_gait1018(N,nq,NMuscle,body_weight,scaling,v_tgt,jointi,nGRF)

%% Qs
guess.Qs = zeros(N,nq.all);
% Pelvis_tx: the model walks at a constant speed (v_tgt) for 0.6s
guess.Qs(:,jointi.pelvis.tx) = linspace(0,0.6*v_tgt,N);
% Pelvis_ty: the model is standing on the ground.
guess.Qs(:,jointi.pelvis.ty) = 0.9385;

%% Qdots
guess.Qdots = zeros(N,nq.all);
% Pelvis_tx: the model walks at a constant speed (v_tgt)
guess.Qdots(:,jointi.pelvis.tx) = v_tgt;

%% Qs and Qdots are intertwined
guess.QsQdots = zeros(N,2*nq.all);
guess.QsQdots(:,1:2:end) = guess.Qs;
guess.QsQdots(:,2:2:end) = guess.Qdots;

%% Qdotdots
guess.Qdotdots = zeros(N,nq.all);

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
