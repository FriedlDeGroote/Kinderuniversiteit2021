% Plots of bounds and initial guesses

figure()
for i = 1:size(bounds.QsQdots.lower,2)
    subplot(5,4,i)
    plot([1,N],[bounds.QsQdots.upper(:,i),bounds.QsQdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.QsQdots.lower(:,i),bounds.QsQdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.QsQdots(:,i),'k','linewidth',2);
end
figure()
for i = 1:size(bounds.Qdotdots.lower,2)
    subplot(4,3,i)
    plot([1,N],[bounds.Qdotdots.upper(:,i),bounds.Qdotdots.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.Qdotdots.lower(:,i),bounds.Qdotdots.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.Qdotdots(:,i),'k','linewidth',2);
end
figure()
for i = 1:nGRF
    subplot(2,2,i)
    plot([1,N],[bounds.GRF.upper(:,i),bounds.GRF.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.GRF.lower(:,i),bounds.GRF.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.GRF(:,i),'k','linewidth',2);
end
figure()
for i = 1:NMuscle
    subplot(5,4,i)
    plot([1,N],[bounds.a.upper(:,i),bounds.a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.a.lower(:,i),bounds.a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a(:,i),'k','linewidth',2);
end
figure()
for i = 1:NMuscle
    subplot(5,4,i)
    plot([1,N],[bounds.vA.upper(:,i),bounds.vA.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.vA.lower(:,i),bounds.vA.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a(:,i),'k','linewidth',2);
end
figure()
for i = 1:NMuscle
    subplot(5,4,i)
    plot([1,N],[bounds.dFTtilde.upper(:,i),bounds.dFTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.dFTtilde.lower(:,i),bounds.dFTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a(:,i),'k','linewidth',2);
end
figure()
for i = 1:NMuscle
    subplot(5,4,i)
    plot([1,N],[bounds.FTtilde.upper(:,i),bounds.FTtilde.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.FTtilde.lower(:,i),bounds.FTtilde.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a(:,i),'k','linewidth',2);
end
figure()
for i = 1:nq.trunk
    plot([1,N],[bounds.a_a.upper(:,i),bounds.a_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.a_a.lower(:,i),bounds.a_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.a_a(:,i),'k','linewidth',2);
end
figure()
for i = 1:nq.trunk
    plot([1,N],[bounds.e_a.upper(:,i),bounds.e_a.upper(:,i)],...
        'b--','linewidth',2);
    hold on
    plot([1,N],[bounds.e_a.lower(:,i),bounds.e_a.lower(:,i)],...
        'r--','linewidth',2);
    plot(guess.e_a(:,i),'k','linewidth',2);
end
% figure()
% scatter(1,bounds.tf.upper,'b');
% hold on
% scatter(1,bounds.tf.lower,'r');
% scatter(1,guess.tf,'k','linewidth',2);