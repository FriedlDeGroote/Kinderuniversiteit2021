% This script plots the results
% clear all
close all
clc

%% User inputs
% Select combination of parameters (index row in params below)
ww  = 1;
% Order
% Weight muscle activation
% Weight torque actuator activation
% Weight GRF
% Weight metabolic energy
% Number of mesh interval
% Power muscle activation
% Power metabolic energy
% Tolerance IPOPT
% Average walking speed
params = [1.33, 4, 50, 1, 1, 1, 1, 3;];

showMainPlotsOnly = 1;
showlegend = 0;

%% Other settings
% Select setup 
setup.ocp = 'PredSim';    
der = 'AD';
% Pre-allocation
legend_case = cell(1,length(ww));
q_opt_GC_unsc_deg = struct('m',[]);
qdot_opt_GC_unsc_deg = struct('m',[]);
qdotdot_opt_GC_unsc_deg = struct('m',[]);
a_opt_GC = struct('m',[]);
FTtilde_opt_GC_unsc = struct('m',[]);
vA_opt_GC_unsc = struct('m',[]);
e_opt_GC = struct('m',[]);
dFTtilde_opt_GC_unsc = struct('m',[]);
GRF_opt_GC_unsc = struct('m',[]);
tau_opt_GC_unsc = struct('m',[]);
mom_opt_GC = struct('m',[]);
% Loop of cases
for k = 1:length(ww)
    v_tgt       = settings(ww(k),1);    % target velocity
    tol_ipopt   = settings(ww(k),2);    % tolerance (means 1e-(tol_ipopt))
    N           = settings(ww(k),3);    % number of mesh intervals
    W.act       = settings(ww(k),4);    % weight muscle activations
    W.back      = settings(ww(k),5);    % weight back torque excitations
    W.GRF       = settings(ww(k),6);    % weight ground reaction forces
    W.acc       = settings(ww(k),7);    % weight ground reaction forces
    exp_act     = settings(ww(k),8);    % power muscle activations
    
    % Select Ipopt tolerance
    setup.tolerance.ipopt = tol_ipopt;    
    % Select derivative supplier
    setup.derivatives = der;

    % load results
    pathmain = pwd;
    [pathrepo,~,~] = fileparts(pathmain);
    pathresults = [pathrepo,'\Results'];
    load([pathresults,'\',setup.ocp,'\SensitivityResults.mat']);
    
%     legend_case{k} = ['Speedx100_',num2str(v_tgt*100), ...
%             'TolIpopt_',num2str(setup.tolerance.ipopt),...
%             'NMesh_',num2str(N),'ExpAct_',num2str(exp_act),...
%             'WeightActx1000_',num2str(1000*W.a),...
%             'WeightGRFx1000_',num2str(1000*W.GRF),...
%             'WeightEbx1000_',num2str(1000*W.e_b), ...
%             'WeightEx1000_',num2str(1000*W.E), ...
%             'ExpE_',num2str(exp_E), ...
%             'Weightux1000_',num2str(1000*W.u), ...
%             'Derivatives_',setup.derivatives];    
    
    % Unstructure data
    q_opt_GC_unsc_deg(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).q_opt_GC;
    qdot_opt_GC_unsc_deg(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).qdot_opt_GC;
    qdotdot_opt_GC_unsc_deg(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).qdotdot_opt_GC;
    a_opt_GC(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).a_opt_GC;
    FTtilde_opt_GC_unsc(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).FTtilde_opt_GC;
    e_opt_GC(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).e_opt_GC;
    GRF_opt_GC_unsc(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).GRF_opt_GC;
    mom_opt_GC(ww(k)).m = ...
    SensitivityResults.(['Speedx10_',num2str(v_tgt_id*100)]). ...                
    (['TolIpopt_',num2str(setup.tolerance.ipopt)]). ...
    (['NMesh_',num2str(N)]). ...
    (['Wact_',num2str(W.act*10)]). ...
    (['Wback_',num2str(W.back*10)]). ...
    (['WGRF_',num2str(W.GRF*10)]). ...            
    (['Wacc_',num2str(W.acc*10)]). ...                
    (['expAct_',num2str(exp_act)]). ...                
    (['Derivatives_',setup.derivatives]).tau_opt_GC; 
end

% Comparison to reference data    
% Joint kinematics: Qs
% pathReferenceData = [pathrepo,'\ReferenceData'];
% load([pathReferenceData '\KinematicRef.mat']);

RefData_str = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
    'pelvis_tx','pelvis_ty','pelvis_tz','hip_flexion','hip_adduction',...
    'hip_rotation','knee_angle','knee_vv','knee_rot','ankle_angle',...
    'subtalar_angle','mtp_angle','lumbar_extension','lumbar_bending',...
    'lumbar_rotation','arm_flex','arm_add','arm_rot','elbow_flex',...
    'pro_sup','wrist_flex','wrist_dev'};

RefData_str_tit = {'pelvis-tilt','pelvis-list','pelvis-rotation',...
    'pelvis-tx','pelvis-ty','pelvis-tz','hip-flexion','hip-adduction',...
    'hip-rotation','knee-angle','knee-vv','knee-rot','ankle-angle',...
    'subtalar-angle','mtp-angle','lumbar-extension','lumbar-bending',...
    'lumbar-rotation','arm-flex','arm-add','arm-rot','elbow-flex',...
    'pro-sup','wrist-flex','wrist-dev'};

idx_q = [1,5,7,9,10];
idx_refq = [1,7,10,13,16];

% col = hsv(length(ww));
figure()

for i = 1:length(idx_q)
    subplot(2,3,i)
    p = gobjects(1,length(ww));
    for k = 1:length(ww)
        x = 1:(100-1)/(size(q_opt_GC_unsc_deg(ww(k)).m,1)-1):100;
        p(k) = plot(x,q_opt_GC_unsc_deg(ww(k)).m(:,idx_q(i)),...
            'linewidth',3);
        hold on;
    end
%     meanPlusSTD = KinematicRef.LR.(RefData_str{idx_refq(i)}).mean + ...
%         2*KinematicRef.LR.(RefData_str{idx_refq(i)}).std;
%     meanMinusSTD = KinematicRef.LR.(RefData_str{idx_refq(i)}).mean - ...
%         2*KinematicRef.LR.(RefData_str{idx_refq(i)}).std;          
%     % Interpolate experimental data
%     stepQ = (size(q_opt_GC_unsc_deg(ww(k)).m,1)-1)/(size(meanPlusSTD,1)-1);
%     intervalQ = 1:stepQ:size(q_opt_GC_unsc_deg(ww(k)).m,1);
%     sampleQ = 1:size(q_opt_GC_unsc_deg(ww(k)).m,1);
%     meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
%     meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
%     hold on
%     fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k');
%     alpha(.25);
    set(gca,'Fontsize',16);
%     title(RefData_str_tit{idx_refq(i)},'Fontsize',16);
    if i == 1 || i == 4
        ylabel('(°)','Fontsize',20);
    end
    if i == 3 || i == 4 || i == 5
        L = get(gca,'XLim');
        NumTicks = 3;
        set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        xlabel('Gait cycle (%)','Fontsize',20);
        set(gca, 'XTickLabel',[0 50 100])
    else
        set(gca,'XTick',[]);
    end
end 
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Joint angles');
set(sp,'Fontsize',24);    

if ~showMainPlotsOnly
% Joint kinematics: Qdots 
figure()
for i = 1:length(idx_q)
    subplot(2,3,i)
    p = gobjects(1,length(ww));
    for k = 1:length(ww)
        x = 1:(100-1)/(size(q_opt_GC_unsc_deg(ww(k)).m,1)-1):100;
        p(k) = plot(x,qdot_opt_GC_unsc_deg(ww(k)).m(:,idx_q(i)),...
            'linewidth',3);
        hold on;
    end
    set(gca,'Fontsize',16)
    title(RefData_str_tit{idx_refq(i)},'Fontsize',16);
    if i == 1 || i == 4
        ylabel('(°/s)','Fontsize',16);
    end    
    if i == 3 || i == 4 || i == 5
        L = get(gca,'XLim');
        NumTicks = 3;
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xlabel('Gait cycle (%)','Fontsize',24);
        set(gca, 'XTickLabel',[0 50 100]);
    else
        set(gca,'XTick',[]);
    end
end
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Joint kinematics: speeds');
set(sp,'Fontsize',20);   

% Joint kinematics: Qdotdots 
figure()
for i = 1:length(idx_q)
    subplot(2,3,i)
    p = gobjects(1,length(ww));
    for k = 1:length(ww)
        x = 1:(100-1)/(size(qdotdot_opt_GC_unsc_deg(ww(k)).m,1)-1):100;
        p(k) = plot(x,qdotdot_opt_GC_unsc_deg(ww(k)).m(:,idx_q(i)),...
            'linewidth',3);
        hold on;
    end
    set(gca,'Fontsize',16)
    title(RefData_str_tit{idx_refq(i)},'Fontsize',16);
    if i == 1 || i == 4
        ylabel('(°/s²)','Fontsize',16);
    end    
    if i == 3 || i == 4 || i == 5
        L = get(gca,'XLim');
        NumTicks = 3;
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xlabel('Gait cycle (%)','Fontsize',24);
        set(gca, 'XTickLabel',[0 50 100]);
    else
        set(gca,'XTick',[]);
    end
end
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Joint kinematics: accelerations');
set(sp,'Fontsize',20);
end

% Ground reaction forces
% load([pathReferenceData '\GRFRef.mat']);
figure()
x = 1:size(GRF_opt_GC_unsc(ww(k)).m,1);
grf_names = {'x','y','z'};
GRF_str = {'Progressional force','Vertical force','Lateral force'};
for i = 1:3
    subplot(1,3,i)
    p = gobjects(1,length(ww));
    for k = 1:length(ww)
        x = 1:(100-1)/(size(GRF_opt_GC_unsc(ww(k)).m,1)-1):100;
        p(k) = plot(x,GRF_opt_GC_unsc(ww(k)).m(:,i),...
            'linewidth',3);
        hold on;  
    end
    hold on;
%     meanPlusSTD = GRFRef.LR.mean.(grf_names{i})(1:100).*100 + ...
%         2*GRFRef.LR.std.(grf_names{i})(1:100).*100;    
%     meanMinusSTD = GRFRef.LR.mean.(grf_names{i})(1:100).*100 - ...
%         2*GRFRef.LR.std.(grf_names{i})(1:100).*100;   
%     % Interpolate
%     stepGRF = (size(GRF_opt_GC_unsc(ww(k)).m,1)-1)/(size(meanPlusSTD,1)-1);
%     intervalGRF = 1:stepGRF:size(GRF_opt_GC_unsc(ww(k)).m,1);
%     sampleGRF = 1:size(GRF_opt_GC_unsc(ww(k)).m,1);
%     meanPlusSTD = interp1(intervalGRF,meanPlusSTD,sampleGRF);
%     meanMinusSTD = interp1(intervalGRF,meanMinusSTD,sampleGRF);
%     hold on
%     fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k');     
%     alpha(.25);
    set(gca,'Fontsize',20);
    if i == 1
        ylabel('Body Weight (%)','Fontsize',20);
    end
    title(GRF_str{i},'Fontsize',24);
    ylim([-50,250])
    L = get(gca,'XLim');
    NumTicks = 3;
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    xlabel('Gait cycle (%)','Fontsize',20);
    set(gca, 'XTickLabel',[0 50 100])
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))       
end
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Ground reaction forces');
set(sp,'Fontsize',24);

% Joint kinetics
% load([pathReferenceData '\KineticRef.mat']);
figure()
x = 1:size(mom_opt_GC(ww(k)).m,1);
for i = 1:length(idx_q)-1
    subplot(2,2,i)
    for k = 1:length(ww)
        x = 1:(100-1)/(size(mom_opt_GC(ww(k)).m,1)-1):100;
        p(k) = plot(x,mom_opt_GC(ww(k)).m(:,idx_q(i+1)),...
            'linewidth',3);
        hold on;
    end
%     meanPlusSTD = KineticRef.LR.(RefData_str{idx_refq(i+1)}).mean + ...
%         2*KineticRef.LR.(RefData_str{idx_refq(i+1)}).std;
%     meanMinusSTD = KineticRef.LR.(RefData_str{idx_refq(i+1)}).mean - ...
%         2*KineticRef.LR.(RefData_str{idx_refq(i+1)}).std;
%     % Interpolate
%     stepID = (size(mom_opt_GC(ww(k)).m,1)-1)/(size(meanPlusSTD,1)-1);
%     intervalID = 1:stepID:size(mom_opt_GC(ww(k)).m,1);
%     sampleID = 1:size(mom_opt_GC(ww(k)).m,1);
%     meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
%     meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID); 
%     hold on
%     fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k');
%     alpha(.25);
    set(gca,'Fontsize',16);
    if i == 1 || i == 4
        ylabel('(Nm/kg)','Fontsize',16);
    end    
    if i == 3 || i == 4 || i == 5
        L = get(gca,'XLim');
        NumTicks = 3;
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xlabel('Gait cycle (%)','Fontsize',24);
        set(gca, 'XTickLabel',[0 50 100]);
    else
        set(gca,'XTick',[]);
    end
    title(RefData_str_tit{idx_refq(i+1)},'Fontsize',16);
    ylim([-2 2])
    NumTicks = 3;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks)) 
end      
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Joint kinetics');
set(sp,'Fontsize',20);

% Muscle activations
muscleNames_str = {'hamstrings_r','bifemsh_r','glut-max_r',...
    'iliopsoas_r','rect-fem_r','vasti_r','gastroc_r','soleus_r',...
    'tib-ant_r'}; 

% load([pathReferenceData,'\EMGrefdata.mat']);

EMGchannel = [4,4,10,99,1,2,9,8,7];
EMGcol = [1,1,1,0,1,1,1,1,1];

x = 1:size(a_opt_GC(ww(k)).m,1);
figure()
for i = 1:size(a_opt_GC(ww(k)).m,2)/2
    subplot(3,3,i)
    p = gobjects(1,length(ww));
    NMuscle = size(a_opt_GC(ww(k)).m,2);
    for k = 1:length(ww)
        x = 1:(100-1)/(size(a_opt_GC(ww(k)).m,1)-1):100;
        p(k) = plot(x,a_opt_GC(ww(k)).m(:,i+NMuscle/2),...
            'linewidth',3);
        hold on;
    end
    title(muscleNames_str{i}(1:end-2),'Fontsize',16);
    set(gca,'Fontsize',16)
    if i == 1 || i == 4 || i == 7
        ylabel('(-)','Fontsize',20);
    end
    if i == 7 || i == 8 || i == 9
        L = get(gca,'XLim');
        NumTicks = 3;
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xlabel('Gait cycle (%)','Fontsize',20);
        set(gca, 'XTickLabel',[0 50 100]);
    else
        set(gca,'XTick',[]);
    end
%     if EMGcol(i)
%         meanPlusSTD = EMGrefdata.mean(:,EMGchannel(i)) + 2*EMGrefdata.std(:,EMGchannel(i));
%         meanMinusSTD = EMGrefdata.mean(:,EMGchannel(i)) - 2*EMGrefdata.std(:,EMGchannel(i));
%         % Interpolate
%         stepa = (size(a_opt_GC(ww(k)).m,1)-1)/(size(meanMinusSTD,1)-1);
%         intervala = 1:stepa:size(a_opt_GC(ww(k)).m,1);
%         samplea = 1:size(a_opt_GC(ww(k)).m,1);
%         meanPlusSTD = interp1(intervala,meanPlusSTD,samplea);
%         meanMinusSTD = interp1(intervala,meanMinusSTD,samplea);     
%         hold on
%         fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k');
%         alpha(.25);            
%     end
    ylim([0,0.7]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
end
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Muscle activations');
set(sp,'Fontsize',24);
 
% Torque actuator (lumbar extension)
% figure()    
% p = gobjects(1,length(ww));
% for k = 1:length(ww)
%     x = 1:(100-1)/(size(tau_opt_GC_unsc(ww(k)).m,1)-1):100;
%     p(k) = plot(x,tau_opt_GC_unsc(ww(k)).m(:,1),'color',...
%         col.(der),'linewidth',3);
%     hold on;
% end
% ylabel('Moment (Nm)','Fontsize',16);
% L = get(gca,'XLim');
% NumTicks = 3;
% set(gca,'XTick',linspace(L(1),L(2),NumTicks));
% xlabel('Gait cycle (%)','Fontsize',24);
% set(gca, 'XTickLabel',[0 50 100]);
% title('Torque actuator: lumbar extension','Fontsize',20);
% set(gca,'Fontsize',16)
% if showlegend
%     l = legend(p,legend_case);
%     set(l,'Fontsize',16)
% end

% Muscle-tendon forces
x = 1:size(FTtilde_opt_GC_unsc(ww(k)).m,1);
figure()
for i = 1:size(FTtilde_opt_GC_unsc(ww(k)).m,2)/2
    subplot(3,3,i)
    p = gobjects(1,length(ww));
    NMuscle = size(FTtilde_opt_GC_unsc(ww(k)).m,2);
    for k = 1:length(ww)
        x = 1:(100-1)/(size(FTtilde_opt_GC_unsc(ww(k)).m,1)-1):100;
        p(k) = plot(x,FTtilde_opt_GC_unsc(ww(k)).m(:,i+NMuscle/2),...
            'linewidth',3);
        hold on;
    end
    title(muscleNames_str{i}(1:end-2),'Fontsize',16);
    set(gca,'Fontsize',16)
    if i == 1 || i == 4 || i == 7
        ylabel('(-)','Fontsize',20);
    end
    if i == 7 || i == 8 || i == 9
        L = get(gca,'XLim');
        NumTicks = 3;
        set(gca,'XTick',linspace(L(1),L(2),NumTicks));
        xlabel('Gait cycle (%)','Fontsize',20);
        set(gca, 'XTickLabel',[0 50 100]);
    else
        set(gca,'XTick',[]);
    end
    
    ylim([0,0.7]);
    NumTicks = 2;
    LY = get(gca,'YLim');
    set(gca,'YTick',linspace(LY(1),LY(2),NumTicks))
end
if showlegend
    l = legend(p,legend_case);
    set(l,'Fontsize',16)
end
sp = suptitle('Normalized tendon forces');
set(sp,'Fontsize',24);
