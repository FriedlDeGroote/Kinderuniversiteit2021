% This function generates polynomials to approximate muscle-tendon lengths
% and moment arms. The code is from Wouter Aerts and is adapted to be used
% with CasADi and external function

function [muscle_spanning_joint_INFO,MuscleInfo] =  PolynomialFit(MuscleData)

    %% Construct the polynomials for the moment arms and muscle length

    muscle_sel=[];
    
    for m_nr = 1:length(MuscleData.muscle_names)
        if strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_r') || strcmp(MuscleData.muscle_names{m_nr}(end-1:end), '_l') % was _l before
            muscle_sel = [muscle_sel m_nr];
        end
    end
    
    muscle_spanning_joint_INFO = squeeze(sum(MuscleData.dM, 1));
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO<=0.0001 & muscle_spanning_joint_INFO>=-0.0001) = 0;
    muscle_spanning_joint_INFO(muscle_spanning_joint_INFO~=0) = 1;
      
    q_all = MuscleData.q;
    
    max_order = 9;
    threshold = 0.003; % 3mm
    nr_samples = length(q_all(:,1));
    
    lMT_all_error = zeros(length(muscle_sel), 1);
    DM_all_error = zeros(length(muscle_sel), length(q_all(1,:)));
    order_all = zeros(length(muscle_sel), 1);
    
    for m_nr=1:length(muscle_sel)
        muscle_index = muscle_sel(m_nr);
        
        index_dof_crossing = find(muscle_spanning_joint_INFO(muscle_index,:)==1);
        nr_dof_crossing = length(index_dof_crossing);
        
        lMT = MuscleData.lMT(:,muscle_index);
        dM = zeros(nr_samples, nr_dof_crossing);
        dM_recon = dM;
        for dof_nr = 1:nr_dof_crossing
            dM(:,dof_nr) = MuscleData.dM(:,muscle_index,index_dof_crossing(dof_nr));
        end
        
        criterion_full_filled = 0;
        order = 3;
        while criterion_full_filled==0
            [mat,diff_mat_q] = n_art_mat_3(q_all(:,index_dof_crossing), order);
            nr_coeffs = length(mat(1,:));
            
            diff_mat_q_all = zeros(nr_samples*nr_dof_crossing, nr_coeffs);
            for dof_nr = 1:nr_dof_crossing
                diff_mat_q_all(nr_samples*(dof_nr-1)+1:nr_samples*dof_nr,:) = -squeeze(diff_mat_q(:,:,dof_nr));
            end
            temp = [mat ; diff_mat_q_all];
            temp2 = [lMT; dM(:)];
            coeff=[mat ; diff_mat_q_all]\[lMT; dM(:)];
            dM_recon = zeros(nr_samples, nr_dof_crossing);
            for dof_nr = 1:nr_dof_crossing
                dM_recon(:,dof_nr) = (-squeeze(diff_mat_q(:,:,dof_nr)))*coeff;
            end
            lMT_recon=mat*coeff;
            
            lMT_error_rms = sqrt(mean((lMT - lMT_recon).^2));
            dm_error_rms = sqrt(mean((dM - dM_recon).^2));
            
            criterion_full_filled = lMT_error_rms<=threshold & max(dm_error_rms)<=threshold;
            if order==max_order
                criterion_full_filled = 1;
            end
            if criterion_full_filled==0
                order = order+1;
            end
        end
        
        MuscleInfo.muscle(m_nr).DOF = MuscleData.dof_names(index_dof_crossing);
        MuscleInfo.muscle(m_nr).m_name = MuscleData.muscle_names{muscle_index};
        MuscleInfo.muscle(m_nr).coeff = coeff;
        MuscleInfo.muscle(m_nr).order = order;
        MuscleInfo.muscle(m_nr).lMT_error_rms = lMT_error_rms;
        MuscleInfo.muscle(m_nr).dm_error_rms = dm_error_rms;
        
        lMT_all_error(m_nr) = lMT_error_rms;
        DM_all_error(m_nr, index_dof_crossing) = dm_error_rms;
        order_all(m_nr) = order;
        
        if m_nr==1023 %10
            dof_nr_index = 1;
            figure;
            hold on;
            plot(q_all(:,index_dof_crossing(dof_nr_index)), dM(:,dof_nr_index), '*b')
            plot(q_all(:,index_dof_crossing(dof_nr_index)), dM_recon(:,dof_nr_index), '*r')
            
            figure;
            hold on;
            plot(q_all(:,index_dof_crossing(dof_nr_index)), lMT, '*b')
            plot(q_all(:,index_dof_crossing(dof_nr_index)), lMT_recon, '*r')
        end       
    end
    
    figure();
%     subplot(2,1,1)
    hold on;
    plot(lMT_all_error)
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold, threshold], 'r', 'linewidth', 2)
    suptitle('RMS error on the approximated muscle-tendon length')
    ylabel('RMS error (m)')
%     set(gca, 'XTickLabel', [])
%     xticklabel_rotate(1:length(muscle_sel),90,MuscleData.muscle_names(muscle_sel))
    
    figure();
%     subplot(2,1,2)
    hold on;
    plot(max(DM_all_error, [], 2))
    xlimits = get(gca, 'XLim');
    plot(xlimits, [threshold, threshold], 'r', 'linewidth', 2)
    suptitle('maximal RMS error on the approximated muscle moment arm')
    ylabel('RMS error (m)')
%     xticklabel_rotate(1:length(muscle_sel),90,MuscleData.muscle_names(muscle_sel))
%     title('maximal RMS error on the approximated muscle moment arm')

    figure();
    hold on;
    plot(order_all)
    ylim([0 max_order+1])
    xlimits = get(gca, 'XLim');
    plot(xlimits, [max_order, max_order], 'r', 'linewidth', 2)
    suptitle('Order of the polynomial approximation')
    ylabel('Order')
%     xticklabel_rotate(1:length(muscle_sel),90,MuscleData.muscle_names(muscle_sel))
    
%     save('MuscleInfo_3D','MuscleInfo');
% end


% if check_implementation
%     % do a muscle analysis
%     dummy_motion = importdata('3D_results/dummy_motion.mot');
%     time = dummy_motion.data(:,1);
%     output_path=fullfile(pwd,'3D_results', 'MuscleAnalysis_poly');mkdir(output_path);
% %     OpenSim_Muscle_Analysis(fullfile(pwd, '3D_results','dummy_motion.mot'),[model_file_name(1:end-5), '_poly.osim'],output_path,[time(1) time(end)])
% 
% %     data_folder = 'C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\MuscleSmoothing\MuscleLengthPolynomials\3D_results';
% %     data_KS = importdata('C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\Data_overground_3D\KS\KS_motion_07.mot');
% %     dof_names = data_KS.colheaders(2:end);
%     
% %     left_leg_DOFs = find(~cellfun('isempty', strfind(dof_names,'_l')));
% %     pelvis_DOFs = find(~cellfun('isempty', strfind(dof_names,'pelvis')));
% %     lumbar_DOFs = find(~cellfun('isempty', strfind(dof_names,'lumbar')));
% %     all_DOFS = [left_leg_DOFs, lumbar_DOFs];
% %     all_DOFS = unique(all_DOFS);
% %     all_DOFS = setdiff(all_DOFS, pelvis_DOFs); % exclude the pelvis DOFS
% %     dof_names = dof_names(all_DOFS);
%     
% %     for dof_nr = 1:length(dof_names)
% %         dM=importdata(fullfile(data_folder,['MuscleAnalysis\dummy_motion_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
% %         MuscleData_normal.dM(:,:,dof_nr) = dM.data(:,2:end);
% %     end
% %     Length=importdata(fullfile(data_folder,'MuscleAnalysis\dummy_motion_MuscleAnalysis_Length.sto'));
% %     MuscleData_normal.lMT = Length.data(:,2:end);
% %     
% %     for dof_nr = 1:length(dof_names)
% %         dM=importdata(fullfile(data_folder,['MuscleAnalysis_poly\dummy_motion_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
% %         MuscleData_poly.dM(:,:,dof_nr) = dM.data(:,2:end);
% %     end
% %     Length=importdata(fullfile(data_folder,'MuscleAnalysis_poly\dummy_motion_MuscleAnalysis_Length.sto'));
% %     MuscleData_poly.lMT = Length.data(:,2:end);
%     
%     data_folder = 'C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\Data_overground_3D\TestMuscleAnalysisPoly\';
%     data_KS = importdata('C:\Users\u0085720\Documents\Doctoraat\Projecten\MultipleShooting\Data_overground_3D\KS\KS_motion_07.mot');
%     dof_names = data_KS.colheaders(2:end);
%     
%     for dof_nr = 1:length(dof_names)
%         dM=importdata(fullfile(data_folder,['Results_normal\p2-scaled_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
%         MuscleData_normal.dM(:,:,dof_nr) = dM.data(:,2:end);
%     end
%     Length=importdata(fullfile(data_folder,'Results_normal\p2-scaled_MuscleAnalysis_Length.sto'));
%     MuscleData_normal.lMT = Length.data(:,2:end);
%     
%     for dof_nr = 1:length(dof_names)
%         dM=importdata(fullfile(data_folder,['Results_poly\p2-scaled_MuscleAnalysis_MomentArm_', dof_names{dof_nr}, '.sto']));
%         MuscleData_poly.dM(:,:,dof_nr) = dM.data(:,2:end);
%     end
%     Length=importdata(fullfile(data_folder,'Results_poly\p2-scaled_MuscleAnalysis_Length.sto'));
%     MuscleData_poly.lMT = Length.data(:,2:end);
% 
%     rms_lMT = sqrt(mean((MuscleData_poly.lMT - MuscleData_normal.lMT).^2));
%     rms_dM = squeeze(sqrt(mean((MuscleData_poly.dM-MuscleData_normal.dM).^2)));
%     
%     figure;
%     plot(rms_lMT)
%     suptitle('RMS error on the approximated muscle-tendon length')
%     ylabel('RMS error (m)')
%     xticklabel_rotate(1:length(Length.data(1,2:end)),90,Length.colheaders(2:end))
%     
%     figure;
%     plot(max(rms_dM, [], 2))
%     suptitle('maximal RMS error on the approximated muscle moment arm')
%     ylabel('RMS error (m)')
%     xticklabel_rotate(1:length(Length.data(1,2:end)),90,Length.colheaders(2:end))
%     
%     right_muscles = find(~cellfun('isempty', strfind(Length.colheaders(2:end),'_r')));
%     left_muscles = setdiff(1:92, right_muscles);
%     
%     figure;
%     plot(rms_lMT(right_muscles))
%     suptitle('Right')
%     
%     figure;
%     plot(rms_lMT(left_muscles))
%     suptitle('Left')
%     
% end
