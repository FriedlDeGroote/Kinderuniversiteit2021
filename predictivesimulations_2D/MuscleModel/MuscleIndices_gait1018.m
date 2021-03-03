% return muscle indices as compared to muscles for gait2392

function musi = MuscleIndices_gait1018(muscleNames)
   
muscleNames_all = {'hamstrings_r','bifemsh_r','glut_max_r',...
    'iliopsoas_r','rect_fem_r','vasti_r','gastroc_r','soleus_r',...
    'tib_ant_r'};
    
count = 1;
musi = zeros(1,length(muscleNames));
for i = 1:length(muscleNames)       
    if (find(strcmp(muscleNames_all,muscleNames{i})) ~= 0)        
        musi(count) = find(strcmp(muscleNames_all,muscleNames{i}));
        count = count + 1;
    end
end

end