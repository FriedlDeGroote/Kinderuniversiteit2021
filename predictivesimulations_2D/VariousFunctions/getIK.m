% This function returns the IK given as inputs the path to the
% file containing the output of IK and the joints of
% interest in the desired order

function Qs = getIK(pathIK,joints)

Qsall = importdata(pathIK);
Qs.time = Qsall.data(:,strcmp(Qsall.colheaders,{'time'}));
Qs.all(:,1) = Qs.time;
Qs.colheaders{1,1} = 'time';
count = 1;
for i = 1:length(joints)
    count = count + 1;
    if strcmp(joints{i},'pelvis_tx') || strcmp(joints{i},'pelvis_ty') || strcmp(joints{i},'pelvis_tz')
        Qs.(joints{i}) = Qsall.data(:,strcmp(Qsall.colheaders,joints{i}));
    else
        Qs.(joints{i}) = Qsall.data(:,strcmp(Qsall.colheaders,joints{i})).*(pi/180);
    end
    Qs.all(:,count) = Qs.(joints{i});
    Qs.colheaders(1,count) = Qsall.colheaders(1,strcmp(Qsall.colheaders,joints{i}));
end

% filter kinematics
order = 4;
cutoff_low = 6;
fs=1/mean(diff(Qs.all(:,1)));
[af,bf] = butter(order/2,cutoff_low./(0.5*fs),'low');
Qs.allfilt = Qs.all;
Qs.allfilt(:,2:end) = filtfilt(af,bf,Qs.allfilt(:,2:end));