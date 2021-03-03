function params = ExtractModelParameters(modelPath,muscleNames)

import org.opensim.modeling.*

model = Model(modelPath);

nom = length(muscleNames);
params = zeros(5,nom);

muscles = model.getMuscles();

for i = 1:nom
   muscle = muscles.get(muscleNames{i});
   params(1,i) = muscle.getMaxIsometricForce();
   params(2,i) = muscle.getOptimalFiberLength();
   params(3,i) = muscle.getTendonSlackLength();
   params(4,i) = muscle.getPennationAngleAtOptimalFiberLength();
   params(5,i) = params(2,i).*10;
end

end
