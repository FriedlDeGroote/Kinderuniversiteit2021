function dadt = BackActivationDynamics(e,a)
tau = 0.035;
dadt = (e-a)./tau;