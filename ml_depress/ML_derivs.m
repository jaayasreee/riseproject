function dydt = ML_derivs(t, y, I, coefs)
% function [dydt] = derivsI(t, y, I)
% The derivs function
persistent V;
persistent n;
persistent Ser;
persistent Sef;
V = y(1,:);
n = y(2,:);
Ser = y(3,:);
Sef = y(4,:);

%dXidt = minf(V)*(1-minf(V)); %MaxRand*sqrt(t-TLast)*randn;
%totalSynapticDrive = gEPSP*((Sef - Ser).*(E_EPSP-V)) + gIPSP*((Sif - Sir).*(E_IPSP-V));

dydt = [1.0./coefs.C*(I-coefs.gL.*(V-coefs.EL) - ...
    coefs.gNa.*(1.0./(1.0 + exp((coefs.Vhalfm - V)/coefs.km ))).*(V-coefs.ENa) - ...
    coefs.gK.*n.*(V-coefs.EK) + ...
    coefs.gEPSP.*(Sef-Ser).*(coefs.E_EPSP-V));...%(coefs.E_EPSP-V))
    (1.0./(1.0+exp((coefs.Vhalfn - V)/coefs.kn ))-n)./exp(-0.07*V-3);...
    -Ser/coefs.tauEPSPr;
    -Sef/coefs.tauEPSPf];

