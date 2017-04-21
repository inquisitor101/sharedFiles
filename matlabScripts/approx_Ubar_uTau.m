function [Ubar, uTau] = approx_Ubar_uTau(isNormalized, ReTau, delta, nu, plotMe)
% % % %
% this is a script for calculating the average 
% velocity given a friction Reynolds number 
% Re_tau, using the relation in Pope's chp 7
% and section 'channel flow.' Commonly referred
% to as Dean's correlation.
%
% INPUT 
%       isNormalized    find nu with all else 1
%       ReTau           friction Reynolds number
%       delta           channel half-width
%       nu              viscosity
%       plotMe          1 (plot solution) 0 (don't) 
% OUTPUT
%       Ubar            bulk velocity 
%       uTau            friction velocity
% 
% This relation is approximated from figure 7.11
% and its expression is:
%                                                
% Re_tau = 0.09 Re^0.88                
%    where,                                       
%           Re     = (2*delta*Ubar) / nu          %
%           Re_tau = (  delta*uTau) / nu          %
%           uTau   = sqrt(tau_Wall/rho)           %
% % % % % % % % % % % % % % % % % % % % % % % % % %

if isNormalized

    temp = (ReTau/0.09)^(1/0.88);   
    Ubar = 2*delta/temp;
    uTau = 1;
    
    disp(' ');
    disp('normalized mode is ON');
else
    
    temp = (ReTau/0.09)^(1/0.88);
    Ubar = (temp*nu)/(2*delta);
    uTau = (ReTau*nu)/delta;
    
    disp(' ');
    disp('normalized mode is OFF');
end

% diagram plot
if plotMe
    Re_all    = 2e3:10:1e6;
    ReTau_all = 0.09*Re_all.^0.88;

    loglog(Re_all, ReTau_all)
    
    hold on 
    Re = (2*delta*Ubar)/nu;
    plot(Re, ReTau, 'ro')
    legend('ReTau','current ReTau')
    grid on 
    xlabel('Re')
    ylabel('ReTau')
    title('Approximate Relation Between Re and ReTau for Channel Flow')
end