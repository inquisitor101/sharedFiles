% start from scratch
clear all
close all
clc

% fixed parameters 
nu_LES    = 1.506e-5;
delta_LES = 0.0134;
Ub_LES    = 6.7123;

% base folder directory
baseFile = '~/thesis/beskowFiles/channelFlow/multiChan';

% available simulation files
allFiles = {'chan610', 'chan610_crc1',...
            'chan610_nz60', 'chan610_nz90',...
            'chan610_y2_nz60', 'chan610r', 'eddsTry'};

% averaging times (uTau)
avgTimes = {'30', '100', '14.5', '14.5',...
            '21', '32.5', '1000'};

% averaged profiles (<uu>)
avgProfs = {'30', '115', '29', '22',...
            '39.5', '32.5', '1243.18'};

% Reynolds stresses
ReStress = {'UPrime2Mean_XX.xy', 'UPrime2Mean_YY.xy',...
            'UPrime2Mean_ZZ.xy', 'UPrime2Mean_XY.xy'};

% load DNS data
fileUpdate = strcat(baseFile, '/DNS_chan300.mat');
load(fileUpdate);


% loop across all 6 simulations
for i=1:7
   
    % current simulation file
    simFile = strcat(baseFile, '/', allFiles{i}, '/postProcessing');
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % get uTau
    %
    % surf until the averaging time
    thisFile = strcat(simFile, '/patchExpression_uTau/');
    thisFile1 = strcat(thisFile, avgTimes{i}, '/bottomWall');
    % get uTau_LES bottom
    temp = importdata(thisFile1, ' ', 1);
    temp = temp.data;
    uTau_vect = temp(:, 2);
    uTau_LES_1 = mean(uTau_vect);
    
    thisFile2 = strcat(thisFile, avgTimes{i}, '/topWall');
    % get uTau_LES top
    temp = importdata(thisFile2, ' ', 1);
    temp = temp.data;
    uTau_vect = temp(:, 2);
    uTau_LES_2 = mean(uTau_vect);
    
    % average overall uTau
    uTau_LES = ( uTau_LES_1 + uTau_LES_2 )/2;
    %
    % get uTau
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % get Reynolds stresses
    %
    % enter collapsedFields directory
    thisFile = strcat(simFile, '/collapsedFields/', avgProfs{i}, '/');

    figure(i)
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    % loop across 4 Reynolds stresses
    for j=1:4
        
        updateFile = strcat(thisFile, ReStress{j});
    
        temp = importdata(updateFile, ' ', 1);
        temp = temp.data;
        
        y_LES = temp(1:end/2, 1);
        
        % eddsTry data ! (simulation #7)
        if i==7
            delta_LES = 1;
            nu_LES    = 1.6526e-4;
            Ub_LES    = 1;
        end
        
        ydelta_LES = y_LES/delta_LES;
        yPlus_LES  = y_LES*uTau_LES/nu_LES;
        
        switch j
            case 1
                uu_LES = temp(1:end/2, 2);
            case 2
                vv_LES = temp(1:end/2, 2);
            case 3
                ww_LES = temp(1:end/2, 2);
            case 4
                uv_LES = temp(1:end/2, 2);
            otherwise
                disp('Unknown input in switch statement');
        end
      

    end
    %
    % get Reynolds stresses
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % plot first subplot data using inner scaling
    subplot(211)

    semilogx(yPlus, u_rms.^2)
    hold on
    semilogx(yPlus_LES, uu_LES/uTau_LES^2, 'b--')

    semilogx(yPlus, v_rms.^2, 'r')
    semilogx(yPlus_LES, vv_LES/uTau_LES^2, 'r--')

    semilogx(yPlus, w_rms.^2, 'g')
    semilogx(yPlus_LES, ww_LES/uTau_LES^2, 'g--')

    semilogx(yPlus, uv, 'm')
    semilogx(yPlus_LES, -uv_LES/uTau_LES^2, 'm--')

    xlabel('yPlus')
    ylabel('<u_i u_j>/uTau^2')
    switch i
        case 1
            title('DNS vs LES for chan610')
        case 2
            title('DNS vs LES for chan610\_crc1')
        case 3
            title('DNS vs LES for chan610\_nz60')
        case 4
            title('DNS vs LES for chan610\_nz90')
        case 5
            title('DNS vs LES for chan610\_y2\_nz60')
        case 6
            title('DNS vs LES for chan610r')
        case 7
            title('DNS vs LES for eddsTry')
        otherwise
            disp('Unknown input in switch statement');
    end
    legend('uu DNS', 'uu LES', 'vv DNS', 'vv LES', 'ww DNS', 'ww LES', '-uv DNS', '-uv LES')
    set(gca,'fontsize',14)
    grid on 
    grid minor

    % plot second subplot data using outer scaling

    % re-scaling into outer scales coefficient using:
    %
    % ReTau / Re_m = #  ==> uTau/(2Ub) = C
    % u_rms.^2 * C^2 = <uu>/Ub^2, outer scaling
    ReTau_DNS = 297.899;
    Re_m_DNS  = 10039.1;
    coef      = 2*ReTau_DNS/Re_m_DNS;

    subplot(212)
    plot(ydelta, (u_rms.^2)*coef^2)
    hold on
    plot(ydelta_LES, uu_LES/Ub_LES^2, 'b--')

    plot(ydelta, (v_rms.^2)*coef^2, 'r')
    plot(ydelta_LES, vv_LES/Ub_LES^2, 'r--')

    plot(ydelta, (w_rms.^2)*coef^2, 'g')
    plot(ydelta_LES, ww_LES/Ub_LES^2, 'g--')

    plot(ydelta, uv*coef^2, 'm')
    plot(ydelta_LES, -uv_LES/Ub_LES^2, 'm--')

    xlabel('y/delta')
    ylabel('<u_i u_j>/Ub^2')

    legend('uu DNS', 'uu LES', 'vv DNS', 'vv LES', 'ww DNS', 'ww LES', '-uv DNS', '-uv LES')
    set(gca,'fontsize',14)
    grid on 
    grid minor

    % u-velocity gradient

    updateFile = strcat(thisFile, 'UMean_X.xy');
    temp = importdata(updateFile, ' ', 1);
    temp = temp.data;

    figure(8)

    uMean_LES  = temp(1:end/2, 2);
    y_grad_LES = temp(1:end/2, 1); 
    yPlus_grad_LES = y_grad_LES*uTau_LES/nu_LES;

    semilogx(yPlus_grad_LES, uMean_LES/uTau_LES)
    hold on 
end



figure(8)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
logLaw = (1/0.38)*log(yPlus) + 4.17;
semilogx(yPlus, u_mean)
semilogx(yPlus, logLaw, 'k')
xlabel('yPlus')
ylabel('uPlus')
legend('chan610', 'chan610\_crc1',...
       'chan610\_nz60', 'chan610\_nz90',...
       'chan610\_y2\_nz60', 'chan610r',...
       'eddsTry', 'DNS', 'log-law');
xlim([1, max(yPlus)])
set(gca,'fontsize',16)
grid on
grid minor
