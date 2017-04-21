
close all
clc

% load OpenFOAM data
temp3        = importdata('~/thesis/beskowFiles/chan_Edmond/postProcessing/patchExpression_uTau/0/bottomWall.csv'); 
temp2 = struct2cell(temp3);
temp1 = cell2mat(temp2(1));

bottomWall(:, 1)  = temp1(:, 1);
bottomWall(:, 2)  = temp1(:, 2);

% data location folder
filename = '~/thesis/beskowFiles/chan_Edmond/Ux';

% total number of time steps
maxTime  = 28;

% specify viscosity value
nu       = 1.506e-5;

% specify 1st cell's distance from wall
dy       = 8.214444444e-05;

% initiate Ux average vector per time
Uavg     = zeros(maxTime, 1);

% extract data from files
for t=0:maxTime-1
    
    % reset Ux
    Ux = []; 
    
    % update file    
    for i=0:6
        % extract 2d data
        fileUpdate = strcat(filename, int2str(i), '.');
        thisFile  = strcat(fileUpdate, int2str(t),'.csv'); 
        
        % check if file is empty
        s = dir(thisFile);
        if s.bytes >1
            temp = csvread(thisFile, 1, 0);
            Ux    = [temp(:, 1); Ux];
        else
            disp(['t = ', int2str(t), ' i = ', int2str(i)]);
        end
    end

    % calculate the mean
    Uavg(t+1) = mean(Ux);
end

% calculate gradient of Ux
gradUx = ( Uavg - 0 )/dy;

% calculate the friction velocity
uTau = sqrt(nu*gradUx);

% plot uTau across time (per second)
plot(uTau)

grid minor
xlabel('time (sec)')
ylabel('u_{tau} (m/sec)')
title('friction velocity as a function of time')

% plot data from OpenFOAM
hold on
dN      = 20e3; % per second sampling
t_OF    = bottomWall(1:dN:end-dN, 1)+1;
uTau_OF = bottomWall(1:dN:end-dN, 2);
plot(t_OF, uTau_OF, 'ro-')



% testing stuff
corrFactor = uTau_OF(2) - uTau(2);

plot([uTau(1); uTau(2:end) + corrFactor], 'kx')


% legend
legend('matlab', 'OpenFOAM', 'correction factor')

% font size
set(gca,'fontsize',14)
