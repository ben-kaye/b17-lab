%% gait.m -- Starter File for Gait Data Analysis
%% B17 Biomechanics -- Hilary Term 2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

% Load data into arrays Xmk (markers) and Xfp (forceplate) for Subject D
mk = csvread('Subject-D-markers.csv');
fp = csvread('Subject-D-forceplate.csv');

% Find the number of rows and columns in the input files
[rmk,cmk] = size(mk);    [rfp,cfp] = size(fp);

% initialize arrays of datapt number in marker file & assign to an array
% note that ':' means 'all the rows (or columns)'
datapt = zeros(rmk,1);
datapt(:,1) = mk(:,1);

% Initialize arrays for marker data (right leg only)
sac = zeros(rmk,3);      % sacrum (SACR)
asi = zeros(rmk,3);      % anterior superior iliac spine (ASIS)
thi = zeros(rmk,3);      % thigh wand marker (THI) 
kne = zeros(rmk,3);      % lateral femoral condyle (KNE)
tib = zeros(rmk,3);      % tibia wand marker (TIB)
ank = zeros(rmk,3);      % later malleolus (LMA)
hee = zeros(rmk,3);      % heel (HEE)
toe = zeros(rmk,3);      % 2nd metatarsal head (TOE)

% Assign xyz coordinates of markers to right side sacrum, asis, thigh, knee,
% ankle, heel, and toe arrays
sac(:,1:3) = mk(:,2:4);
asi(:,1:3) = mk(:,8:10);             
thi(:,1:3) = mk(:,29:31);      
kne(:,1:3) = mk(:,32:34); 
tib(:,1:3) = mk(:,35:37);
ank(:,1:3) = mk(:,38:40);      
hee(:,1:3) = mk(:,41:43);      
toe(:,1:3) = mk(:,44:46);    

%%%%%%%%%%%%% YOU NEED TO CONTINUE THE CODE FROM HERE 

% Plot yz trajectories
     figure(1)
     plot(sac(:,2),sac(:,3))
     hold on
     text(sac(rmk,2),sac(rmk,3),'SACRUM')
     
     plot(asi(:,2),asi(:,3))
     text(asi(rmk,2),asi(rmk,3),'ILIAC')
     
% continue with rest of markers     
    plot(thi(:,2), thi(:,3))
    text(thi(rmk,2), thi(rmk,3), 'THIGH')
    
    plot(kne(:,2), kne(:,3))
    text(kne(rmk, 2),kne(rmk,3), 'KNEE')
    
    plot(tib(:,2), tib(:,3))
    text(tib(rmk, 2),tib(rmk,3),'TIBIA')
    
    plot(ank(:,2), ank(:,3))
    text(ank(rmk, 2), ank(rmk,3), 'ANKLE')
    
    plot(hee(:,2), hee(:,3))
    text(hee(rmk, 2), hee(rmk,3), 'HEEL')
    
    plot(toe(:,2), toe(:,3))
    text(toe(rmk, 2), toe(rmk,3), 'TOE')
    hold off
    
     xlabel('Y (Posterior - Anterior)')
     ylabel('Z (Inferior - Superior)')
     title('Subject D Marker Trajectories')
     axis equal
     
% compute com velocity (by approximating the com to the sacrum)
sample_rate = 100; % [Hz]
vel_com = (sac(2:end, :) - sac(1:end-1, :))*sample_rate; % first order approximation to velocity


time = (datapt - datapt(1))/sample_rate; % calculate time from the marker index


figure(2)
plot(time(1:end-1), vel_com(:,2))
avg_forvel = mean(vel_com(:,2));
hold on
yline(avg_forvel, '--')
hold off
legend('COM', 'avg')
xlabel('time (s)')
ylabel('forward velocity (mms^{-1})')

% compute shank lengths in the YZ plane and in cartesian coordinates
shank_vec = ank - kne;
shank_length_3d = vecnorm(shank_vec');
shank_length_yz = vecnorm(shank_vec(:,2:3)');
figure(3)
plot(time, shank_length_3d) 
hold on
plot(time, shank_length_yz)
legend('3D', 'YZ');
ylabel('length (mm)')
xlabel('time (s)')
   hold off
   
% plot foot section vertical trajcetory
figure(4)


% indices of events
FF = 384; % foot flat
HS = 364; % heel strike
HO = 404; % heel off
TO = 426; % retrospectively included

gait_pos = [ FF, HS, HO, TO ];
gait_idx = gait_pos - datapt(1) + 1;
% get indices


plot(datapt, ank(:,3));
hold on
plot(datapt, hee(:,3), '-o', 'MarkerIndices', gait_idx([2,3]));
plot(datapt, toe(:,3), '-o', 'MarkerIndices', gait_idx([1,4]));
text(FF, toe(gait_idx(1),3), 'FF')
text(HS, hee(gait_idx(2),3), 'HS')  
text(HO, hee(gait_idx(3),3), 'HO')
text(TO, toe(gait_idx(4),3), 'TO')
hold off
legend('Ankle', 'Heel', 'Toe')
xlabel('Marker index')
ylabel('Z position (mm)')

% force plate data processing
datafp = fp(:,1);
Fz = fp(:,4);

figure(5)
plot(datafp, Fz)
xlabel('fp data point')
ylabel('Vertical reaction force (N)')

% calculate knee angle

figure(6)
thigh_vec_yz = kne(:, 2:3) - asi(:, 2:3);
thigh_length_yz = vecnorm(thigh_vec_yz');
knee_angle = acos(sum(thigh_vec_yz.*shank_vec(:,2:3),2)'./thigh_length_yz./shank_length_yz)*180/pi; % use dot product angle equation
plot(datapt, knee_angle);
xlabel('marker index')
ylabel('Knee angle, (^\circ) +ve flexion');

% calculate stride length
HS1 = gait_idx(1);
HS2 = 458 - datapt(1) + 1; % get indices of strikes
stride_length = hee(HS2, 2) - hee(HS1, 2) % compute heel travel
cadence = 60*2*100/(HS2 - HS1) % compute steps/min 

