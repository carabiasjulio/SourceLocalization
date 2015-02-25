DB_directory='./BBDD_Bach/';

%% Room Dimensions
room = zeros(14,22);

%% Room space lower to upper
lsb = [1,1,1];
usb = [14,22,5];

%% Mic position (case 4 mics)
mic_loc = [7,20,2;  ...
           7,16,2;  ...
           3,15,2;  ...   
           11,15,2];

%% Load signals
%load([DB_directory 'bassoon4_echo_av.mat'],'y','fs');   
load([DB_directory 'clarinet4_anechoid.mat'],'y','fs');  
%load([DB_directory 'saxphone4_echo_av.mat'],'y','fs');  
%load([DB_directory 'violin4_echo_av.mat'],'y','fs');

%% Mic position (case 5 mics)
% mic_loc = [7,20,2;  ...
%           7,16,2;  ...
%           3,15,2;  ...
%           11,15,2; ...
%           7,9,2];

%% Load signals
%load('bassoon5_anechoid.mat','y','fs');   
%load('clarinet5_anechoid.mat','y','fs');  
%load('saxphone5_anechoid.mat','y','fs');  
%load('violin5_anechoid.mat','y','fs');

pos_ini = 10000;
pos_fin = fs+17270;

[finalpos,finalsrp]=srpphat(y(pos_ini:pos_fin,:), mic_loc, fs, lsb, usb);
%[finalpos,finalsrp]=srpphat_sg(y(pos_ini:pos_fin,:), mic_loc, fs, lsb, usb);

for k=1:numel(finalsrp),
    room(round(finalpos(k,1)),round(finalpos(k,2)))=max(room(round(finalpos(k,1)),round(finalpos(k,2))),finalsrp(k));
end;


figure;  imagesc(room');axis xy;


