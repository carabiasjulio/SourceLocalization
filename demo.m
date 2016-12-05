lsb = [1,1,1];          % Minimum room x,y,x dimensions
usb = [14,22,5];        % Maximum room x,y,x dimensions

% Mics location (xyz, 1 mic per row)
mic_loc = [7,20,2; ...  
          7,16,2;  ...
          3,15,2;  ...
          11,15,2; ...
          7,9,2];
        
load('./Examples/bassoon_5mics.mat','y','fs');    % Ground-Truth (4,18,2)
%load('./Examples/clarinet_5mics.mat','y','fs');   % Ground-Truth (10,18,2)
%load('./Examples/saxphone_5mics.mat','y','fs');   % Ground-Truth (4,13,2)
%load('./Examples/violin_5mics.mat','y','fs');     % Ground-Truth (10,13,2)

% Take a chuck from the signal (fast computation for demo)
pos_ini = 10000;
pos_fin = pos_ini + 2*fs; % 2 seconds from pos_ini

% Compute SRP_PHAT
[finalpos,finalsrp]=srpphat(y(pos_ini:pos_fin,:), mic_loc);%, fs, lsb, usb);

% Display the source location, just xy coord to simplify
room = zeros(usb(1),usb(2));
for k=1:numel(finalsrp),
    room(round(finalpos(k,1)),round(finalpos(k,2)))=max(room(round(finalpos(k,1)),round(finalpos(k,2))),finalsrp(k));
end;
figure;  imagesc(room');axis xy;
