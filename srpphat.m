function [finalpos,finalsrp]=srpphat(x, mic_loc, fs, lsb, usb)
%% This function uses SRP-PHAT
%% Inputs:
%%% 1) x is the multi-channel input data (samples x channels)
%%% 2) mic_loc is the microphone 3D-locations (M x 3) ( in meters)
%%% 3) fs: sampling rate (Hz)
%%% 4) lsb: a row-vector of the lower rectangular search boundary, e.g., [1 1 1] (meters)
%%% 5) usb: row-vector of the upper rectangular search boundary, e.g., [14 22 5] (m)
%% Outputs:
%%% 1) finalpos: estimated location of the source
%%% 2) finalsrp: srp-phat value of the point source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Julio Carabias
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ngrid = 1810;       % Number of TDoA to consider
win = 4096;         % window size (for the spectrogram)
sp_resolution= .5;  % In meters, spatial resolution for the searching procedure
                    %                   (low value = more computational demand)

%% Initialize variables:
if nargin < 5, usb=max(mic_loc); end
if nargin < 4, lsb=[1 1 1]; end
if nargin < 3, fs=44100; end

M = size(mic_loc,1); %%% number of microphones (also number of channels of the input X).
np=M*(M-1)/2;        %%% number of independent pairs

%% Determine the maximum end-fire length (in samples) of a microphone pairs:
mdist=pdist(mic_loc);
efmax=max(mdist);


%% Doing the GCC-PHAT:

X = stft_multi(x',win);
[nbin,nfram,nchan] = size(X);

tau_grid = linspace(-efmax/343,efmax/343,ngrid);
f = linspace(0,(fs/2),nbin)';

R = zeros(ngrid,np);
p=1;
for ii=1:M-1
    X1 = X(:,:,ii);
    for jj=ii+1:M,
        spec = zeros(nbin,nfram,ngrid);
    
        X2 = X(:,:,jj);      
        P = X1.*conj(X2);
        P = P./(abs(P)+eps);
        
        for ind = 1:ngrid,
            EXP = exp(-2*1i*pi*tau_grid(ind)*f);
            spec(:,:,ind) = real(bsxfun(@times,P,EXP));
        end
        R(:,p) = shiftdim(max(sum(spec,1),[],2));
        p=p+1;
    end;
end;


%% Initialize to do SRC:
bstart = lsb;
bend = usb;

if bstart==bend  %%%check if the lower search-boundary is the same as the upper one, i.e., V_{search}=0
    fprintf('Search volume is 0! Please expand it');
end

%%Doing SRC:
iter = 1;
niter=numel(bstart(1):sp_resolution:bend(1))*numel(bstart(2):sp_resolution:bend(2));
yval = zeros(niter,1);
position = zeros(niter,3);
for pos_x=bstart(1):sp_resolution:bend(1),
    for pos_y=bstart(2):sp_resolution:bend(2),
        for pos_z=bstart(3):sp_resolution:bend(3),
            % pos_z=2;
            [yval(iter),position(iter,:)] = fe([pos_x,pos_y,pos_z],efmax,mic_loc,R);
            iter=iter+1;
        end;
    end;
end;

% Savinf outputs
[yval_sort,x_sort] = sort(yval,'descend');
finalpos=position(x_sort(1:end),:);  %%%Final source location estimate
finalsrp=yval_sort(1:end);  %%%Final source's SRP-PHAT value (normalized by number of pairs)

return;


%% Functional Evaluation sub-routine (calculate SRP-PHAT value for a point in the search space):

function [yval1,position1] = fe(x,efmax,mic_loc,R)
%%%This function evaluates an 'fe' in the search space, gives back the 'fe' position and its SRP-PHAT value.
M=size(mic_loc,1); %%%number of mics
ngrid=size(R,1);
np=size(R,2);

%%%Find the distances from M-microphones to the point:
a1=ones(M,1);
xx1=a1*x;
xdiff1=xx1-mic_loc;
dists=sqrt(sum(xdiff1.*xdiff1,2));

%%%%Differences in distances:
ddiffs_ones=ones(M,1)*dists';
ddm=ddiffs_ones-ddiffs_ones';

%%% Calculate the TDOA index:
% v=nonzeros(tril(ddm,0))';
v=nonzeros(tril(ddm+eps,-1))'; % MODIFICADO

% Pasa las distancias a indices
interidx = round(-v*(ngrid-1)/(2*efmax) + (ngrid-1)/2 + 1);

%%%Pull out the GCC-PHAT value corresponding to that TDOA:
v1=R(sub2ind([ngrid np],interidx,1:np));

%%% SRP-PHAT value of the point:
yval1=sum(v1);
position1=x;

return;
