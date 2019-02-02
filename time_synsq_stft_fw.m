function [STFT,Vx,time,freqr]=time_synsq_stft_fw(sig,fs,sigma,m)
% Time-reassigned Synchrosqueezing Transform
%
% Input:
%  sig: vector of signal samples
%  fs: sampling frequency
%  sigma: standard deviation of the gaussin window
%  m: hop size of window
%    
% Output:
%  STFT: Short-time Fourier Transform of sig
%  Vx: Time-reassigned Synchrosqueezing Transform of x (columns associated with time t)
%  time: time axis of Vx
%  freqr: frequency axis of Vx
% 
% Modified by Dong He
% Date: 2017.07
% Email: hedong@stu.xjtu.edu.cn
%
% He D, Cao H, Wang S, et al. Time-reassigned synchrosqueezing transform: 
% The algorithm and its applications in mechanical signal processing[J]. 
% Mechanical Systems and Signal Processing, 2019, 117: 255-279.

N = length(sig);

if mod(N,2)
    sig(end)=[];   % NΪ����
    N=length(sig);
end

time=(0:m:N-1)/fs;
tsig=time'.*sig;
f = [zeros(N,1); sig; zeros(N,1)];  %���غ���ź�x(t)
% �����ź� t * x(t) 
tf= [zeros(N,1); tsig; zeros(N,1)];  %���غ���ź�t*x(t)

freqr=(0:N-1)*fs/N;    %����Ƶ���Ữ�ֽ��
% Initialize output matrix,
nw      = floor(N ./ m);
% ����N��һ�봦�Ƕ���
halfN=N/2+1;   % NΪż��

STFT    = zeros(halfN,nw);
tSTFT   = zeros(halfN,nw);
% Gaussian window function
w = N/2;  
ix     = ((-w):w);
t_win = ix/fs;
win = (pi*sigma^2)^(-0.25).* exp(-(t_win/sigma).^2/2);
win = win(:);

for u=1:nw
    tstart = 1+ (u-1)*m;
    tim = N + tstart + ix;
    seg = f(tim);
    tseg = tf(tim);
    
    seg_STFT = seg.*win;
    totseg_STFT = zeros(1,3*N);
    totseg_STFT(tim) = seg_STFT;
    localspec_STFT = fft(totseg_STFT(N+1:2*N))*2/N;%
    STFT(:,u)  = localspec_STFT(1:halfN);      %STFT(u,ksi)

    seg_tSTFT = tseg.*win;
    totseg_tSTFT = zeros(1,3*N);
    totseg_tSTFT(tim) = seg_tSTFT;
    localspec_tSTFT = fft(totseg_tSTFT(N+1:2*N))*2/N;%
    
    tSTFT(:,u) = localspec_tSTFT(1:halfN);   
end

freqr=freqr(1:halfN);   % Ϊ��ʡ�ռ䣬Ƶ�������
gamma = 1e-6;       % ��ֵ��С�����ֵ���ԣ���ΪNAN������ʽ2-53��
CandidateGD =  real(tSTFT ./STFT  ) ;   

clear tSTFT localspec_tSTFT localspec_STFT

% Synchrosqueezing 
Vx = SynchroSqueezing(STFT,CandidateGD,time);

end

% ʱ�����ų���
function Vx = SynchroSqueezing(STFT,delay,time)
% STFT-based Synchrosqueezing 
%   input:
%       STFT: the TF representation by STFT
%       w: Candidate instantaneous frequency
%       freqr: the frequency associated with STFT
%   output:
%       Tx: `the synchrosqueezing result

if nargin ~= 3
    error('you should check the input parameter of SynSqu_STFT function');
end

[STFT_rows,STFT_cols] = size(STFT);
Vx = zeros(size(STFT));       % ���Ž����ռλ��

delta_t = time(2)-time(1);  % STFT_rows��Ƶ�ʵ���,STFT_cols��ʱ�����
k=zeros(size(delay));
for u=1:STFT_cols
    for fi=1:STFT_rows
        if (~isnan(delay(fi, u)) && (delay(fi,u)>0))
            % Find w_l nearest to w(a_i,b)         
            k(fi,u) = 1 + round(delay(fi,u)/delta_t);
            if ~isnan(k(fi,u)) && k(fi,u)>0 && k(fi,u) <= STFT_cols
                Vx(fi,k(fi,u)) = Vx(fi,k(fi,u))  + STFT(fi, u); % ����ֵ�ۼ�
            end
        end
    end 
end 
end