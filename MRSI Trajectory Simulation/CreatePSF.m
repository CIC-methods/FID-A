prompt1 = 'enter k(s) values for your function';
prompt2 = 'enter s values for your function';
s = input(prompt2);
k_s = input(prompt1);
figure(1);
hold on;
subplot(2,1,1)
plot(s,k_s)
hold off
MakePSF(k_s, s);   
function MakePSF(varargin)
    if(nargin==2)
        k_s = varargin{1};
        s = varargin{2};
        dT = mean(diff(s));
        f_n = 1/dT
        SW = 2*f_n;
        n = length(s)
        w = (-n/2:n/2-1)*(SW/n);
        y = fft(k_s)
        y = real(y);
        yshift = fftshift(y);
        hold on;
        subplot(2,1,2);
        plot(w,yshift)
        hold off;
    elseif(nargin==3)
        
        
        
    end
end

function

