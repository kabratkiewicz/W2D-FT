function Ren = Renyi_Entropy_4D(W2DFT,t,r,w,n,alpha)
%   
%   Inputs:
%	W2DFT : (M,N,K,L) 4-D W2DFT function.
%	t : a time vector   (default : (1:K)).
%   r : a space vector (default : (1:L)).
%   w : a frequency vector (default : (1:M)).
%	n : a spatial frequency vector	(default : (1:N)).	
%	ALPHA : Renyi measure order	(default : 3).
%   
%   Outputs:
%   Ren : Alpha-order Renyi entropy
%   ALPHA = 1: limit case, the outcomes will be Shannon entropy


if (nargin == 0)
 error('At least one parameter required');
end

[M,N,K,L] = size(W2DFT);
if (nargin == 1)
 t = 1:K; 
 r = 1:L; 
 w = (1:M)'; 
 n = (1:N)'; 
 alpha=3;
end

w = sort(w); %sort frequency vector in ascending order such that the first 
n = sort(n);
%row TFR must correspond to the lower frequencies
%W2DFT dimensions: (NFFT_omega, NFFT_eta, T, R)
W2DFT = W2DFT./trapz(w,trapz(n,trapz(t,trapz(r,W2DFT,4),3),2));
% Normalisation W2DFT;
%trapz function is used to calculate 2D integral of matrix being a cut of W2DFT according
%to abscissa X and ordinate Y

if alpha == 1 % limit case case: Shannon entropy
 if (min(min(W2DFT))<0)
     error('distribution with negative values => alpha=1 not allowed');
 else
     Ren=-trapz(w,trapz(n,trapz(t,trapz(r,W2DFT.*log2(W2DFT+eps),4),3),2));
 end
else % Renyi entropy
    Ren=1/(1-alpha)*log2(trapz(w,trapz(n,trapz(t,trapz(r,W2DFT.^alpha,4),3),2))+eps);
end



