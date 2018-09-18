function Y = Wendland( R, d, k )
% Wendland function
% Formula from Corollary 9.14 in Wendland's book (2005).

% Inputs
% R: Distance matrix. It can be a vector.
% d: dimension of the Eulclidian space.
% k: smoothness of the function. k = 1,2,3.

% Output
% Y: matrix containing evaluations of the Wendland function with respect to
% mat_r


if min(R) < 0
    warning('Input must be non-negative')
end

ell = floor(d/2) + k + 1;

switch k
    case 0

        Y = max(1-R,0).^(floor(d/2)+1);          
        
        
    case 1
        
        C1 = ell + 1;
        C0 = 1;

        Y = ( max(1-R,0).^(ell+1) ) .* (C1*R + C0);        
        
        
    case 2

        C2 = ell^2 + 4*ell + 3;
        C1 = 3*ell + 6;
        C0 = 3;

        Y = ( max(1-R,0).^(ell+2) ) .* (C2*R.^2 + C1*R + C0);

        
    case 3

        C3 = ell^3 + 9*ell^2 + 23*ell + 15;
        C2 = 6*ell^2 + 36*ell + 45;
        C1 = 15*ell + 45;
        C0 = 15;

        Y = ( max(1-R,0).^(ell+3) ) .* (C3*R.^3 + C2*R.^2 + C1*R + C0);

        
    otherwise
        error('k must be 0, 1, 2, or 3.')
end



end

% % Test 
% N = 10;
% X = (-2*N:2*N)'/N;
% delta = 1;
% k = 0;
% 
% R = pdist2(X,X);
% 
% % Y = Wendland(abs(X)/delta,1,k);
% % plot(X,Y);
% 
% G = Wendland(R,1,k);
% [~,D] = eig(G);
% diag(D)

