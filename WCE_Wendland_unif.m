function WCE = WCE_Wendland_unif( X, W, delta, a, b, k ) 
% Compute the worst case error when the kenel is Wendland and the
% distribution is uniform 



d = 1;
ell = floor(d/2) + k + 1;

switch k
    case 0
        
        B = ( 2*delta * (b-a-2*delta) ) / ( (b-a)^2 * (floor(d/2)+2) ) ...
            + ( (2*delta^2) / ((b-a)^2 * (floor(d/2)+2)) ) *  (2 - 1/(floor(d/2)+3) );
        
        % Distance matrix
        R = pdist2(X,X); 

        % Gram matrix: (wendland kernel)
        G = sparse( Wendland(R/delta,d,k) );

        % Evaluations of the kernel mean
        Z = kmeanval_Wendland_unif(X,delta,a,b,k);

        WCE = sqrt( W'*G*W - 2*W'*Z + B );           
    
    case 1
        
        C1 = ell + 1;
        C0 = 1;

        A1 = (-C1)/(ell+3) + (C1+C0)/(ell+2);
        A2 = (-C1)/((ell+3)*(ell+4)) + (C1+C0)/((ell+2)*(ell+3));
        B = ( (2*delta*(b-a-2*delta))/(b-a)^2 ) * A1 + ( (2*delta^2)/(b-a)^2 ) * ( 2*A1 - A2 );    
        
        % Distance matrix
        R = pdist2(X,X); 

        % Gram matrix: (wendland kernel)
        G = sparse( Wendland(R/delta,d,k) );

        % Evaluations of the kernel mean
        Z = kmeanval_Wendland_unif(X,delta,a,b,k);

        WCE = sqrt( W'*G*W - 2*W'*Z + B );             
    
    case 2
        
        C2 = ell^2 + 4*ell + 3;
        C1 = 3*ell + 6;
        C0 = 3;

        A1 = (C2)/(ell+5) + (-2*C2-C1)/(ell+4) + (C2+C1+C0)/(ell+3);
        A2 = (C2)/((ell+5)*(ell+6)) + (-2*C2-C1)/((ell+4)*(ell+5)) + (C2+C1+C0)/((ell+3)*(ell+4));
        B = ( (2*delta*(b-a-2*delta))/(b-a)^2 ) * A1 + ( (2*delta^2)/(b-a)^2 ) * ( 2*A1 - A2 );
        
        % Distance matrix
        R = pdist2(X,X); 

        % Gram matrix: (wendland kernel)
        G = sparse( Wendland(R/delta,d,k) );

        % Evaluations of the kernel mean
        Z = kmeanval_Wendland_unif(X,delta,a,b,k);

        WCE = sqrt( W'*G*W - 2*W'*Z + B );        
        
    case 3
        
        % Constants
        C3 = ell^3 + 9*ell^2 + 23*ell + 15;
        C2 = 6*ell^2 + 36*ell + 45;
        C1 = 15*ell + 45;
        C0 = 15;        
        
        A1 = (-C3)/(ell+7) + (3*C3+C2)/(ell+6) + (-3*C3-2*C2-C1)/(ell+5) + (C3+C2+C1+C0)/(ell+4);
        A2 = (-C3)/((ell+7)*(ell+8)) + (3*C3+C2)/((ell+6)*(ell+7)) + (-3*C3-2*C2-C1)/((ell+5)*(ell+6)) + (C3+C2+C1+C0)/((ell+4)*(ell+5));
        B = ( (2*delta*(b-a-2*delta))/(b-a)^2 ) * A1 + ( (2*delta^2)/(b-a)^2 ) * ( 2*A1 - A2 );
        
        % Distance matrix
        R = pdist2(X,X); 

        % Gram matrix: (wendland kernel)
        G = sparse( Wendland(R/delta,d,k) );

        % Evaluations of the kernel mean
        Z = kmeanval_Wendland_unif(X,delta,a,b,k);

        WCE = sqrt( W'*G*W - 2*W'*Z + B );
        
    otherwise
        error('k must be 0, 1, 2, or 3!!!')
end

end

