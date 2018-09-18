function Z = kmeanval_Wendland_unif( X, delta, a, b, k )
% Returns evaluations of the kernel mean for a uniform distribution on [a,b] using a
% Wendland kernel with smoothness k.

% Input
% X: column vector containing the points to be evalauted
% delta: bandwidth of the Wendland kernel
% a,b: edges of the interval [a,b] on which the uniform distribuiton is
% defiend
% k: smoothness of the Wendland kernel. k = 0,1,2,3.

% Output
% z: column bector containing the evaluations of the kernel mean

if delta > (b-a)/2
    error('the bandwidth is too large')
end

% Constants
d = 1; % dimensionality of the input space
ell = floor(d/2) + k + 1;


switch k
    case 0
        n = length(X);

        % Three cases: 
        % 1: a <= x-delta  and x+delta <= b
        % 2: x - delta < a
        % 3: b < x + delta

        Z = zeros(n,1);

        % Case 1

        Ind1 = ( (a <= X-delta) & (X+delta <= b) );
        Z(Ind1) = (2*delta) / ((b-a) * (floor(d/2)+2));


        % Case 2
        Ind2 = (X-delta < a);
        V = (delta / ( (b-a) * (floor(d/2)+2) ) ) * (2 - (1 - (X-a)/delta).^(floor(d/2)+2) );
        Z(Ind2) = V(Ind2);


        % Case 3
        Ind3 = (b < X + delta);
        V = (delta / ( (b-a) * (floor(d/2)+2) ) ) * (2 - (1 - (b-X)/delta).^(floor(d/2)+2) );      
        Z(Ind3) = V(Ind3);
        
        
    case 1
        
        C1 = ell + 1;
        C0 = 1;

        A = (-C1)/(ell+3) + (C1+C0)/(ell+2);

        n = length(X);

        % Three cases: 
        % 1: a <= x-delta  and x+delta <= b
        % 2: x - delta < a
        % 3: b < x + delta

        Z = zeros(n,1);

        % Case 1

        Ind1 = ( (a <= X-delta) & (X+delta <= b) );
        Z(Ind1) = A * (2*delta)/(b-a);


        % Case 2
        Ind2 = (X-delta < a);

        V1 = ( (-C1)/(ell+3) ) * ( 1 - (1-(X-a)/delta).^(ell+3) );
        V0 = ( (C1+C0)/(ell+2) ) * ( 1 - (1-(X-a)/delta).^(ell+2) );
        V = V1+V0;

        Z(Ind2) = (delta/(b-a)) * (A + V(Ind2));



        % Case 3
        Ind3 = (b < X + delta);

        V1 = ( (-C1)/(ell+3) ) * ( 1 - (1-(b-X)/delta).^(ell+3) );
        V0 = ( (C1+C0)/(ell+2) ) * ( 1 - (1-(b-X)/delta).^(ell+2) );
        V = V1+V0;

        Z(Ind3) = (delta/(b-a)) * (A + V(Ind3));
        
               
        
        
    case 2

        C2 = ell^2 + 4*ell + 3;
        C1 = 3*ell + 6;
        C0 = 3;

        A = (C2)/(ell+5) + (-2*C2-C1)/(ell+4) + (C2+C1+C0)/(ell+3);

        n = length(X);

        % Three cases: 
        % 1: a <= x-delta  and x+delta <= b
        % 2: x - delta < a
        % 3: b < x + delta

        Z = zeros(n,1);

        % Case 1

        Ind1 = ( (a <= X-delta) & (X+delta <= b) );
        Z(Ind1) = A * (2*delta)/(b-a);



        % Case 2
        Ind2 = (X-delta < a);

        V2 = ( (C2)/(ell+5) ) * ( 1 - (1-(X-a)/delta).^(ell+5) );
        V1 = ( (-2*C2-C1)/(ell+4) ) * ( 1 - (1-(X-a)/delta).^(ell+4) );
        V0 = ( (C2+C1+C0)/(ell+3) ) * ( 1 - (1-(X-a)/delta).^(ell+3) );
        V = V2+V1+V0;

        Z(Ind2) = (delta/(b-a)) * (A + V(Ind2));



        % Case 3
        Ind3 = (b < X + delta);

        V2 = ( (C2)/(ell+5) ) * ( 1 - (1-(b-X)/delta).^(ell+5) );
        V1 = ( (-2*C2-C1)/(ell+4) ) * ( 1 - (1-(b-X)/delta).^(ell+4) );
        V0 = ( (C2+C1+C0)/(ell+3) ) * ( 1 - (1-(b-X)/delta).^(ell+3) );
        V = V2+V1+V0;

        Z(Ind3) = (delta/(b-a)) * (A + V(Ind3));
        
                
        
    case 3

        C3 = ell^3 + 9*ell^2 + 23*ell + 15;
        C2 = 6*ell^2 + 36*ell + 45;
        C1 = 15*ell + 45;
        C0 = 15;

        A = ( -C3/(ell+7) + (3*C3+C2)/(ell+6) + (-3*C3-2*C2-C1)/(ell+5) + (C3+C2+C1+C0)/(ell+4) );

        n = length(X);

        % Three cases: 
        % 1: a <= x-delta  and x+delta <= b
        % 2: x - delta < a
        % 3: b < x + delta

        Z = zeros(n,1);

        % Case 1

        Ind1 = ( (a <= X-delta) & (X+delta <= b) );
        Z(Ind1) = A * (2*delta)/(b-a);



        % Case 2
        Ind2 = (X-delta < a);

        V3 = ( -C3/(ell+7) ) * ( 1 - (1-(X-a)/delta).^(ell+7) );
        V2 = ( (3*C3+C2)/(ell+6) ) * ( 1 - (1-(X-a)/delta).^(ell+6) );
        V1 = ( (-3*C3-2*C2-C1)/(ell+5) ) * ( 1 - (1-(X-a)/delta).^(ell+5) );
        V0 = ( (C3+C2+C1+C0)/(ell+4) ) * ( 1 - (1-(X-a)/delta).^(ell+4) );
        V = V3+V2+V1+V0;

        Z(Ind2) = (delta/(b-a)) * (A + V(Ind2));



        % Case 3
        Ind3 = (b < X + delta);

        V3 = ( -C3/(ell+7) ) * ( 1 - (1-(b-X)/delta).^(ell+7) );
        V2 = ( (3*C3+C2)/(ell+6) ) * ( 1 - (1-(b-X)/delta).^(ell+6) );
        V1 = ( (-3*C3-2*C2-C1)/(ell+5) ) * ( 1 - (1-(b-X)/delta).^(ell+5) );
        V0 = ( (C3+C2+C1+C0)/(ell+4) ) * ( 1 - (1-(b-X)/delta).^(ell+4) );
        V = V3+V2+V1+V0;

        Z(Ind3) = (delta/(b-a)) * (A + V(Ind3));
        
        
    otherwise
        error('k must be eathier of 0, 1, 2, or 3.')

        
end

end

