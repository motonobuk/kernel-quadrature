% Plot results of quadrature experiments in the following paper:
% M. Kanagawa, B. Sriperumbudur and K. Fukumizu. 
% Convergence analysis of deterministic kernel-based quadrature rules in misspecified settings.
% arXiv preprint, arXiv:1709.00147 [math.NA] 2017.

clear all
close all;
seed=1;
randn('state',seed);
rand('state',seed);
tol=1e-24;

d = 1; % Dimension: currently the following implementation can only be used for d = 1
if rem(d,2) == 0
   error('The dimension d must be odd.'); 
end

r = 4;  % Smoothness of the Wendland kernel used for computing quadrature weights
k = r - (d+1)/2;
if k < 0
   error('r should be larger than or equal to (d+1)/2'); 
end

%DesignPoints = 'Equidist'; % Equally spaced design points
DesignPoints = 'Irregular'; % Irregularly spaced design points

delta = 0.1; % Bandwidth of the Wendland kernel

a = 0; % Support of the distribution [a,b];
b = 1;

list_N = ceil(logspace(log10(50), log10(400), 100)); % Number of design points
len_N=length(list_N);

list_s = [1,2,3,4]; % Smoothness of test integrands
len_s = length(list_s);

if min(list_s) < (d+1)/2
   error('s should be larger than or equal to (d+1)/2'); 
end

errors = zeros(len_s,len_N); % Worst case errors for test integrands
WCEs = zeros(1,len_N); % Worst case errors in the correctly specified RKHS
AbsSumWeights = zeros(1,len_N); % Sums of absolute weights
FillDists = zeros(1,len_N); % Fill distances
SepRadius = zeros(1,len_N); % Separation radius

for iN=1:len_N  
    N=list_N(iN);
    switch DesignPoints
        case 'Equidist'
            X = (b-a) * (0:N-1)'/(N-1) + a;     
        case 'Irregular'
            X = (b-a) * (0:N-1)'/(N-1) + a;
            ind_even = (mod(1:N,2) == 0);
            X(ind_even) = X(ind_even) - 1/(N-1) + 1/(N+1)^2;
    end
    R = pdist2(X,X);   % Distance matrix
    G = sparse( Wendland(R/delta,d,k) );  % Gram matrix: (wendland kernel)
    Z = kmeanval_Wendland_unif(X,delta,a,b,k); % Evaluations of the kernel mean
    W = G\Z;     % Weight vector for quadrature
    WCEs(iN) = WCE_Wendland_unif(X, W, delta, a, b, k );
    AbsSumWeights(iN) = sum(abs(W));
    R_inf = R + diag( Inf* ones(1, size(R,1) ) );
    SepRadius(iN) = min(min(R_inf)) / 2;
    FillDists(iN) = max(diag(R,1)) / 2;
    for is=1:len_s      
        k_test = list_s(is) - (d+1)/2;
        errors(is,iN) = WCE_Wendland_unif(X, W, delta, a, b, k_test );
    end
end




figure;
Col=['r','b','g','m'];
for i=1:len_s  

    loglog(list_N,errors(i,:),'linewidth',2,'Color',Col(i));    
    if i==1    
        set(gca,'FontName','Arial');
        set(gca,'FontSize',20);
        hold on;
    end
end
loglog(list_N, WCEs, 'linewidth',2, 'Color', 'k'); hold on;
axis([list_N(1) list_N(len_N) min(min(errors)) max(max(errors))]); 
list_a = zeros(len_s,1);
for i=1:len_s 
    ID_nz = errors(i,:)>tol;
    p = polyfit( log(list_N(ID_nz)), log(errors(i,ID_nz)), 1);
    a = p(1);
    b = p(2);
    loglog(list_N,exp(b)*list_N.^a,'-.','linewidth',1.75, 'Color',Col(i)); hold on;
    list_a(i) = a;
end
list_a = round(list_a,3);

ID_nz = WCEs>tol;
p = polyfit( log(list_N(ID_nz)), log(WCEs(ID_nz)), 1 );
a_WCE = p(1);
b_WCE = p(2);
loglog(list_N, exp(b_WCE)*list_N.^a_WCE,'k-.','linewidth',1.75); hold on
a_WCE = round(a_WCE,3);

ID_nz = SepRadius > tol;
p = polyfit( log(list_N(ID_nz)), log(SepRadius(ID_nz)), 1 );
a_Sep = p(1);

ID_nz = FillDists > tol;
p = polyfit( log(list_N(ID_nz)), log(FillDists(ID_nz)), 1 );
a_Fil = p(1);

p = polyfit( log(list_N), log(AbsSumWeights), 1 );
a_weights = p(1);

legend({strcat('s = 1:', char(8239), num2str(list_a(1))),...
    strcat('s = 2:', char(8239), num2str(list_a(2))),...
    strcat('s = 3:', char(8239), num2str(list_a(3))),...
    strcat('s = 4:', char(8239), num2str(list_a(4))),...
    strcat('r = ', char(8239), char(8239), num2str(r), ':', char(8239), num2str(a_WCE))},...    
    'Location','SouthWest','FontSize', 20);
xlabel('Sample size');
ylabel('Worst case error');

title( ['Fill Dist: ', num2str(round(a_Fil,2)), '. Sep Rad: ',...
    num2str(round(a_Sep,2)), '. Weights: ', num2str(round(a_weights,2)), '.' ] )
hold on;

