%% Return the radial wave function
% 
function result = radial_wavfn(n, l, h, XorR)

    X = getX(n,l,h);
    if XorR == 1
        R = XtoR(X,n,h);
        result = R;
        %int = check_norm(result,n,h,1)
    else
        result = X;
        %int = check_norm(result,n,h,0)
    end
    
end

%% Check if input wave function is normalized
function result = check_norm(wavfn,n,h,XorR)
    int = 0;
    if XorR == 1
        for i = 1:((ceil(log(2*n*(n+15)))/h)+1)
            int = int + (wavfn(i)^2)*exp(2*h*(i-1))*exp(i*h)*...
                (1-exp(-h));
        end       
    else
        for i = 1:((ceil(log(2*n*(n+15)))/h)+1)
            int = int + (wavfn(i)^2)*exp(2*h*(i-1))*h;
        end
    end
    result = int;
end

%% Compute g(x)
function G = g(x,n,l)
    G = 2*exp(2*x)*(-1/exp(x) + (2*n^2)^(-1)) + (l+1/2)^2;
end

%% Return normalized wavefunction X
function ret_val = getX(n,l,h)
    rs = 2*n*(n+15);
    xs = ceil(log(rs));
    
    % Initialize the arrays containing the graph of each function
    X_array = zeros(1,xs/h +1);
    X_array(end) = 1e-5;
    X_array(end-1) = -1e-4;
    
    g_array = zeros(1,xs/h +1);
    g_array(end) = g(xs,n,l);
    g_array(end-1) = g(xs-h,n,l);
    
%    for i = 2:((xs/h)) %<-- without stopping criteria
        
    if l < 2
        end_x = 2*h;
    else
        end_x = 0.6*log(n*(n-sqrt(n^2-l*(l+1))));
    end
    ex = xs;
    i = 2;
    while ex >= end_x
        
        g_array(end-i) = g(xs-i*h,n,l);
        
        % Numerov Algorithm
        X_array(end-i) = (X_array(end-i+2)*(g_array(end-i+2) ...
           - 12/h^2) + X_array(end-i+1)*(10*g_array(end-i+1) ...
           + 24/h^2)) / (12/h^2 - g_array(end-i));
        
        i = i+1;
        ex = ex-h;
    end
    
    % Normalize
    A = 0;
    for i = 1:((xs/h)+1)
        A = A + h*exp(2*h*(i-1))*(X_array(i)^2);
    end
    A = 1/sqrt(A);
    ret_val = A*X_array;
end

%% Map X to R coordinates
function result = XtoR(X,n,h)
    
    rs = 2*n*(n+15);
    xs = ceil(log(rs));
    
    R_array = zeros(1,xs/h +1);
    R_array(end) = X(end)*exp(-(xs)/2);
    R_array(end-1) = X(end-1)*exp(-(xs-h)/2);
    
    for i = 2:(xs/h)
        R_array(end-i) = X(end-i)*exp(-(xs-i*h)/2);
    end
    
    % Return the graph of R
    result = R_array;
end
