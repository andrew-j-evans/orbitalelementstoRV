format longG
% Declaration of variables
R = [-7953.8073 -4174.5370 -1008.9496];
V = [3.6460035 -4.9118820 -4.9193608];

%Converseion to orbital elements
[a,e,i,RA,w,nu,mu] = RV2COE(R,V)

%Conversion to RV
[R2,V2] = COE2RV(a,e,i,RA,w,nu,mu)

%Function to convert classical orbital elements to position and velocity
%vectors
function [a,e,i,RA,w,nu,mu] = RV2COE(R,V)
    mu = 398600;

    %Magnitude
    r = norm(R);
    v = norm(V);
    
    % Determination of Classical Orbital Elements
    % 1) Semi Major Axis
        epsillon = (v*v) / 2 - mu / r;
        a = - mu / (2 * epsillon);

    % 2) Eccentricity
        RVDot = dot(R,V);
        ebar = (1 / mu) * (R*(v*v-mu/r) - RVDot*V);
        e = norm(ebar);
    
    % 3) True Anomaly
        hbar = cross(R,V);
        h = norm(hbar);
        i = acosd(hbar(3)/h);
   
   % 4) Right Ascention of the Ascending Node
        kbar = [0,0,1];
        nbar = cross(kbar,hbar);
        n = norm(nbar);
        RA = acosd(nbar(1)/n);

        if nbar(2) < 0
            RA = 360 - RA;
        end

   % 5) Argument of Perigee
        nedot = dot(nbar,ebar);
        w = acosd(nedot/(n*e));
        if ebar(3) < 0
            w = 360 - w;
        end

   % 6) True Anomaly

        eRdot = dot(ebar, R);
        nu = acosd(eRdot/(e*r));   
        
        if RVDot < 0
            nu = 360 - nu;
        end


end

%Function to convert position and velocity vectors into their classical
%orbital elements
function [R,V] = COE2RV(a,e,i,RA,w,nu,mu)
    if(i == 0 || e == 0)
        disp("Warning!!! This is a special case! i or e = 0!")
    end
    P = a*(1-e*e);
    r = P / (1 + e*cosd(nu));

    snu = sind(nu);
    cnu = cosd(nu);

    Rpqw = [r*cnu; r*snu; 0];
    mup = sqrt(mu / P);
    Vpqw = [mup*(-snu); mup*(e+cnu); 0];
    
    cRA = cosd(RA);
    sRA = sind(RA);
    cw = cosd(w);
    sw = sind(w);
    ci = cosd(i);
    si = sind(i);
 
    A = [cRA*cw - sRA*sw*ci, -cRA*sw-sRA*cw*ci, sRA*si;
         sRA*cw+cRA*sw*ci, -sRA*sw+cRA*cw*ci, -cRA*si;
         sw*si, cw*si, ci
    ];

    R = A * Rpqw;
    V = A * Vpqw;

end
