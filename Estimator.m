function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    % Initial state mean
    posEst = [0, 0]; % 1x2 matrix
    linVelEst = [0, 0]; % 1x2 matrix
    oriEst = [0]; % 1x1 matrix
    windEst = [0]; % 1x1 matrix
    driftEst = [0]; % 1x1 matrix
    
    % Initial state variance
    posVar = [power(estConst.StartRadiusBound, 4)*pi/4, power(estConst.StartRadiusBound, 4)*pi/4]; % 1x2 matrix
    linVelVar = [0, 0]; % 1x2 matrix
    oriVar = [estConst.RotationStartBound/3]; % 1x1 matrix
    windVar = [estConst.WindAngleStartBound/3]; % 1x1 matrix
    driftVar = [estConst.GyroDriftStartBound]; % 1x1 matrix
    
    % Estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar(1), posVar(2), linVelEst(1), linVelEst(1), oriEst(1), windEst(1), driftEst(1)]);
    % Estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst];
    % Time of last update
    estState.tm = tm;
end

if (tm ~= 0)
    %% Estimator iteration.
    % get time since last estimator update
    dt = tm - estState.tm;
    estState.tm = tm; % update measurement update time
    tspan = linspace(tm-dt, tm, 50);
%     tspan = [tm-dt, tm];

    % prior update
    x0 = estState.xm;
    [t, x] = ode45(@(t, x) xodefcn(t, x, actuate(1), actuate(2), estConst), tspan, x0');
    xp = x(end, :);
 
    P0 = estState.Pm;
    [T, P] = ode45(@(T, P) podefcn(T, P, t, x, actuate(1), actuate(2), estConst), tspan, P0(:));
    Pp = P(end, :);
    Pp = reshape(Pp, [7, 7]);

    % measurement update
    H = zeros(5, 7);
    H(4, 5) = 1;
    H(4, 7) = 1;
    H(5, 5) = 1;
    H(1, 1) = (xp(1) - estConst.pos_radioA(1))/sqrt(power(xp(1)-estConst.pos_radioA(1), 2) + power(xp(2)-estConst.pos_radioA(2), 2));
    H(1, 2) = (xp(2) - estConst.pos_radioA(2))/sqrt(power(xp(1)-estConst.pos_radioA(1), 2) + power(xp(2)-estConst.pos_radioA(2), 2));
    H(2, 1) = (xp(1) - estConst.pos_radioB(1))/sqrt(power(xp(1)-estConst.pos_radioB(1), 2) + power(xp(2)-estConst.pos_radioB(2), 2));
    H(2, 2) = (xp(2) - estConst.pos_radioB(2))/sqrt(power(xp(1)-estConst.pos_radioB(1), 2) + power(xp(2)-estConst.pos_radioB(2), 2));
    H(3, 1) = (xp(1) - estConst.pos_radioC(1))/sqrt(power(xp(1)-estConst.pos_radioC(1), 2) + power(xp(2)-estConst.pos_radioC(2), 2));
    H(3, 2) = (xp(2) - estConst.pos_radioC(2))/sqrt(power(xp(1)-estConst.pos_radioC(1), 2) + power(xp(2)-estConst.pos_radioC(2), 2));

    M = eye(5);
    R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.DistNoiseC, estConst.GyroNoise, estConst.CompassNoise]);
    hk = [sqrt(power(xp(1)-estConst.pos_radioA(1), 2) + power(xp(2)-estConst.pos_radioA(2), 2)); 
          sqrt(power(xp(1)-estConst.pos_radioB(1), 2) + power(xp(2)-estConst.pos_radioB(2), 2));
          sqrt(power(xp(1)-estConst.pos_radioC(1), 2) + power(xp(2)-estConst.pos_radioC(2), 2));
          xp(5) + xp(7);
          xp(5)];

    z = sense';
    if sense(3) == inf
        H(3, :) = [];
        M = eye(4);
        R = diag([estConst.DistNoiseA, estConst.DistNoiseB, estConst.GyroNoise, estConst.CompassNoise]);
        z(3, :) = [];
        hk(3, :) = [];
    end
    
    K = Pp*H' / (H*Pp*H'+M*R*M');
    xm = xp' + K*(z - hk);
    Pm = (eye(7)-K*H)*Pp;

    % Get resulting estimates and variances
    % Output quantities
    estState.xm = xm;
    estState.Pm = Pm;

    posEst = xm(1:2);
    linVelEst = xm(3:4);
    oriEst = xm(5);
    windEst = xm(6);
    driftEst = xm(7);

    posVar = [Pm(1, 1), Pm(2, 2)];
    linVelVar = [Pm(3, 3), Pm(4, 4)];
    oriVar = [Pm(5, 5)];
    windVar = [Pm(6, 6)];
    driftVar = [Pm(7, 7)];
end

end

function dxdt = xodefcn(t, x, ut, ur, estConst)

dxdt = zeros(7, 1);
dxdt(1) = x(3);
dxdt(2) = x(4);
dxdt(3) = cos(x(5)) * (tanh(ut) - estConst.dragCoefficientHydr*(x(3)*x(3)+x(4)*x(4))) ...,
    - estConst.dragCoefficientAir*(x(3)-estConst.windVel*cos(x(6))) * sqrt(power(x(3)-estConst.windVel*cos(x(6)), 2) + power(x(4)-estConst.windVel*sin(x(6)), 2));
dxdt(4) = sin(x(5)) * (tanh(ut) - estConst.dragCoefficientHydr*(x(3)*x(3)+x(4)*x(4))) ...,
    - estConst.dragCoefficientAir*(x(4)-estConst.windVel*sin(x(6))) * sqrt(power(x(3)-estConst.windVel*cos(x(6)), 2) + power(x(4)-estConst.windVel*sin(x(6)), 2));
dxdt(5) = estConst.rudderCoefficient*ur;
dxdt(6) = 0;
dxdt(7) = 0;

end

function dpdt = podefcn(t, p, tx, x, ut, ur, estConst)
x = interp1(tx, x, t);

Sx = x(3);
Sy = x(4);
phi = x(5);
rou = x(6);
Cdh = estConst.dragCoefficientHydr;
Cda = estConst.dragCoefficientAir;
Cw = estConst.windVel;
Cr = estConst.rudderCoefficient;
squareroot = sqrt((Sx-Cw*cos(rou))*(Sx-Cw*cos(rou)) + (Sy-Cw*sin(rou))*(Sy-Cw*sin(rou)));

A = zeros(7, 7);
A(1, 3) = 1;
A(2, 4) = 1;
A(3, 3) = -2*cos(phi)*Cdh*Sx - Cda*(squareroot + (Sx-Cw*cos(rou))*(Sx-Cw*cos(rou))/squareroot);
A(4, 3) = -2*sin(phi)*Cdh*Sx - Cda*(squareroot + (Sx-Cw*sin(rou))*(Sx-Cw*cos(rou))/squareroot);
A(3, 4) = -2*cos(phi)*Cdh*Sy - Cda*(Sx-Cw*cos(rou))*(Sy-Cw*sin(rou))/squareroot;
A(4, 4) = -2*sin(phi)*Cdh*Sy - Cda*(Sx-Cw*cos(rou))*(Sy-Cw*sin(rou))/squareroot;
A(3, 5) = -sin(phi)*tanh(ut) + sin(phi)*Cdh*(Sx*Sx+Sy*Sy);
A(4, 5) = cos(phi)*tanh(ut) - cos(phi)*Cdh*(Sx*Sx+Sy*Sy);
A(3, 6) = -Cda*Cw*sin(phi)*squareroot - Cda*(Sx-Cw*cos(rou))*((Sx-Cw*cos(rou))*Cw*sin(rou)-(Sy-Cw*sin(rou))*Cw.*cos(rou))/squareroot;
A(4, 6) = Cda*Cw*cos(phi)*squareroot - Cda*(Sx-Cw*sin(rou))*((Sx-Cw*cos(rou))*Cw*sin(rou)-(Sy-Cw*sin(rou))*Cw*cos(rou))/squareroot;

L = zeros(7, 4);
L(3, 1) = -cos(phi)*Cdh*(Sx*Sx+Sy*Sy);
L(4, 1) = -sin(phi)*Cdh*(Sx*Sx+Sy*Sy);
L(5, 2) = Cr*ur;
L(6, 3) = 1;
L(7, 4) = 1;

Q = diag([estConst.DragNoise, estConst.RudderNoise, estConst.WindAngleNoise, estConst.GyroDriftNoise]);

P = reshape(p, [7, 7]);
dP = A*P + P*A' + L*Q*L';

dpdt = dP(:);
end