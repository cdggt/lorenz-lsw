function x = lorenz_rk4(x,dt)
%LORENZ_RK4 4th-order Runge-Kutta time stepping scheme for integrating the
%1963 Lorenz dynamics.
%
% Inputs:
%
%   x    : 3x1 tuple containing the state of the system
%   dt   : time step

v1 = timederivative(x);
v2 = timederivative(x + dt/2 * v1); 
v3 = timederivative(x + dt/2 * v2); 
v4 = timederivative(x + dt * v3);

x = x + (v1 + 2*v2 + 2*v3 + v4) * dt/6;

end

function v = timederivative(x)

v = [
    10*(x(2)-x(1));
    x(1)*(28-x(3))-x(2); 
    x(1)*x(2)-8/3*x(3);
];
    
end