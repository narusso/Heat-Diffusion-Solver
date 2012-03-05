function [] = BE_MG()
  L      = 2;                  % length of 1D object in m
  nx     = 255;                % number of divisions
  nsteps = 175000;             % number of time steps
  alpha  = 1.1234e-4;          % diffusivity of copper in m^2/s
  
  dx     = L/(nx+1);           % length of one segment in m
  dt     = .003;               % length of one time step in seconds
  C      = alpha*dt/dx^2;      % convenient constant
  x      = linspace(0,1,nx+2); % nx+2 divisions of [0,1] inclusive (simplification of [0,L])
  
  W    = .5*sin(x*13*pi) + .2*sin(x*50*pi); % wiggle function
  T    = exp(-(5*x-2.5).^2)+W; % initial temps: gaussian distribution with wiggle
  T(1) = 0; T(nx+2) = 0;       % temp at edges is started (and later fixed) at 0
  
  plot(x,T);                   % show initial temperatures along the length
  ylim([-.5 1.5]);
  pause(1);                    % pause while the figure is foregrounded and observed
  
  for n=1:nsteps
    % Kill high frequencies in Ah*x = bh
    Tnew = Jacobi(T,C,5);      % a few sweeps of (unweighted for now) Jacobi

    % Restriction
    A = spdiags([-C*ones(nx+2,1) (2*C*ones(nx+2,1) + 1) -C*ones(nx+2,1)],[-1 0 1], nx+2, nx+2);
    rh = A*Tnew' - T'; % Ae = A*x - b
    r2h = rh(1:2:end)';        % take every other value

    % Solve A2h*e2h = r2h using sweeps to convergence or whatever.
    e2h = Jacobi(r2h,C,500);   % Jacobi puts zeroes in at boundaries. Should it?

    % Interpolate e2h to fine grid. eh = I(h.2h)*e2h, add eh + uh
    eh = interp1(1:2:nx+2,e2h,1:nx+1);
    T(1:nx+1) = T(1:nx+1) + eh;% improve T

    % A few more Jacobi sweeps
    T = Jacobi(T,C,5);         % a few more sweeps of (unweighted for now) Jacobi

    if (~mod(n,500))           % if n is a multiple of 500
        n
        plot(x,T);             % plot the graph
        ylim([-.5 1.5]);       % keep the y-axis steady
        pause(.001);           % pause
    end
  end
end

function xnew = Jacobi(xold,C,MAX_ITER)
  S = size(xold);
  x = xold; % first guess from old solution
  xnew = zeros(S);
  for m=1:MAX_ITER
    for i=2:max(S)-1
      xnew(i) = C/(2*C+1)*(x(i-1) + x(i+1)) + 1/(2*C+1)*xold(i);
    end
    if mean(abs(xnew-x)) < 1.e-9
      break;
    end
    x = xnew;
  end
end
