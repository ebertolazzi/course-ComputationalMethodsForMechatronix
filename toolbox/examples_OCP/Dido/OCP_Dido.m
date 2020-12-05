%--------------------------------------------------------------------------%
%                                                                          %
%  Copyright (C) 2018                                                      %
%                                                                          %
%         , __                 , __                                        %
%        /|/  \               /|/  \                                       %
%         | __/ _   ,_         | __/ _   ,_                                %
%         |   \|/  /  |  |   | |   \|/  /  |  |   |                        %
%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       %
%                           /|                   /|                        %
%                           \|                   \|                        %
%                                                                          %
%      Enrico Bertolazzi                                                   %
%      Dipartimento di Ingegneria Industriale                              %
%      Universita` degli Studi di Trento                                   %
%      email: enrico.bertolazzi@unitn.it                                   %
%                                                                          %
%--------------------------------------------------------------------------%

classdef OCP_Dido < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x0;
    x1;
    L;
  end

  methods ( Hidden = true )
    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, T )
      u = U(1);
      v = U(2);
      RES = [ u; v; sqrt(u^2+v^2) ];
    end
    %
    function RES = JAC( self, nseg, t, X, U, T )
      RES = sparse( 3, 5 );
      u = U(1);
      v = U(2);
      tmp = sqrt(u^2+v^2+eps);
      RES(1:3,4:5) = [ 1, 0; 0, 1; u/tmp, v/tmp; ];
    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( 3, 5 );
      RES(1:3,4:5) = 1;
    end
    %
    function H = HESS( self, nseg, t, X, U, T, L )
      H   = sparse(5,5);
      u   = U(1);
      v   = U(2);
      tmp = 1/sqrt(u^2+v^2+eps)^3;
      H(4:5,4:5) = tmp*[ v^2, -u*v; -u*v, u^2 ];
    end
    %
    function H = HESS_pattern( self )
      H = sparse(5,5);
      H(4:5,4:5) = 1;
    end
  end

  methods

    function self = OCP_Dido( L, x0, x1 )
      nx  = 3; % number of states
      nu  = 2; % number of controls
      np  = 0; % number of free parameters
      nbc = 4; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
      self.x0 = x0;
      self.x1 = x1;
      self.L  = L;
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
    end

    function info = solve( self )

      utot = (self.N-self.nseg)*self.nu;
      xtot = self.N*self.nx;

      xones = ones(1,xtot);
      options.lb = [ -xones*Inf, -Inf*ones(1,utot) ]; % Lower bound on the variables.
      options.ub = [  xones*Inf,  Inf*ones(1,utot) ]; % Upper bound on the variables.

      % The constraint functions are bounded to zero
      dim = (self.N-self.nseg)*self.nx + (self.nseg-1)*self.njmp + self.nbc;
      options.cl = zeros(1,dim); %  constraints
      options.cu = zeros(1,dim);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant   = 'no';
      options.ipopt.hessian_constant = 'no';
      options.ipopt.mu_strategy      = 'adaptive';
      options.ipopt.max_iter         = 400;
      options.ipopt.tol              = 1e-10;%
      options.ipopt.linear_solver    = 'mumps';

      % The callback functions.
      funcs.objective         = @(Z) self.NLP_target(Z);
      funcs.gradient          = @(Z) self.NLP_target_gradient(Z);

      funcs.constraints       = @(Z) self.NLP_constraints(Z);
      funcs.jacobian          = @(Z) self.NLP_constraints_jacobian(Z);
      funcs.jacobianstructure = @()  self.NLP_constraints_jacobian_pattern();

      if false
        options.ipopt.derivative_test = 'second-order';
        funcs.hessian                 = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure        = @() self.NLP_hessian_pattern();
      else
        options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      xguess  = self.x0+(self.x1-self.x0)*self.nodes;
      yguess  = ones(size(xguess));
      zguess  = self.L*ones(size(xguess));
      uaguess = zeros(1,self.N-self.nseg);
      ubguess = zeros(1,self.N-self.nseg);

      x0 = [ reshape( [ xguess, yguess, zguess ], xtot ,1 ); ...
             zeros( utot, 1 ) ];

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      T     = self.parameters();

      subplot( 3, 1, 1 );
      plot( nodes, X(1,:), 'Linewidth', 2 );
      title('x');

      subplot( 3, 1, 2 );
      plot( nodes, X(2,:), 'Linewidth', 2 );
      title('v');

      subplot( 3, 1, 3 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Unodes, UU(1,:), Unodes, UU(2,:), 'Linewidth', 2 );
      title('ua,ub');

    end

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( ~, ~, tL, tR, XL, XR, UC, ~ )
      xL = XL(1); yL = XL(2);
      xR = XR(1); yR = XR(2);
      u  = UC(1); v  = UC(2);
      x  = (xL+xR)/2;
      y  = (yL+yR)/2;
      L  = (tR-tL) * (x*v-y*u);
    end

    %
    function gradL = lagrange_gradient( ~, ~, tL, tR, XL, XR, UC, ~ )
      xL = XL(1); yL = XL(2);
      xR = XR(1); yR = XR(2);
      u  = UC(1); v  = UC(2);
      x  = (xL+xR)/2;
      y  = (yL+yR)/2;
      gradL = (tR-tL) * [ v/2, -u/2, 0, v/2, -u/2, 0, -y, x ];
    end

    %
    function hessL = lagrange_hessian( ~, ~, tL, tR, XL, XR, UC, ~ )
      hessL = sparse(8,8);
      tmp = (tR-tL)/2;
      hessL(1,8) = tmp;
      hessL(2,7) = -tmp;
      hessL(4,8) = tmp;
      hessL(5,7) = -tmp;
      hessL(7,2) = -tmp;
      hessL(7,5) = -tmp;
      hessL(8,1) = tmp;
      hessL(8,4) = tmp;
    end

    %
    function hessL = lagrange_hessian_pattern( ~ )
      hessL = sparse(8,8);
      hessL(1,8) = 1;
      hessL(2,7) = 1;
      hessL(4,8) = 1;
      hessL(5,7) = 1;
      hessL(7,2) = 1;
      hessL(7,5) = 1;
      hessL(8,1) = 1;
      hessL(8,4) = 1;
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( ~, tL, tR, XL, XR, ~ )
      M = 0;
    end

    %
    function gradM = mayer_gradient( ~, tL, tR, XL, XR, ~ )
      gradM = zeros(1,6);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( ~, tL, tR, XL, XR, ~ )
      hessM = sparse(6,6);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( ~ )
      hessM = sparse(6,6);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, nseg, tL, tR, XL, XR, UC, T )
      C = self.midpoint_ds( nseg, tL, tR, XL, XR, UC, T, @self.RHS );
    end
    %
    function CJ = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, T )
      CJ = self.midpoint_ds_jacobian( nseg, tL, tR, XL, XR, UC, T, @self.JAC );
    end
    %
    function CJ = ds_jacobian_pattern( self )
      CJ = self.midpoint_ds_jacobian_pattern( @self.JAC_pattern );
    end
    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, T, L )
      H = self.midpoint_ds_hessian( nseg, tL, tR, XL, XR, UC, T, L, @self.HESS );
    end
    %
    function H = ds_hessian_pattern( self )
      H = self.midpoint_ds_hessian_pattern( @self.HESS_pattern );
    end

    %              _   _                           _             _           _
    %  _ __   __ _| |_| |__     ___ ___  _ __  ___| |_ _ __ __ _(_)_ __  ___| |_
    % | '_ \ / _` | __| '_ \   / __/ _ \| '_ \/ __| __| '__/ _` | | '_ \/ __| __|
    % | |_) | (_| | |_| | | | | (_| (_) | | | \__ \ |_| | | (_| | | | | \__ \ |_
    % | .__/ \__,_|\__|_| |_|  \___\___/|_| |_|___/\__|_|  \__,_|_|_| |_|___/\__|
    % |_|
    %

    % Path constraints
    function C = pc( self, t, X, PARS )
      C = zeros(0,1);
    end

    function CJ = pc_jacobian( self, t, X, PARS )
      CJ = sparse(0,self.nx);
    end

    function CJ = pc_jacobian_pattern( self )
      CJ = sparse(0,self.nx);
    end

    function CH = pc_hessian( self, t, X, PARS, L )
      CH = sparse(self.nx,self.nx);
    end

    function CH = pc_hessian_pattern( self )
      CH = sparse(self.nx,self.nx);
    end

    %     _
    %  _ | |_  _ _ __  _ __
    % | || | || | '  \| '_ \
    %  \__/ \_,_|_|_|_| .__/
    %                 |_|
    %
    function JMP = jump( self, nsegL, t, XL, XR, P )
      JMP = self.jump_standard( nsegL, t, XL, XR, P );
    end

    %
    function JAC = jump_jacobian( self, nsegL, t, XL, XR, P )
      JAC = self.jump_standard_jacobian( nsegL, t, XL, XR, P );
    end

    %
    function JAC = jump_jacobian_pattern( self )
      JAC = self.jump_standard_jacobian_pattern();
    end

    %
    function H = jump_hessian( self, nsegL, t, XL, XR, P, L )
      H = self.jump_standard_hessian( nsegL, t, XL, XR, P, L );
    end

    %
    function H = jump_hessian_pattern( self )
      H = self.jump_standard_hessian_pattern( );
    end

    %  ___  ___
    % | _ )/ __|
    % | _ \ (__
    % |___/\___|
    %
    function bc = bc( self, tL, tR, XL, XR, ~ )
      xL = XL(1); yL = XL(2);
      xR = XR(1); yR = XR(2);
      bc = [ xL - self.x0; yL; xR - self.x1; yR ];
    end

    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, ~ )
      Jac = sparse(4,6);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,4) = 1;
      Jac(4,5) = 1;
    end

    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse(4,6);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,4) = 1;
      Jac(4,5) = 1;
    end

    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, L )
      Hess = sparse(6,6);
    end

    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = sparse(6,6);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    %
    function tvU = TVU( self, dt, UCL, UCR )
      tvU = self.TVU_zero( dt, UCL, UCR );
    end
    %
    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_zero_gradient( dt, UCL, UCR );
    end
    %
    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_zero_hessian( dt, UCL, UCR );
    end
    %
    function tvH = TVU_hessian_pattern( self )
      tvH = self.TVU_zero_hessian_pattern();
    end
  end
end
