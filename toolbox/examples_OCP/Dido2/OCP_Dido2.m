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

classdef OCP_Dido2 < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x0;
    x1;
    L;
    epsi;
  end

  methods ( Hidden = true )
    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    function RES = RHS( self, nseg, t, X, U, T )
      theta = U(1);
      RES = [ cos(theta); sin(theta) ];
    end
    %
    function RES = JAC( self, nseg, t, X, U, T )
      RES      = sparse( 2, 3 );
      theta    = U(1);
      RES(1,3) = -sin(theta);
      RES(2,3) = cos(theta);
    end
    %
    function RES = JAC_pattern( self )
      RES      = sparse( 2, 3 );
      RES(1,3) = 1;
      RES(2,3) = 1;
    end
    %
    function H = HESS( self, nseg, t, X, U, T, L )
      H      = sparse(3,3);
      theta  = U(1);
      H(3,3) = -(L(1)*cos(theta)+L(2)*sin(theta));
    end
    %
    function H = HESS_pattern( self )
      H      = sparse(3,3);
      H(3,3) = 1;
    end
    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lag( ~, ~, t, X, U, ~ )
      x     = X(1);
      y     = X(2);
      theta = U(1);
      u     = cos(theta);
      v     = sin(theta);
      L     = x*v-y*u;
    end

    %
    function gradL = lag_gradient( ~, ~, t, X, U, ~ )
      x     = X(1);
      y     = X(2);
      theta = U(1);
      u     = cos(theta);
      v     = sin(theta);
      gradL = [ v, -u, u*x+v*y ];
    end

    %
    function hessL = lag_hessian( ~, ~, t, X, U, ~ )
      hessL = sparse(3,3);
      x     = X(1);
      y     = X(2);
      theta = U(1);
      u     = cos(theta);
      v     = sin(theta);
      hessL(1,3) = u;
      hessL(2,3) = v;
      hessL(3,1) = u;
      hessL(3,2) = v;
      hessL(3,3) = y*u-x*v;
    end

    %
    function hessL = lag_hessian_pattern( ~ )
      hessL = sparse(3,3);
      hessL(1,3) = 1;
      hessL(3,1) = 1;
      hessL(2,3) = 1;
      hessL(3,2) = 1;
      hessL(3,3) = 1;
    end
  end

  methods

    function self = OCP_Dido2( L, x0, x1 )
      nx  = 2; % number of states
      nu  = 1; % number of controls
      np  = 0; % number of free parameters
      nbc = 4; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
      self.x0   = x0;
      self.x1   = x1;
      self.L    = L;
      self.epsi = 1e-8;
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
    end

    function info = solve( self )

      utot = (self.N-self.nseg)*self.nu;
      xtot = self.N*self.nx;

      xones = ones(1,xtot);
      options.lb = [ -xones*Inf, -200*pi*ones(1,utot) ]; % Lower bound on the variables.
      options.ub = [  xones*Inf,  200*pi*ones(1,utot) ]; % Upper bound on the variables.

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

      if true
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian                 = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure        = @() self.NLP_hessian_pattern();
      else
        %options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      zeta   = (0:self.N-1)/(self.N-1);
      xguess = self.x0+(self.x1-self.x0).*zeta;
      yguess = 4*self.L*zeta.*(1-zeta);
      uguess = atan2(diff(yguess),diff(xguess));
      %plot(xguess,yguess);
      %stop;

      x0 = [ reshape( [ xguess, yguess ], xtot ,1 ); uguess(:) ];

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      T     = self.parameters();

      subplot( 2, 1, 1 );
      plot( X(1,:), X(2,:), '-o', 'Linewidth', 2 );
      axis equal;
      title('x/y');

      subplot( 2, 1, 2 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Unodes, UU(1,:), 'Linewidth', 2 );
      title('theta');

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
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, P )
      L = self.midpoint_lagrange(nseg, tL, tR, XL, XR, UC, P, @self.lag );
    end
    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, P )
      gradL = self.midpoint_lagrange_gradient( ...
        nseg, tL, tR, XL, XR, UC, P, @self.lag_gradient ...
      );
    end
    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, P )
      hessL = self.midpoint_lagrange_hessian( ...
        nseg, tL, tR, XL, XR, UC, P, @self.lag_hessian ...
      );
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      hessL = self.midpoint_lagrange_hessian_pattern( ...
        @self.lag_hessian_pattern ...
      );
    end
    
    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( self, tL, tR, XL, XR, P )
      M = self.mayer_zero( tL, tR, XL, XR, P );
    end
    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, P )
      gradM = self.mayer_zero_gradient( tL, tR, XL, XR, P );
    end
    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( self, tL, tR, XL, XR, P )
      hessM = self.mayer_zero_hessian( tL, tR, XL, XR, P );
    end
    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( self )
      hessM = self.mayer_zero_hessian_pattern();
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
      Jac = sparse(4,4);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,4) = 1;
    end
    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse(4,4);
      Jac(1,1) = 1;
      Jac(2,2) = 1;
      Jac(3,3) = 1;
      Jac(4,4) = 1;
    end
    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, L )
      Hess = sparse(4,4);
    end
    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = sparse(4,4);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    %
    function tvU = TVU( self, dt, UCL, UCR )
      tvU = self.TVU_zero( dt, UCL, UCR );
      %tvU = self.TVU_reg( self.epsi, dt, UCL, UCR );
    end

    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_zero_gradient( dt, UCL, UCR );
      %tvG = self.TVU_reg_gradient( self.epsi, dt, UCL, UCR );
    end

    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_zero_hessian( dt, UCL, UCR );
      %tvH = self.TVU_reg_hessian( self.epsi, dt, UCL, UCR );
    end

    function tvH = TVU_hessian_pattern( self )
      tvH = self.TVU_zero_hessian_pattern();
      %tvH = self.TVU_reg_hessian_pattern();
    end
  end
end
