%-------------------------------------------------------------------------%
%                                                                         %
%  Copyright (C) 2018                                                     %
%                                                                         %
%         , __                 , __                                       %
%        /|/  \               /|/  \                                      %
%         | __/ _   ,_         | __/ _   ,_                               %
%         |   \|/  /  |  |   | |   \|/  /  |  |   |                       %
%         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                      %
%                           /|                   /|                       %
%                           \|                   \|                       %
%                                                                         %
%      Enrico Bertolazzi                                                  %
%      Dipartimento di Ingegneria Industriale                             %
%      Universita` degli Studi di Trento                                  %
%      email: enrico.bertolazzi@unitn.it                                  %
%                                                                         %
%-------------------------------------------------------------------------%

classdef OCP_Stirred_Tank < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x1_i
    x1_f
    x2_i
    x2_f
    u_f
    a
    theta
    k
    En
    Tc
    Tf
  end

  methods

    function self = OCP_Stirred_Tank( )
      nx  = 2; % number of states
      nu  = 1; % number of controls
      np  = 0; % number of free parameters
      nbc = 4; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.a     = 0.000195;%*600;
      self.theta = 20;
      self.k     = 300;
      self.En    = 5;
      self.Tc    = 0.38158;
      self.Tf    = 0.3947;
      self.x1_i  = 0.98;
      self.x2_i  = 0.39;
      self.x1_f  = 0.26;
      self.x2_f  = 0.65;
      self.u_f   = 0.76;
    end

    function info = solve( self )

      utot = (self.N-self.nseg)*self.nu;
      xtot = self.N*self.nx;

      xones      = ones(1,xtot);
      uones      = ones(1,utot);
      options.lb = [ 0*xones, 0*uones ];        % Lower bound on the variables.
      options.ub = [   xones, 2*uones ]; % Upper bound on the variables.

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
      ss      = (self.nodes-self.nodes(1))./(self.nodes(end)-self.nodes(1));
      x1guess = self.x1_i+(self.x1_f-self.x1_i)*ss;
      x2guess = self.x2_i+(self.x2_f-self.x2_i)*ss;
      uguess  = 0*ones(1,self.N-self.nseg);

      x0 = [ reshape( [ x1guess, x2guess], self.N*self.nx ,1 ); ...
             reshape( uguess, (self.N-self.nseg)*self.nu ,1 ) ];

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      nodes = self.nodes;
      X     = self.states();
      %T     = self.parameters();

      subplot( 2, 1, 1 );
      plot( nodes, X(1,:), nodes, X(2,:), 'Linewidth', 2 );
      title('x1,x2');

      subplot( 2, 1, 2 );
      [UU,Unodes] = self.controls_for_plot();
      plot( Unodes, UU(1,:), 'Linewidth', 2 );
      title('u');

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
    function L = lagrange( self, ~, tL, tR, XL, XR, UC, ~ )
      u   = UC(1);
      x1L = XL(1); x2L = XL(2);
      x1R = XR(1); x2R = XR(2);
      x1m = (x1L+x1R)/2;
      x2m = (x2L+x2R)/2;
      dx1 = x1m - self.x1_f;
      dx2 = x2m - self.x2_f;
      du  = u - self.u_f;
      L   = (tR-tL) * ( dx1^2 + dx2^2 + du^2 );
    end

    %
    function gradL = lagrange_gradient( self, ~, tL, tR, XL, XR, UC, ~ )
      u   = UC(1);
      x1L = XL(1); x2L = XL(2);
      x1R = XR(1); x2R = XR(2);
      x1m = (x1L+x1R)/2;
      x2m = (x2L+x2R)/2;
      dx1 = x1m - self.x1_f;
      dx2 = x2m - self.x2_f;
      du  = u - self.u_f;
      gradL = (tR-tL) * [ dx1, dx2, dx1, dx2, 2*du];
    end

    %
    function hessL = lagrange_hessian( ~, ~, tL, tR, XL, XR, UC, ~ )
      hessL = sparse(0.5*(tR-tL)*[ ...
        1, 0, 1, 0, 0; ...
        0, 1, 0, 1, 0; ...
        1, 0, 1, 0, 0; ...
        0, 1, 0, 1, 0; ...
        0, 0, 0, 0, 4; ...
      ]);
    end

    %
    function hessL = lagrange_hessian_pattern( ~ )
      hessL = sparse([ ...
        1, 0, 1, 0, 0; ...
        0, 1, 0, 1, 0; ...
        1, 0, 1, 0, 0; ...
        0, 1, 0, 1, 0; ...
        0, 0, 0, 0, 1; ...
      ]);
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
      gradM = zeros(1,4);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian( ~, tL, tR, XL, XR, ~ )
      hessM = sparse(4,4);
    end

    % [ M, gradM, hessianM ]
    function hessM = mayer_hessian_pattern( ~ )
      hessM = sparse(4,4);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, ~, tL, tR, XL, XR, UC, ~ )
      x1L = XL(1); x2L = XL(2);
      x1R = XR(1); x2R = XR(2);
      x1m = (x1L+x1R)/2;
      x2m = (x2L+x2R)/2;
      u   = UC(1);
      % ----------
      DT  = tR - tL;
      % ----------
      C     = zeros(2,1);
      react = self.k * x1m * exp( - self.En/x2m );
      ctrl  = self.a * u * ( x2m - self.Tc );
      C(1)  = (x1R - x1L)/DT - (1-x1m)/self.theta + react;
      C(2)  = (x2R - x2L)/DT - (self.Tf-x2m)/self.theta - react + ctrl;
    end

    %
    function JAC = ds_jacobian( self, ~, tL, tR, XL, XR, UC, ~ )
      x1L = XL(1); x2L = XL(2);
      x1R = XR(1); x2R = XR(2);
      x1m = (x1L+x1R)/2;
      x2m = (x2L+x2R)/2;
      u   = UC(1);
      % ----------
      DT = tR - tL;
      % ----------
      exp1  = self.k * exp( - self.En/x2m );
      %react = x1m * exp1;
      % -----------
      exp1_2  = 0.5*(self.En/x2m^2)*exp1;
      react_1 = 0.5*exp1;
      react_2 = x1m * exp1_2;
      ctrl_2  = 0.5*self.a * u;
      ctrl_u  = self.a * ( x2m - self.Tc );
      tmp     = 0.5/self.theta;

      JAC = sparse(2,5);

      JAC(1,1) = -1/DT+tmp+react_1;
      JAC(1,2) = react_2;
      JAC(1,3) = 1/DT+tmp+react_1;
      JAC(1,4) = react_2;

      JAC(2,1) = -react_1;
      JAC(2,2) = -1/DT+tmp-react_2+ctrl_2;
      JAC(2,3) = -react_1;
      JAC(2,4) = 1/DT+tmp-react_2+ctrl_2;
      JAC(2,5) = ctrl_u;
    end

    %
    function JAC = ds_jacobian_pattern( self )
      JAC = sparse(2,5);
      JAC(1,1) = 1;
      JAC(1,2) = 1;
      JAC(1,3) = 1;
      JAC(1,4) = 1;
      JAC(2,1) = 1;
      JAC(2,2) = 1;
      JAC(2,3) = 1;
      JAC(2,4) = 1;
      JAC(2,5) = 1;
    end

    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, P, L )
      H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, P, L );
    end

    %
    function H = ds_hessian_pattern( self )
      H = sparse(ones(5,5));
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
      x1L = XL(1); x2L = XL(2);
      x1R = XR(1); x2R = XR(2);
      bc = [ x1L - self.x1_i; ...
             x1R - self.x1_f; ...
             x2L - self.x2_i; ...
             x2R - self.x2_f ];
    end

    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, ~ )
      Jac = sparse(4,4);
      Jac(1,1) = 1;
      Jac(2,3) = 1;
      Jac(3,2) = 1;
      Jac(4,4) = 1;
    end

    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = sparse(4,4);
      Jac(1,1) = 1;
      Jac(2,3) = 1;
      Jac(3,2) = 1;
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
    % tvU
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
