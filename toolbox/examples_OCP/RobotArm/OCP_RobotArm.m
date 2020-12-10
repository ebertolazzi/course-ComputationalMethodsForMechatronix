%--------------------------------------------------------------------------%
%                                                                          %
%  Copyright (C) 2020                                                      %
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

classdef OCP_RobotArm < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    % PARAMETERS OF THE PROBLEM
    L % maximum lenght of the arm
    %
    epsi
  end

  methods ( Hidden = true )

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %
    function res = Iphi( self, rho )
      L   = self.L;
      res = ((L-rho)^3+rho^3)/3;
    end
    %%
    function res = Iphi_D( self, rho )
      L   = self.L;
      res = L * (2 * rho - L);
    end
    %%
    function res = Iphi_DD( self, rho )
      L   = self.L;
      res = 2*L;
    end
    %
    function res = Itheta( self, rho, phi )
      res = self.Iphi( rho )*sin(phi)^2;
    end
    %
    function res = Itheta_r( self, rho, phi )
      res = self.Iphi_D( rho )*sin(phi)^2;
    end
    %
    function res = Itheta_p( self, rho, phi )
      res = 2*self.Iphi( rho )*sin(phi)*cos(phi);
    end
    %
    function res = Itheta_rr( self, rho, phi )
      res = self.Iphi_DD( rho )*sin(phi)^2;
    end
    %
    function res = Itheta_rp( self, rho, phi )
      res = 2*self.Iphi_D( rho )*sin(phi)*cos(phi);
    end
    %
    function res = Itheta_pp( self, rho, phi )
      res = self.Iphi( rho )*(4*cos(phi)^2-2);
    end
    %   ___  ___  ___
    %  / _ \|   \| __|
    % | (_) | |) | _|
    %  \___/|___/|___|
    %
    %  rho, phi, theta, v_rho, v_phi, v_theta
    %
    %
    function RES = RHS( self, nseg, t, X, U, P )
      rho     = X(1);
      phi     = X(2);
      theta   = X(3);
      v_rho   = X(4);
      v_phi   = X(5);
      v_theta = X(6);
      u_rho   = U(1);
      u_phi   = U(2);
      u_theta = U(3);
      T       = P(1);
      L       = self.L;
      RES = [    ...
        T*v_rho;   ...
        T*v_phi;   ...
        T*v_theta; ...
        T*u_rho/L; ...
        T*u_phi/self.Iphi(rho); ...
        T*u_theta/self.Itheta(rho,phi); ...
      ];
    end
    %
    function RES = JAC( self, nseg, t, X, U, P )
      rho     = X(1);
      phi     = X(2);
      theta   = X(3);
      v_rho   = X(4);
      v_phi   = X(5);
      v_theta = X(6);
      u_rho   = U(1);
      u_phi   = U(2);
      u_theta = U(3);
      T       = P(1);
      L       = self.L;
      % 6 x, 3 u, 1 p
      RES = sparse( self.nx, self.nx+self.nu+self.np );
      % row 1
      RES(1,4) = T;
      RES(1,10) = v_rho;
      % row 2
      RES(2,5) = T;
      RES(2,10) = v_phi;
      % row 3
      RES(3,6) = T;
      RES(3,10) = v_theta;
      % row 4
      RES(4,7)  = T/L;
      RES(4,10) = u_rho/L;
      % row 5
      RES(5,1)  = -T*u_phi*self.Iphi_D(rho)/self.Iphi(rho)^2;
      RES(5,8)  = T/self.Iphi(rho);
      RES(5,10) = u_phi/self.Iphi(rho);
      % row 6
      RES(6,1)  = -T / self.Itheta(rho, phi)^2 * u_theta * self.Itheta_r(rho, phi);
      RES(6,2)  = -T / self.Itheta(rho, phi)^2 * u_theta * self.Itheta_p(rho, phi);
      RES(6,9)  = T / self.Itheta(rho, phi);
      RES(6,10) = u_theta/self.Itheta(rho,phi);
    end
    %
    function RES = JAC_pattern( self )
      RES = sparse( self.nx, self.nx+self.nu+self.np );
      RES(1,4)  = 1;
      RES(1,10) = 1;
      % row 2
      RES(2,5)  = 1;
      RES(2,10) = 1;
      % row 3
      RES(3,6)  = 1;
      RES(3,10) = 1;
      % row 4
      RES(4,7)  = 1;
      RES(4,10) = 1;
      % row 5
      RES(5,1)  = 1;
      RES(5,8)  = 1;
      RES(5,10) = 1;
      % row 6
      RES(6,1)  = 1;
      RES(6,2)  = 1;
      RES(6,9)  = 1;
      RES(6,10) = 1;
    end
  end

  methods

    function self = OCP_RobotArm( )
      nx  = 6; % number of states
      nu  = 3; % number of controls
      np  = 1; % number of free parameters
      nbc = 12; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      % setup the parameters of the model
      self.L     = 5;
      self.epsi  = 0; % no penalization of total variations of the control
    end

    function info = solve( self )

      % load discretization parmeters
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      N    = self.N;
      nseg = self.nseg;
      njmp = self.njmp;
      nbc  = self.nbc;

      L    = self.L;

      utot = (N-nseg)*nu;
      xtot = N*nx;

      % set the lower bound of the states (nx x N)
      x_lb = [ zeros(1,N); zeros(1,N); -pi*ones(1,N); ...
               -Inf*ones(1,N); -Inf*ones(1,N); -Inf*ones(1,N) ];

      % set the upper bound of the states (nx x N)
      x_ub = [ L*ones(1,N); pi*ones(1,N); pi*ones(1,N); ...
               Inf*ones(1,N); Inf*ones(1,N); Inf*ones(1,N) ];

      % set the lower bound of the controls (nu x (N-nseg))
      u_lb = -ones(3, N-nseg);

      % set the upper bound of the controls (nu x (N-nseg))
      u_ub = ones(3, N-nseg);

      % set the lower bound of the optimization parameters (1 x np)
      p_lb = 0;

      % set the upper bound of the optimization parameters (1 x np)
      p_ub = 10;

      options.lb = self.pack( x_lb, u_lb, p_lb ); % Lower bound on the variables.
      options.ub = self.pack( x_ub, u_ub, p_ub ); % Upper bound on the variables.

      % The constraint functions are bounded to zero
      dim = (N-nseg)*nx + (nseg-1)*njmp + nbc;
      options.cl = zeros(1,dim); %  constraints
      options.cu = zeros(1,dim);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant   = 'no';
      options.ipopt.hessian_constant = 'no';
      options.ipopt.mu_strategy      = 'adaptive';
      options.ipopt.max_iter         = 1000;
      options.ipopt.tol              = 1e-10;%
      options.ipopt.linear_solver    = 'mumps';

      % The callback functions.
      funcs.objective         = @(Z) self.NLP_target(Z);
      funcs.gradient          = @(Z) self.NLP_target_gradient(Z);

      funcs.constraints       = @(Z) self.NLP_constraints(Z);
      funcs.jacobian          = @(Z) self.NLP_constraints_jacobian(Z);
      funcs.jacobianstructure = @() self.NLP_constraints_jacobian_pattern();

      if false
        %options.ipopt.derivative_test = 'first-order';
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian           = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure  = @() self.NLP_hessian_pattern();
      else
        options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % GUESS SOLUTION
      Xguess = [ 4.5*ones(1,N); pi/4*ones(1,N); 2*pi/3*self.nodes.^2; ...
                 zeros(1,N); zeros(1,N); zeros(1,N) ];
      Uguess = [ zeros(1,N-nseg); zeros(1,N-nseg); zeros(1,N-nseg) ];
      Pguess = 1;

      x0 = self.pack( Xguess, Uguess, Pguess ); % Lower bound on the variables.
      
      tic
      [self.sol, info] = ipopt(x0,funcs,options);

      self.sol(end)

    end

    function plot( self )
      % example of pliot solution
      nodes = self.nodes;
      T = self.sol(end);

      [UU,Unodes] = self.controls_for_plot();
      subplot(3,1,1);
      plot( T*Unodes, UU(1,:), 'o-r', 'Linewidth', 2 );
      subplot(3,1,2);
      plot( T*Unodes, UU(2,:), 'o-g', 'Linewidth', 2 );
      subplot(3,1,3);
      plot( T*Unodes, UU(3,:), 'o-b', 'Linewidth', 2 );
    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( self, nseg, tL, tR, XL, XR, UC, PARS )
      L = self.lagrange_zero( nseg, tL, tR, XL, XR, UC, PARS );
    end
    %
    function gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, PARS )
      gradL = self.lagrange_zero_gradient( nseg, tL, tR, XL, XR, UC, PARS );
    end
    %
    function hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, PARS )
      hessL = self.lagrange_zero_hessian( nseg, tL, tR, XL, XR, UC, PARS );
    end
    %
    function hessL = lagrange_hessian_pattern( self )
      hessL = self.lagrange_zero_hessian_pattern();
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer( self, tL, tR, XL, XR, P )
      M = P(1);
    end
    %
    function gradM = mayer_gradient( self, tL, tR, XL, XR, P )
      gradM = [ zeros(1,12), 1];
    end
    %
    function hessM = mayer_hessian( self, tL, tR, XL, XR, P )
      hessM = sparse(13,13);
    end
    %
    function hessM = mayer_hessian_pattern( self )
      hessM = sparse(13,13);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, nseg, tL, tR, XL, XR, UC, PARS )
      C = self.midpoint_ds( nseg, tL, tR, XL, XR, UC, PARS, @self.RHS );
    end
    %
    function CJ = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, PARS )
      CJ = self.midpoint_ds_jacobian( nseg, tL, tR, XL, XR, UC, PARS, @self.JAC );
    end
    %
    function CJ = ds_jacobian_pattern( self )
      CJ = self.midpoint_ds_jacobian_pattern( @self.JAC_pattern );
    end
    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, PARS, L )
      H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, PARS, L );
    end
    %
    function H = ds_hessian_pattern( self )
      dim = 2*self.nx+self.nu+self.np;
      H   = sparse(ones(dim,dim));
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
    %
    function CJ = pc_jacobian( self, t, X, PARS )
      CJ = zeros(0,self.nx);
    end
    %
    function CJ = pc_jacobian_pattern( self )
      CJ = sparse(0,self.nx);
    end
    %
    function CH = pc_hessian( self, t, X, PARS, L )
      CH = sparse(self.nx,self.nx);
    end
    %
    function CH = pc_hessian_pattern( self )
      CH = sparse(self.nx,self.nx);
    end

    %     _
    %  _ | |_  _ _ __  _ __
    % | || | || | '  \| '_ \
    %  \__/ \_,_|_|_|_| .__/
    %                 |_|
    %
    function ODE = jump( self, nsegL, t, XL, XR, P )
      ODE = self.jump_standard( nsegL, t, XL, XR, P );
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
      H = self.jump_standard_hessian_pattern();
    end

    %  ___  ___
    % | _ )/ __|
    % | _ \ (__
    % |___/\___|
    %
    function bc = bc( self, tL, tR, XL, XR, P )
      rho_L     = XL(1);
      phi_L     = XL(2);
      theta_L   = XL(3);
      v_rho_L   = XL(4);
      v_phi_L   = XL(5);
      v_theta_L = XL(6);
      rho_R     = XR(1);
      phi_R     = XR(2);
      theta_R   = XR(3);
      v_rho_R   = XR(4);
      v_phi_R   = XR(5);
      v_theta_R = XR(6);
      bc = [ ...
        rho_L-4.5; ...
        rho_R-4.5; ...
        phi_L-pi/4; ...
        phi_R-pi/4; ...
        theta_L; ...
        theta_R-2*pi/3; ...
        v_rho_L; ...
        v_phi_L; ...
        v_theta_L; ...
        v_rho_R; ...
        v_phi_R; ...
        v_theta_R; ...
      ];
    end
    %
    function Jac = bc_jacobian( self, tL, tR, XL, XR, P )
      % rho, phi, theta, v_tho, v_phi, v_theta (LEFT)
      % rho, phi, theta, v_tho, v_phi, v_theta (RIGHT)
      % T
      dim = 2*self.nx+self.np;
      Jac = sparse(12,dim);
      Jac(1,1)   = 1;
      Jac(2,7)   = 1;
      Jac(3,2)   = 1;
      Jac(4,8)   = 1;
      Jac(5,3)   = 1;
      Jac(6,9)   = 1;
      Jac(7,4)   = 1;
      Jac(8,5)   = 1;
      Jac(9,6)   = 1;
      Jac(10,10) = 1;
      Jac(11,11) = 1;
      Jac(12,12) = 1;
    end
    %
    function Jac = bc_jacobian_pattern( self )
      % rho, phi, theta, v_tho, v_phi, v_theta (LEFT)
      % rho, phi, theta, v_tho, v_phi, v_theta (RIGHT)
      % T
      dim = 2*self.nx+self.np;
      Jac = sparse(12,dim);
      Jac(1,1)   = 1;
      Jac(2,7)   = 1;
      Jac(3,2)   = 1;
      Jac(4,8)   = 1;
      Jac(5,3)   = 1;
      Jac(6,9)   = 1;
      Jac(7,4)   = 1;
      Jac(8,5)   = 1;
      Jac(9,6)   = 1;
      Jac(10,10) = 1;
      Jac(11,11) = 1;
      Jac(12,12) = 1;
    end
    %
    function Hess = bc_hessian( self, tL, tR, XL, XR, P, L )
      dim = 2*self.nx+self.np;
      Hess = sparse(dim,dim);
    end
    %
    function Hess = bc_hessian_pattern( self )
      dim = 2*self.nx+self.np;
      Hess = sparse(dim,dim);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    % tvU
    function tvU = TVU( self, dt, UCL, UCR )
      tvU = self.TVU_reg( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvG = TVU_gradient( self, dt, UCL, UCR )
      tvG = self.TVU_reg_gradient( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvH = TVU_hessian( self, dt, UCL, UCR )
      tvH = self.TVU_reg_hessian( self.epsi*ones(self.nu,1), dt, UCL, UCR );
    end
    function tvH = TVU_hessian_pattern( self )
      tvH = self.TVU_reg_hessian_pattern();
    end

  end
end
