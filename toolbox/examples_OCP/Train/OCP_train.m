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

classdef OCP_train < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    x_i
    x_f
    v_i
    v_f
    alpha
    beta
    gm
    uepsi
    epsilon
    ss
    zz
    uaMax
    ubMax
    epsia % TV epsi ctrl a
    epsib % TV epsi ctrl b
  end

  methods

    function self = OCP_train( )
      nx  = 2; % number of states
      nu  = 2; % number of controls
      np  = 0; % number of free parameters
      nbc = 4; % number of boundary conditions
      self@OCP_NLP( nx, nu, np, nbc );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
      self.x_i     = 0;
      self.x_f     = 6;
      self.v_i     = 0;
      self.v_f     = 0;
      self.alpha   = 0.3;
      self.beta    = 0.14;
      self.gm      = 0.16;
      self.uepsi   = 1e-7;
      self.epsilon = 0.05;
      self.ss      = [ -2, 0, 2 ];
      self.zz      = [ 2, 4 ];
      self.uaMax   = 10;
      self.ubMax   = 2;
      self.epsia   = 1e-8;
      self.epsib   = 1e-8;
    end

    function info = solve( self )

      utot = (self.N-self.nseg)*self.nu;
      xtot = self.N*self.nx;

      xones = ones(1,xtot);
      umax  = [self.uaMax;self.ubMax]*ones(1,self.N-self.nseg);

      options.lb = [ -xones*Inf, zeros(1,utot) ];        % Lower bound on the variables.
      options.ub = [  xones*Inf, reshape(umax,1,utot) ]; % Upper bound on the variables.

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
      funcs.jacobianstructure = @() self.NLP_constraints_jacobian_pattern();

      if true
        %options.ipopt.derivative_test = 'second-order';
        funcs.hessian           = @( Z, sigma, lambda ) self.NLP_hessian( Z, sigma, lambda );
        funcs.hessianstructure  = @() self.NLP_hessian_pattern();
      else
        %options.ipopt.derivative_test            = 'first-order';
        options.ipopt.hessian_approximation      = 'limited-memory';
        options.ipopt.limited_memory_update_type = 'bfgs'; % {bfgs}, sr1 = 6; % {6}
      end

      % Run IPOPT.
      xguess  = self.x_i+(self.x_f-self.x_i)*self.nodes./4.8;
      vguess  = ones(1,self.N);
      uaguess = zeros(1,self.N-self.nseg);
      ubguess = zeros(1,self.N-self.nseg);

      x0 = [ reshape( [ xguess,   vguess], self.N*self.nx ,1 ); ...
             reshape( [ uaguess, ubguess], (self.N-self.nseg)*self.nu ,1 ) ];

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
    function res = hfun( self, x )
      res = 0;
      for j=1:2
        res = res + (self.ss(j+1)-self.ss(j))*atan((x-self.zz(j))/self.epsilon);
      end
      res = res / pi;
    end

    %
    function res = hfun_D( self, x )
      res = 0;
      for j=1:2
        res = res + (self.ss(j+1)-self.ss(j))/(1+((x-self.zz(j))/self.epsilon)^2);
      end
      res = res / pi / self.epsilon;
    end

    %
    function res = hfun_DD( self, x )
      res = 0;
      for j=1:2
        dz  = x-self.zz(j);
        dz2 = (dz/self.epsilon)^2;
        res = res + (self.ss(j+1)-self.ss(j))*dz/(1+dz2)^2;
      end
      res = -2*res / pi / self.epsilon^3;
    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( ~, ~, tL, tR, XL, XR, UC, ~ )
      ua = UC(1); ub = UC(2);
      xL = XL(1); vL = XL(2);
      xR = XR(1); vR = XR(2);
      v  = (vL+vR)/2;
      L  = (tR-tL) * ua * v;
    end

    %
    function gradL = lagrange_gradient( ~, ~, tL, tR, XL, XR, UC, ~ )
      ua = UC(1); ub = UC(2);
      xL = XL(1); vL = XL(2);
      xR = XR(1); vR = XR(2);
      v  = (vL+vR)/2;
      gradL = (tR-tL) * [ 0, ua/2, 0, ua/2,  v, 0];
    end

    %
    function hessL = lagrange_hessian( ~, ~, tL, tR, XL, XR, UC, ~ )
      hessL = sparse(6,6);
      tmp   = (tR-tL) *0.5;
      hessL(2,5) = tmp;
      hessL(4,5) = tmp;
      hessL(5,2) = tmp;
      hessL(5,4) = tmp;
    end

    %
    function hessL = lagrange_hessian_pattern( ~ )
      hessL = sparse(6,6);
      hessL(2,5) = 1;
      hessL(4,5) = 1;
      hessL(5,2) = 1;
      hessL(5,4) = 1;
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
      xL = XL(1); vL = XL(2);
      xR = XR(1); vR = XR(2);
      ua = UC(1); ub = UC(2);
      % ----------
      DT = tR - tL;
      xM = (xR+xL)/2;
      vM = (vR+vL)/2;
      % ----------
      C    = zeros(2,1);
      acc  = self.hfun(xM) - ( self.alpha + self.beta * vM + self.gm * vM^2 ); 
      C(1) = (xR - xL)/DT - vM;
      C(2) = (vR - vL)/DT - acc - ua + ub;
    end

    %
    function JAC = ds_jacobian( self, ~, tL, tR, XL, XR, UC, ~ )
      xL = XL(1); vL = XL(2);
      xR = XR(1); vR = XR(2);
      ua = UC(1); ub = UC(2);
      % ----------
      DT = tR - tL;
      xM = (xR+xL)/2;
      vM = (vR+vL)/2;
      % ----------
      tmp_x = -0.5 * self.hfun_D(xM);
      tmp_v =  0.5 * self.beta + self.gm * vM;
      JAC = sparse([ -1/DT,       -0.5,  1/DT,       -0.5,  0, 0; ...
                     tmp_x, tmp_v-1/DT, tmp_x, tmp_v+1/DT, -1, 1 ] );
    end

    %
    function JAC = ds_jacobian_pattern( self )
      JAC = sparse([ 1, 1, 1, 1, 0, 0; ...
                     1, 1, 1, 1, 1, 1 ] );
    end

    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, P, L )
      if false
        H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, P, L );
      else
        xL = XL(1); vL = XL(2);
        xR = XR(1); vR = XR(2);
        ua = UC(1); ub = UC(2);
        % ----------
        DT = tR - tL;
        xM = (xR+xL)/2;
        vM = (vR+vL)/2;
        % ----------
        nx = 2;
        nu = 2;
        tmp_xx = -0.25 * self.hfun_DD(xM);
        tmp_vv =  0.5 * self.gm;
        Ahess  = [ tmp_xx, 0; 0, tmp_vv ];
        H      = sparse( L(2) * [ Ahess, Ahess, zeros(nx,nu); ...
                                  Ahess, Ahess, zeros(nx,nu); ...
                                  zeros(nu,2*nx+nu) ] );
      end
    end

    %
    function H = ds_hessian_pattern( self )
      Ahess  = speye(2);
      H      = sparse( [ Ahess, Ahess, zeros(self.nx,self.nu); ...
                         Ahess, Ahess, zeros(self.nx,self.nu); ...
                         zeros(self.nu,2*self.nx+self.nu) ] );
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
      xL = XL(1); vL = XL(2);
      xR = XR(1); vR = XR(2);
      bc = [ xL - self.x_i; ...
             xR - self.x_f; ...
             vL - self.v_i; ...
             vR - self.v_f ];
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
    %
    function tvU = TVU( self, dt, UCL, UCR )
      %tvU = self.TVU_zero( dt, UCL, UCR );
      tvU = self.TVU_reg( [self.epsia;self.epsib], dt, UCL, UCR );
    end

    function tvG = TVU_gradient( self, dt, UCL, UCR )
      %tvG = self.TVU_zero_gradient( dt, UCL, UCR );
      tvG = self.TVU_reg_gradient( [self.epsia;self.epsib], dt, UCL, UCR );
    end

    function tvH = TVU_hessian( self, dt, UCL, UCR )
      %tvH = self.TVU_zero_hessian( dt, UCL, UCR );
      tvH = self.TVU_reg_hessian( [self.epsia;self.epsib], dt, UCL, UCR );
    end

    function tvH = TVU_hessian_pattern( self )
      %tvH = self.TVU_zero_hessian_pattern();
      tvH = self.TVU_reg_hessian_pattern();
    end
  end
end
