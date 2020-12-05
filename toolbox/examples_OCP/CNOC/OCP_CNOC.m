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

classdef OCP_CNOC < OCP_NLP

  properties (SetAccess = private, Hidden = true)
    kappa_data
    kappa_break
    nominalFeed
    deltaFeed
  end

  methods

    function self = OCP_CNOC( )
      nx   = 7; % number of states
      nu   = 2; % number of controls
      np   = 0; % number of free parameters
      nbc  = 12; % number of boundary conditions
      npth = 0;
      njmp = 8;
      self@OCP_NLP( nx, nu, np, nbc, npth, njmp );
    end

    function setup( self, nodes )
      setup@OCP_NLP( self, nodes );
    end

    function info = solve( self )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      nbc  = self.nbc;
      %ntol = self.ntol;

      XLB = reshape( [ self.x1_min; self.x2_min; self.x3_min; self.x4_min ] * ones(1,N), 1, N*nx );
      XUB = reshape( [ self.x1_max; self.x2_max; self.x3_max; self.x4_max ] * ones(1,N), 1, N*nx );

      ULB = self.u_min * ones(1,N-self.nseg);
      UUB = self.u_max * ones(1,N-self.nseg);

      options.lb = [ XLB, ULB ];
      options.ub = [ XUB, UUB ];

      options.lb = [-ones(1,N*nx)*Inf, zeros( 1, (N-1)*nu ) ];  % Lower bound on the variables.
      options.ub = [ ones(1,N*nx)*Inf, reshape( [uaMax;ubMax]*ones(1,N-1), 1, nu*(N-1)) ];  % Upper bound on the variables.

      % The constraint functions are bounded to zero
      options.cl = zeros(1,(N-1)*nx+nbc); %  constraints
      options.cu = zeros(1,(N-1)*nx+nbc);

      % Set the IPOPT options.
      options.ipopt.jac_d_constant   = 'no';
      options.ipopt.hessian_constant = 'no';
      options.ipopt.mu_strategy      = 'adaptive';
      options.ipopt.max_iter         = 400;
      options.ipopt.tol              = 1e-10;
      options.ipopt.linear_solver    = 'mumps';
 
      % The callback functions.
      funcs.objective         = @(Z) self.NLP_target(Z);
      funcs.gradient          = @(Z) self.NLP_target_gradient(Z);

      funcs.constraints       = @(Z) self.NLP_constraints(Z);
      funcs.jacobian          = @(Z) self.NLP_constraints_jacobian(Z);
      funcs.jacobianstructure = @()  self.NLP_constraints_jacobian_pattern();

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
      vguess  = ones(1,N);
      uaguess = zeros(1,N-1);
      ubguess = zeros(1,N-1);

      x0 = [ reshape( [ xguess,   vguess], 2*N ,1 ); ...
             reshape( [ uaguess, ubguess], 2*(N-1) ,1 ) ];

      tic
      [self.sol, info] = ipopt(x0,funcs,options);
      elapsed = toc;

    end

    function plot( self )
      N     = self.N;
      x     = self.sol(1:2:2*N);
      v     = self.sol(2:2:2*N);
      ua    = self.sol(2*N+(1:2:2*N-2));
      ub    = self.sol(2*N+(2:2:2*N-2));
      nodes = self.nodes;

      subplot( 3, 1, 1 );
      plot( nodes, x, 'Linewidth', 2 );
      title('x');

      subplot( 3, 1, 2 );
      plot( nodes, v, 'Linewidth', 2 );
      title('v');

      subplot( 3, 1, 3 );
      plot( nodes(1:end-1), ua, nodes(1:end-1), ub, 'Linewidth', 2 );
      title('ua,ub');

    end

    %                      __              _   _
    %  _  _ ___ ___ _ _   / _|_  _ _ _  __| |_(_)___ _ _  ___
    % | || (_-</ -_) '_| |  _| || | ' \/ _|  _| / _ \ ' \(_-<
    %  \_,_/__/\___|_|   |_|  \_,_|_||_\__|\__|_\___/_||_/__/
    %
    function res = kappa( self, s )
      res = 0;
      for k=1:length(self.kappa_break)
        if s < self.kappa_break(k)
          res = self.kappa_data(k);
          break;
        end
      end
    end

    function res = kappa_D( self, s )
      res = 0;
    end

    function res = kappa_DD( self, s )
      res = 0;
    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange( ~, nseg, tL, tR, XL, XR, UC, ~ )
      sL   = XL(1); sR   = XR(1);
      nL   = XL(2); nR   = XR(2);
      vsL  = XL(3); vsR  = XR(3);
      vnL  = XL(4); vnR  = XR(4);
      asL  = XL(5); asR  = XR(5);
      anL  = XL(6); anR  = XR(6);
      coVL = XL(7); coVR = XR(7);

      js = UC(1);
      jn = UC(2);

      coV = (coVL+coVR)/2;
      vs  = (vsL+vsR)/2;
      vn  = (vnL+vnR)/2;

      L  = (tR-tL) * (coV*((hypot(vs,vn)-self.nominalFeed)/self.deltaFeed)^2);
    end

    %
    function gradL = lagrange_gradient( ~, nseg, tL, tR, XL, XR, UC, ~ )
      sL   = XL(1); sR   = XR(1);
      nL   = XL(2); nR   = XR(2);
      vsL  = XL(3); vsR  = XR(3);
      vnL  = XL(4); vnR  = XR(4);
      anL  = XL(6); anR  = XR(6);
      coVL = XL(7); coVR = XR(7);

      js = UC(1);
      jn = UC(2);

      coV = (coVL+coVR)/2;
      vs  = (vsL+vsR)/2;
      vn  = (vnL+vnR)/2;

      hp  = hypot(vs,vn);

      tmp = coV*((1-self.nominalFeed/hp)/self.deltaFeed^2);
      t1  = tmp*vs;
      t2  = tmp*vn;
      t3  = ((hp-self.nominalFeed)/self.deltaFeed)^2;

      gradL = (tR-tL) * [ 0, 0, t1, t2, 0, 0, t3, 0, 0, t1, t2, 0, 0, t3, 0, 0 ];
    end

    %
    function hessL = lagrange_hessian( ~, nseg, tL, tR, XL, XR, UC, ~ )
      sL   = XL(1); sR   = XR(1);
      nL   = XL(2); nR   = XR(2);
      vsL  = XL(3); vsR  = XR(3);
      vnL  = XL(4); vnR  = XR(4);
      asL  = XL(5); asR  = XR(5);
      anL  = XL(6); anR  = XR(6);
      coVL = XL(7); coVR = XR(7);

      js = UC(1);
      jn = UC(2);

      coV = (coVL+coVR)/2;
      vs  = (vsL+vsR)/2;
      vn  = (vnL+vnR)/2;

      tmp   = hypot(vs,vn);
      tvsvs = 0.5*(1-vn^2*nominalFeed/tmp^3)*coV/self.deltaFeed^2;
      tvnvn = 0.5*(1-vs^2*nominalFeed/tmp^3)*coV/self.deltaFeed^2;
      tvsvn = 0.5*self.nominalFeed*coV*vs*vn/(self.deltaFeed^2*tmp^3);
      tC    = 0.5*(1-nominalFeed/tmp)/self.deltaFeed^2;
      tvsC  = tC*vs;
      tvnC  = tC*vn;

      hessL = sparse( 16, 16 );
      hessL(3,3) = (tR-tL) * tvsvs;
      hessL(3,4) = (tR-tL) * tvsvn; hessL(4,3) = hessL(3,4);
      hessL(4,4) = (tR-tL) * tvnvn;

      hessL(7,3) = (tR-tL) * tvsC; hessL(3,7) = hessL(7,3);
      hessL(7,4) = (tR-tL) * tvnC; hessL(4,7) = hessL(7,4);

    end

    %
    function hessL = lagrange_hessian_pattern( ~ )
      hessL = sparse( 16, 16 );
      hessL(3,3) = 1;
      hessL(3,4) = 1;
      hessL(4,3) = 1;
      hessL(4,4) = 1;
      hessL(7,3) = 1;
      hessL(3,7) = 1;
      hessL(7,4) = 1;
      hessL(4,7) = 1;
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
      gradM = sparse(1,14);
    end
    %
    function hessM = mayer_hessian( ~, tL, tR, XL, XR, ~ )
      hessM = sparse(14,14);
    end
    %
    function hessM = mayer_hessian_pattern( ~ )
      hessM = sparse(14,14);
    end

    %   ___  ___  ___   _____   _   ___
    %  / _ \|   \| __| / /   \ /_\ | __|
    % | (_) | |) | _| / /| |) / _ \| _|
    %  \___/|___/|___/_/ |___/_/ \_\___|
    %
    function C = ds( self, nseg, tL, tR, XL, XR, UC, ~ )
      sL   = XL(1); sR   = XR(1);
      nL   = XL(2); nR   = XR(2);
      vsL  = XL(3); vsR  = XR(3);
      vnL  = XL(4); vnR  = XR(4);
      asL  = XL(5); asR  = XR(5);
      anL  = XL(6); anR  = XR(6);
      coVL = XL(7); coVR = XR(7);

      js = UC(1);
      jn = UC(2);

      DT  = tR - tL;
      coV = (coVL+coVR)/2;
      s   = (sL+sR)/2;
      n   = (nL+nR)/2;
      vs  = (vsL+vsR)/2;
      vn  = (vnL+vnR)/2;
      ks  = self.kappa(nseg);

      % ----------
      C    = zeros(7,1);
      C(1) = (sR - sL)/DT - vs*coV/(1-n*ks);
      C(2) = (nR - nL)/DT - vn*coV;
      C(3) = (vsR - vsL)/DT - as - vn*kn * (sR - sL)/DT;
      C(4) = (vnR - vnL)/DT - an + vs*kn * (sR - sL)/DT;
      C(5) = (asR - asL)/DT - js - an*kn * (sR - sL)/DT;
      C(6) = (anR - anL)/DT - jn + as*kn * (sR - sL)/DT;
      C(7) = (coVR - coVL)/DT;
    end
    %
    function JAC = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, ~ )
      sL   = XL(1); sR   = XR(1);
      nL   = XL(2); nR   = XR(2);
      vsL  = XL(3); vsR  = XR(3);
      vnL  = XL(4); vnR  = XR(4);
      asL  = XL(5); asR  = XR(5);
      anL  = XL(6); anR  = XR(6);
      coVL = XL(7); coVR = XR(7);

      js = UC(1);
      jn = UC(2);

      DT = tR - tL;
      coV = (coVL+coVR)/2;
      s   = (sL+sR)/2;
      n   = (nL+nR)/2;
      vn  = (vnL+vnR)/2;
      ks  = self.kappa(nseg);

      % ----------
      JAC       = sparse(7,16);
      JAC(1,1)  = -1/DT;
      JAC(1,8)  =  1/DT;
      JAC(1,2)  = -0.5*vs*coV*ks/(1-n*ks)^2;
      JAC(1,3)  = -0.5*coV/(1-n*ks);
      JAC(1,7)  = -0.5*vs/(1-n*ks);
      JAC(1,9)  = JAC(1,2);
      JAC(1,10) = JAC(1,3);
      JAC(1,14) = JAC(1,7);

      JAC(2,2)  = -1/DT;
      JAC(2,9)  =  1/DT;
      JAC(2,4)  = -0.5*coV;
      JAC(2,7)  = -0.5*vn;
      JAC(2,11) = JAC(2,4);
      JAC(2,14) = JAC(2,7);

      JAC(3,3)  = -1/DT;
      JAC(3,10) =  1/DT;
      JAC(3,5)  = -0.5;
      JAC(3,12) = -0.5;
      JAC(3,4)  = -0.5 * kn * (sR - sL)/DT;
      JAC(3,11) = JAC(3,4);
      JAC(3,1)  = vn*kn/DT;
      JAC(3,8)  = -JAC(3,1);

      JAC(4,4)  = -1/DT;
      JAC(4,11) =  1/DT;
      JAC(4,6)  = -0.5;
      JAC(4,13) = -0.5;
      JAC(4,3)  = 0.5 * kn * (sR - sL)/DT;
      JAC(4,10) = JAC(4,3);
      JAC(4,1)  = vs*kn/DT;
      JAC(4,8)  = -JAC(4,1);

      JAC(5,5)  = -1/DT;
      JAC(5,12) =  1/DT;
      JAC(5,15) = -1;
      JAC(5,6)  = 0.5 * kn * (sR - sL)/DT;
      JAC(5,13) = JAC(5,6);
      JAC(5,1)  = an*kn/DT;
      JAC(5,8)  = -JAC(5,1);

      JAC(6,6)  = -1/DT;
      JAC(6,13) =  1/DT;
      JAC(6,16) = -1;
      JAC(6,5)  = 0.5 * kn * (sR - sL)/DT;
      JAC(6,13) = JAC(6,5);
      JAC(6,1)  = as*kn/DT;
      JAC(6,8)  = -JAC(6,1);

      JAC(7,7)  = -1/DT;
      JAC(7,14) =  1/DT;

    end
    %
    function JAC = ds_jacobian_pattern( self )
      JAC       = sparse(7,16);
      JAC(1,1)  = 1;
      JAC(1,8)  = 1;
      JAC(1,2)  = 1;
      JAC(1,3)  = 1;
      JAC(1,7)  = 1;
      JAC(1,9)  = 1;
      JAC(1,10) = 1;
      JAC(1,14) = 1;

      JAC(2,2)  = 1;
      JAC(2,9)  = 1;
      JAC(2,4)  = 1;
      JAC(2,7)  = 1;
      JAC(2,11) = 1;
      JAC(2,14) = 1;

      JAC(3,3)  = 1;
      JAC(3,10) = 1;
      JAC(3,5)  = 1;
      JAC(3,12) = 1;
      JAC(3,4)  = 1;
      JAC(3,11) = 1;
      JAC(3,1)  = 1;
      JAC(3,8)  = 1;

      JAC(4,4)  = 1;
      JAC(4,11) = 1;
      JAC(4,6)  = 1;
      JAC(4,13) = 1;
      JAC(4,3)  = 1;
      JAC(4,10) = 1;
      JAC(4,1)  = 1;
      JAC(4,8)  = 1;

      JAC(5,5)  = 1;
      JAC(5,12) = 1;
      JAC(5,15) = 1;
      JAC(5,6)  = 1;
      JAC(5,13) = 1;
      JAC(5,1)  = 1;
      JAC(5,8)  = 1;

      JAC(6,6)  = 1;
      JAC(6,13) = 1;
      JAC(6,16) = 1;
      JAC(6,5)  = 1;
      JAC(6,13) = 1;
      JAC(6,1)  = 1;
      JAC(6,8)  = 1;

      JAC(7,7)  = 1;
      JAC(7,14) = 1;

    end
    %
    function H = ds_hessian( self, nseg, tL, tR, XL, XR, UC, P, L )
      H = self.FD_ds_hessian( nseg, tL, tR, XL, XR, UC, P, L );
    end
    %
    function H = ds_hessian_pattern( self )
      dim = 2*self.nx+self.nu+self.np;
      H   = sparse(ones(dim,dim));
    end
    %     _
    %  _ | |_  _ _ __  _ __
    % | || | || | '  \| '_ \
    %  \__/ \_,_|_|_|_| .__/
    %                 |_|
    %
    function JMP = jump( self, nsegL, t, XL, XR, P )
      sL  = XL(1); sR  = XR(1);
      nL  = XL(2); nR  = XR(2);
      vsL = XL(3); vsR = XR(3);
      vnL = XL(4); vnR = XR(4);
      asL = XL(5); asR = XR(5);
      anL = XL(6); anR = XR(6);

      dtheta = self.theta_L(nsegL+1)-self.theta_R(nsegL);

      S = sin(dtheta);
      C = cos(dtheta);

      JMP = [ ...
        sL-t, ...
        sR-t, ...
        nL, ...
        nR, ...
        vsR + S*vnL - C*vsL, ...
        vnR - C*vnL - S*vsL, ...
        asR + S*anL - C*asL, ...
        anR - C*anL - S*asL ...
      ];
    end

    %
    function JAC = jump_jacobian( self, nsegL, t, XL, XR, P )

      dtheta = self.theta_L(nsegL+1)-self.theta_R(nsegL);

      S = sin(dtheta);
      C = cos(dtheta);

      JAC = zeros(self.njmp,2*self.nx+self.np);
      JAC(1,1)  = 1;
      JAC(2,8)  = 1;
      JAC(3,2)  = 1;
      JAC(4,9)  = 1;

      JAC(5,10) = 1;
      JAC(5,4)  = S;
      JAC(5,3)  = -C;

      JAC(6,11) = 1;
      JAC(6,4)  = -C;
      JAC(6,3)  = -S;

      JAC(7,12) = 1;
      JAC(7,6)  = S;
      JAC(7,5)  = -C;

      JAC(8,13) = 1;
      JAC(8,6)  = -C;
      JAC(8,5)  = -S;
    end

    %
    function JAC = jump_jacobian_pattern( self )
      JAC = zeros(self.njmp,2*self.nx+self.np);
      JAC(1,1)  = 1;
      JAC(2,8)  = 1;
      JAC(3,2)  = 1;
      JAC(4,9)  = 1;

      JAC(5,10) = 1;
      JAC(5,4)  = 1;
      JAC(5,3)  = 1;

      JAC(6,11) = 1;
      JAC(6,4)  = 1;
      JAC(6,3)  = 1;

      JAC(7,12) = 1;
      JAC(7,6)  = 1;
      JAC(7,5)  = 1;

      JAC(8,13) = 1;
      JAC(8,6)  = 1;
      JAC(8,5)  = 1;
    end

    %
    function H = jump_hessian( self, nsegL, t, XL, XR, P, L )
      dim = 2*self.nx+self.np;
      H   = sparse(dim,dim);
    end

    %
    function H = jump_hessian_pattern( self )
      dim = 2*self.nx+self.np;
      H   = sparse(dim,dim);
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
      CJ = zeros(0,self.nx);
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

    %  ___  ___
    % | _ )/ __|
    % | _ \ (__
    % |___/\___|
    %
    function bc = bc( self, tL, tR, XL, XR, ~ )
      sL   = XL(1); sR   = XR(1);
      nL   = XL(2); nR   = XR(2);
      vsL  = XL(3); vsR  = XR(3);
      vnL  = XL(4); vnR  = XR(4);
      asL  = XL(5); asR  = XR(5);
      anL  = XL(6); anR  = XR(6);
      coVL = XL(7); coVR = XR(7);

      bc = [ sL  - self.s_i;  ...
             nL  - self.n_i;  ...
             vsL - self.vs_i; ...
             vnL - self.vn_i; ...
             asL - self.as_i; ...
             anL - self.an_i; ...
             sR  - self.s_f;  ...
             vsR - self.vs_f; ...
             vnR - self.vn_f; ...
             asR - self.as_f; ...
             anR - self.an_f ];
    end
    %
    function Jac = bc_jacobian( ~, tL, tR, XL, XR, ~ )
      Jac = zeros(12,14);
      Jac(1,1)   = 1;
      Jac(2,2)   = 1;
      Jac(3,3)   = 1;
      Jac(4,4)   = 1;
      Jac(5,5)   = 1;
      Jac(6,6)   = 1;
      Jac(7,8)   = 1;
      Jac(8,9)   = 1;
      Jac(9,10)  = 1;
      Jac(10,11) = 1;
      Jac(11,12) = 1;
      Jac(12,13) = 1;
    end
    %
    function Jac = bc_jacobian_pattern( ~ )
      Jac = zeros(12,14);
      Jac(1,1)   = 1;
      Jac(2,2)   = 1;
      Jac(3,3)   = 1;
      Jac(4,4)   = 1;
      Jac(5,5)   = 1;
      Jac(6,6)   = 1;
      Jac(7,8)   = 1;
      Jac(8,9)   = 1;
      Jac(9,10)  = 1;
      Jac(10,11) = 1;
      Jac(11,12) = 1;
      Jac(12,13) = 1;
    end
    %
    function Hess = bc_hessian( ~, tL, tR, XL, XR, ~, L )
      Hess = zeros(14,14);
    end
    %
    function Hess = bc_hessian_pattern( ~ )
      Hess = zeros(14,14);
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
