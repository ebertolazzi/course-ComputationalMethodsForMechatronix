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

classdef (Abstract) OCP_NLP < handle

  properties (SetAccess = protected, Hidden = true)
    nodes
    N    % total number of nodes

    nx   % number of states
    nu   % number of controls
    np   % number of parameters
    nbc  % number of boundary conditions

    npth % number of path constraints
    njmp % number of interface equations
    nseg % number of segments

    epsiReg
    kappa

    sol % stored soulution
  end

  methods (Abstract)
    % Lagrange target
    % nseg = number of the segment
    % tL   = left node
    % tR   = right node
    % XL   = left state
    % XR   = right state
    % UC   = control like state
    % PARS = free parameters
    L     = lagrange( self, nseg, tL, tR, XL, XR, UC, PARS )
    gradL = lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, PARS )
    hessL = lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, PARS )
    pat   = lagrange_hessian_pattern( self )

    % Mayer target
    M     = mayer( self, tL, tR, XL, XR, PARS )
    gradM = mayer_gradient( self, tL, tR, XL, XR, PARS )
    hessM = mayer_hessian( self, tL, tR, XL, XR, PARS )
    pat   = mayer_hessian_pattern( self )

    % Dynamical system part
    C   = ds( self, nseg, tL, tR, XL, XR, UC, PARS )
    CJ  = ds_jacobian( self, nseg, tL, tR, XL, XR, UC, PARS )
    pat = ds_jacobian_pattern( self )
    CH  = ds_hessian( self, nseg, tL, tR, XL, XR, UC, PARS, L )
    pat = ds_hessian_pattern( self )

    % Path constraints
    C   = pc( self, t, X, U, PARS )
    CJ  = pc_jacobian( self, t, X, U, PARS )
    pat = pc_jacobian_pattern( self )
    CH  = pc_hessian( self, t, X, U, PARS, L )
    pat = pc_hessian_pattern( self )

    % Jump condition
    % nsegL = number of the left segment, the right segment is nsegL+1
    jmp  = jump( self, nsegL, t, XL, XR, PARS )
    jmpJ = jump_jacobian( self, nsegL, t, XL, XR, PARS )
    pat  = jump_jacobian_pattern( self )
    jmpH = jump_hessian( self, nsegL, t, XL, XR, PARS, L )
    pat  = jump_hessian_pattern( self )

    % Boundary conditions
    bcf = bc( self, tL, tR, XL, XR, PARS )
    bcJ = bc_jacobian( self, tL, tR, XL, XR, PARS )
    pat = bc_jacobian_pattern( self )
    bcH = bc_hessian( self, tL, tR, XL, XR, PARS, L )
    pat = bc_hessian_pattern( self )

    % tvU
    tvU = TVU( self, dL, UCL, UCR )
    tvG = TVU_gradient( self, dL, UCL, UCR )
    tvH = TVU_hessian( self, dL, UCL, UCR )
    pat = TVU_hessian_pattern( self )
  end

  methods (Hidden = true)
    function check_pos_int( ~, msg, n )
      if ~ (isinteger(n) || floor(n) == n)
        error('%s: argument must be an integer, found %s', msg, class(n) );
      end
      if n < 0
        error('%s: argument must be a non negative integer, found %d', msg, n);
      end
    end

    function res = absreg( self, x )
      res = x*erf(self.kappa); % hypot( x, self.epsiReg ) - self.epsiReg;
    end

    function res = absreg_D( self, x )
      res = erf(self.kappa*x)+2*x.*exp(-self.kappa^2*x.^2)*self.kappa/sqrt(pi); % x ./ hypot( x, self.epsiReg );
    end

    function res = absreg_DD( self, x )
      res = (-4*self.kappa^3*x.^2+4*self.kappa).*exp(-self.kappa^2*x.^2)/sqrt(pi); % self.epsiReg^2 ./ hypot( x, self.epsiReg ).^3;
    end
  end

  methods
    %
    % nx  = number of states
    % nu  = number of controls
    % np  = number of free parameters
    % npc = number of path constraints
    % nbc = number of boundary conditions
    %
    function self = OCP_NLP( nx, nu, np, nbc, varargin )

      self.check_pos_int( 'nx (number of states)', nx );
      self.check_pos_int( 'nu (number of controls)', nu );
      self.check_pos_int( 'np (number of free parameters)', np );
      self.check_pos_int( 'nbc (number of boundary conditions)', nbc );

      self.nx  = nx;
      self.nu  = nu;
      self.np  = np;
      self.nbc = nbc;
      if nargin > 4
        self.check_pos_int( 'npth (number of path constraints)', varargin{1} );
        self.npth = varargin{1};
      else
        self.npth = 0;
      end
      if nargin > 5
        self.check_pos_int( 'njmp (number of interface conditions)', varargin{2} );
        self.njmp = varargin{2};
      else
        self.njmp = nx;
      end
      self.nseg    = 1;
      self.epsiReg = 0.0001;
      self.kappa   = 100;
    end

    %           _
    %  ___  ___| |_ _   _ _ __
    % / __|/ _ \ __| | | | '_ \
    % \__ \  __/ |_| |_| | |_) |
    % |___/\___|\__|\__,_| .__/
    %                    |_|
    function setup( self, nodes )
      self.nodes = nodes;
      self.N     = length(nodes);
      self.nseg  = 1;
      for k=2:self.N
        if nodes(k-1) == nodes(k)
          self.nseg = self.nseg+1;
        end
      end
    end

    %                   _
    %  _ __   __ _  ___| | __
    % | '_ \ / _` |/ __| |/ /
    % | |_) | (_| | (__|   <
    % | .__/ \__,_|\___|_|\_\
    % |_|
    %
    function XUP = pack( self, X, U, P )
      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      if size(X,1) ~= nx || size(X,2) ~= N
        error( ...
          'in function pack, X must be %i x %i, found %i x %i', ...
          nx, N, size(X,1), size(X,2) ...
        );
      end
      if size(U,1) ~= nu || size(U,2) ~= N-nseg
        error( ...
          'in function pack, U must be %i x %i, found %i x %i', ...
          nu, N, size(U,1), size(U,2) ...
        );
      end
      if min(size(P)) ~= 1 || max(size(P)) ~= np
        error( ...
          'in function pack, P must be %i x 1 or 1 x 1 found %i x %i', ...
          np, np, size(P,1), size(P,2) ...
        );
      end
      XUP = [ reshape( X, nx*N, 1 ); ...
              reshape( U, nu*(N-nseg), 1 ); ...
              reshape( P, np, 1 ) ];
    end
    %                               _
    %  _   _ _ __  _ __   __ _  ___| | __
    % | | | | '_ \| '_ \ / _` |/ __| |/ /
    % | |_| | | | | |_) | (_| | (__|   <
    %  \__,_|_| |_| .__/ \__,_|\___|_|\_\
    %             |_|
    %
    function [X,U,P] = unpack( self, XUP )
      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      totx = N*nx;
      totu = (N-nseg)*nu;
      llll = totx+totu+np;
      if size(UXP,1) ~= llll || self(X,2) ~= N
        error( ...
          'in function pack, X must be %i x %i, found %i x %i', ...
          nx, N, size(X,1), self(X,2) ...
        );
      end
      idx = 1:totx;
      X   = reshape( XUP(idx), nx, N );

      idx = totx+(1:totu);
      U   = reshape( XUP(idx), nu, N-nseg );

      if np > 0
        idx = totx+totu+(1:np);
        P   = reshape( XUP(idx), np, 1 );
      else
        P   = zeros(0,1);
      end
    end
    %      _        _
    %  ___| |_ __ _| |_ ___  ___
    % / __| __/ _` | __/ _ \/ __|
    % \__ \ || (_| | ||  __/\__ \
    % |___/\__\__,_|\__\___||___/
    %
    function X = states( self )
      N  = self.N;
      nx = self.nx;
      X  = reshape( self.sol(1:N*nx), nx, N );
    end

    %                  _             _
    %   ___ ___  _ __ | |_ _ __ ___ | |___
    %  / __/ _ \| '_ \| __| '__/ _ \| / __|
    % | (_| (_) | | | | |_| | | (_) | \__ \
    %  \___\___/|_| |_|\__|_|  \___/|_|___/
    %
    function U = controls( self )
      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      totx = N*nx;
      totu = (N-nseg)*nu;
      U    = reshape( self.sol(totx+1:totx+totu), nu, N - nseg );
    end

    function [UU,Unodes] = controls_for_plot( self )
      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      nseg = self.nseg;

      nn     = unique(self.nodes);
      dim    = 2*(N-nseg);
      Unodes = reshape( [ nn(1:end-1); nn(2:end)], 1, dim );
      U      = self.controls();
      UU     = zeros( nu, dim );
      UU(:,1:2:end-1) = U;
      UU(:,2:2:end)   = U;
    end
    %                                       _
    %  _ __   __ _ _ __ __ _ _ __ ___   ___| |_ ___ _ __ ___
    % | '_ \ / _` | '__/ _` | '_ ` _ \ / _ \ __/ _ \ '__/ __|
    % | |_) | (_| | | | (_| | | | | | |  __/ ||  __/ |  \__ \
    % | .__/ \__,_|_|  \__,_|_| |_| |_|\___|\__\___|_|  |___/
    % |_|
    %
    function P = parameters( self )
      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      nseg = self.nseg;
      totx = N*nx;
      totu = (N-nseg)*nu;
      P    = self.sol(totx+totu+1:end);
    end

    %  _                     _
    % | |_ __ _ _ _ __ _ ___| |_
    % |  _/ _` | '_/ _` / -_)  _|
    %  \__\__,_|_| \__, \___|\__|
    %              |___/
    %
    function res = NLP_target( self, Z )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      idx  = 1:nx;
      idx1 = (totx-nx)+(1:nx);
      idp  = (totx+totu)+(1:np);

      res = self.mayer( self.nodes(1), self.nodes(end), Z(idx), Z(idx1), Z(idp) );
      idu = totx+(1:nu);
      nsg = 1; % segment number
      cc  = zeros(1,ncu);
      kk  = 0;
      for k=1:N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk < nk1
          res    = res + self.lagrange( nsg, nk, nk1, Z(idx), Z(idx1), Z(idu), Z(idp) );
          idu    = idu + self.nu;
          kk     = kk+1;
          cc(kk) = nk1-nk;
        else
          nsg = nsg + 1; % next segment
        end
        idx = idx1;
      end

      % variation for controls
      idu = totx+(1:nu);
      for k=2:ncu
        idu1 = idu + nu;
        res  = res + self.TVU( cc(k)+cc(k-1), Z(idu), Z(idu1) );
        idu  = idu1;
      end

    end

    %
    function g = NLP_target_gradient( self, Z )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      idx  = 1:nx;
      idx1 = (totx-nx)+(1:nx);
      idp  = (totx+totu)+(1:np);

      g = zeros( 1, totx + totu + np );
      tmp = self.mayer_gradient( self.nodes(1), self.nodes(end), Z(idx), Z(idx1), Z(idp) );
      if size(tmp,1) ~= 1 || size(tmp,2) ~= 2*nx+np
        error( ...
          'call to ``mayer_gradient`` return a matrix of size %i x %i, expected 1 x %i\n', ...
          size(tmp,1), size(tmp,2), 2*nx+np ...
        );
      end
      g([idx,idx1,idp]) = tmp;
      idu = totx+(1:nu);
      nsg = 1; % segment number
      cc  = zeros(1,ncu);
      kk  = 0;

      for k=1:self.N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk < nk1
          id  = [ idx, idx1, idu, idp ];
          tmp = self.lagrange_gradient( nsg, nk, nk1, Z(idx), Z(idx1), Z(idu), Z(idp) );
          if size(tmp,1) ~= 1 || size(tmp,2) ~= 2*nx+nu+np
            error( ...
              'call to ``lagrange_gradient`` return a matrix of size %i x %i, expected 1 x %i\n', ...
              size(tmp,1), size(tmp,2), 2*nx+nu+np ...
            );
          end
          g(id)  = g(id) + tmp;
          idu    = idu + nu;
          kk     = kk+1;
          cc(kk) = nk1-nk;
        else
          nsg = nsg + 1; % next segment
        end
        idx = idx1;
      end

      % variation for controls
      idu = totx+(1:nu);
      for k=2:ncu
        idu1  = idu + nu;
        id    = [ idu, idu1 ];
        g(id) = g(id) + self.TVU_gradient( cc(k)+cc(k-1), Z(idu), Z(idu1) );
        idu   = idu1;
      end

    end

    %  _  _           _
    % | || |___ _____(_)__ _ _ _
    % | __ / -_|_-<_-< / _` | ' \
    % |_||_\___/__/__/_\__,_|_||_|
    %
    function H = NLP_hessian( self, Z, sigma, lambda )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      H    = sparse( totx+totu+np, totx+totu+np );

      idx  = 1:nx;
      idx1 = (totx-nx)+(1:nx);
      idp  = (totx+totu)+(1:np);
      idl  = (nx*(N-1))+(1:nbc);

      imap = [ idx, idx1, idp ];
      n1   = self.nodes(1);
      ne   = self.nodes(end);

      H(imap,imap) = sigma * self.mayer_hessian( n1, ne, Z(idx), Z(idx1), Z(idp) ) + ...
                     self.bc_hessian( n1, ne, Z(idx), Z(idx1), Z(idp), lambda(idl) );
      idu = totx+(1:nu);
      nsg = 1; % segment number
      cc  = zeros(1,ncu);
      kk  = 0;
      for k=1:self.N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk < nk1
          imap = [ idx, idx1, idu, idp ];
          H(imap,imap) = H(imap,imap) + ...
                         sigma * self.lagrange_hessian( nsg, nk, nk1, Z(idx), Z(idx1), Z(idu), Z(idp) ) + ...
                         self.ds_hessian( nsg, nk, nk1, Z(idx), Z(idx1), Z(idu), Z(idp), lambda(idx) );
          idu    = idu + nu;
          kk     = kk+1;
          cc(kk) = nk1-nk;
        else
          nsg = nsg + 1; % next segment
        end
        idx = idx1;
      end

      % variation for controls
      idu = totx+(1:nu);
      for k=2:ncu
        idu1     = idu + nu;
        id       = [ idu, idu1 ];
        H(id,id) = H(id,id) + sigma * self.TVU_hessian( cc(k)+cc(k-1), Z(idu), Z(idu1) );
        idu      = idu1;
      end
      H = tril(H);
    end

    function H = NLP_hessian_pattern( self )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      H    = sparse( totx+totu+self.np, totx+totu+self.np );

      idx  = 1:nx;
      idx1 = (totx-nx)+(1:nx);
      idp  = (totx+totu)+(1:np);

      imap         = [ idx, idx1, idp ];
      H(imap,imap) = self.mayer_hessian_pattern() + self.bc_hessian_pattern();
      idu          = totx+(1:nu);
      PAT          = self.lagrange_hessian_pattern() + self.ds_hessian_pattern();
      for k=1:self.N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk < nk1
          imap         = [ idx, idx1, idu, idp ];
          H(imap,imap) = H(imap,imap) + PAT;
          idu          = idu + nu;
        end
        idx = idx1;
      end

      % variation for controls
      if nu > 0
        idu = totx+(1:nu);
        PAT = self.TVU_hessian_pattern();
        for k=2:ncu
          idu1     = idu + nu;
          id       = [ idu, idu1 ];
          H(id,id) = H(id,id) + PAT;
          idu      = idu1;
        end
      end

      H = tril(H);

    end

    %   ___             _            _     _
    %  / __|___ _ _  __| |_ _ _ __ _(_)_ _| |_ ___
    % | (__/ _ \ ' \(_-<  _| '_/ _` | | ' \  _(_-<
    %  \___\___/_||_/__/\__|_| \__,_|_|_||_\__/__/
    %
    function C = NLP_constraints( self, Z )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      idx  = 1:nx;
      idu  = totx+(1:nu);
      idp  = (totx+totu)+(1:np);

      C = zeros( (N-nseg)*nx + (nseg-1)*njmp + nbc, 1 );

      nsg = 1; % segment number
      ic  = 0;
      for k=1:N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk1 > nk
          tmp = self.ds( nsg, nk, nk1, Z(idx), Z(idx1), Z(idu), Z(idp) );
          if size(tmp,1) ~= nx || size(tmp,2) ~= 1
            error( ...
              'call to ``constraints`` return a matrix of size %i x %i, expected %i x 1\n', ...
              size(tmp,1), size(tmp,2), nx ...
            );
          end
          C(ic+(1:nx)) = tmp;
          idu = idu + nu;
          ic  = ic + nx;
        else
          tmp = self.jump( nsg, nk, Z(idx), Z(idx1), Z(idp) );
          if size(tmp,1) ~= njmp || size(tmp,2) ~= 1
            error( ...
              'call to ``jump`` return a matrix of size %i x %i, expected %i x 1\n', ...
              size(tmp,1), size(tmp,2), nx ...
            );
          end
          C(ic+(1:njmp)) = tmp;
          ic  = ic + njmp;
          nsg = nsg+1;
        end
        idx = idx1;
      end

      idx = 1:nx;
      tmp = self.bc( self.nodes(1), self.nodes(end), Z(idx), Z(idx1), Z(idp) );
      if size(tmp,1) ~= self.nbc || size(tmp,2) ~= 1
        error( ...
          'call to ``bc`` return a matrix of size %i x %i, expected %i x 1\n', ...
          size(tmp,1), size(tmp,2), nbc ...
        );
      end
      C(ic+(1:nbc)) = tmp;
    end

    % -------
    function Jac = NLP_constraints_jacobian( self, Z )
      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      idx  = 1:nx;
      idu  = totx+(1:nu);
      idp  = (totx+totu)+(1:np);

      dimC = (N-nseg)*nx + (nseg-1)*njmp + nbc;
      dimZ = totx + totu + np;
      Jac  = sparse( dimC, dimZ );

      nsg = 1; % segment number
      ic  = 0;
      for k=1:N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk < nk1
          J = self.ds_jacobian( nsg, nk, nk1, Z(idx), Z(idx1), Z(idu), Z(idp) );
          if size(J,1) ~= nx || size(J,2) ~= 2*nx+nu+np
            error( ...
              'call to ``ds_jacobian`` return a matrix of size %i x %i, expected %i x %i\n', ...
              size(J,1), size(J,2), nx, 2*nx+nu+np ...
            );
          end
          Jac( ic+(1:nx), [ idx, idx1, idu, idp ] ) = J;
          idu = idu + nu;
          ic  = ic + nx;
        else
          J  = self.jump_jacobian( nsg, nk, Z(idx), Z(idx1), Z(idp) );
          if size(J,1) ~= nx || size(J,2) ~= 2*nx+np
            error( ...
              'call to ``jump_jacobian`` return a matrix of size %i x %i, expected %i x %i\n', ...
              size(J,1), size(J,2), nx, 2*nx+np ...
            );
          end
          Jac( ic+(1:njmp), [ idx, idx1, idp ] ) = J;
          ic  = ic + njmp;
          nsg = nsg+1;
        end
        idx = idx1;
      end
      idx  = 1:nx;
      J    = self.bc_jacobian( self.nodes(1), self.nodes(end), Z(idx), Z(idx1), Z(idp) );
      if size(J,1) ~= nbc || size(J,2) ~= 2*nx+np
        error( ...
          'call to ``bc_jacobian`` return a matrix of size %i x %i, expected %i x %i\n', ...
          size(J,1), size(J,2), nbc, 2*nx+np ...
        );
      end
      imap = ic + (1:nbc);
      jmap = [ idx, idx1, idp ];
      Jac(imap,jmap) = J;
    end

    % -------
    function Jac = NLP_constraints_jacobian_pattern( self )

      N    = self.N;
      nx   = self.nx;
      nu   = self.nu;
      np   = self.np;
      nseg = self.nseg;
      nbc  = self.nbc;
      njmp = self.njmp;

      totx = N*nx;
      ncu  = N-nseg;
      totu = ncu*nu;

      idx  = 1:nx;
      idu  = totx+(1:nu);
      idp  = (totx+totu)+(1:np);

      dimC = (N-nseg)*nx + (nseg-1)*njmp + nbc;
      dimZ = totx + totu + np;
      Jac  = sparse( dimC, dimZ );

      nsg = 1; % segment number
      ic  = 0;
      P1  = self.ds_jacobian_pattern();
      P2  = self.jump_jacobian_pattern();

      if size(P1,1) ~= nx || size(P1,2) ~= 2*nx+nu+np
        error( ...
          'call to ``ds_jacobian_pattern`` return a matrix of size %i x %i, expected %i x %i\n', ...
          size(P1,1), size(P1,2), nx, 2*nx+nu+np ...
        );
      end
      if size(P2,1) ~= nx || size(P2,2) ~= 2*nx+np
        error( ...
          'call to ``jump_jacobian_pattern`` return a matrix of size %i x %i, expected %i x %i\n', ...
          size(P2,1), size(P2,2), nx, 2*nx+nu+np ...
        );
      end

      for k=1:self.N-1
        nk   = self.nodes(k);
        nk1  = self.nodes(k+1);
        idx1 = idx + nx;
        if nk < nk1
          Jac( ic+(1:nx), [ idx, idx1, idu, idp ] ) = P1;
          idu = idu + nu;
          ic  = ic + nx;
        else
          Jac( ic+(1:njmp), [ idx, idx1, idp ] ) = P2;
          ic  = ic + njmp;
          nsg = nsg+1;
        end
        idx = idx1;
      end
      idx = 1:nx;
      imap = ic + (1:nbc);
      jmap = [ idx, idx1, idp ];
      tmp = self.bc_jacobian_pattern();
      if size(tmp,1) ~= nbc || size(tmp,2) ~= 2*nx+np
        error( ...
          'call to ``jump_jacobian_pattern`` return a matrix of size %i x %i, expected %i x %i\n', ...
          size(tmp,1), size(tmp,2), nbc, 2*nx+np ...
        );
      end
      Jac(imap,jmap) = tmp;
    end

    %  _   _ _   _ _
    % | | | | |_(_) |___
    % | |_| |  _| | (_-<
    %  \___/ \__|_|_/__/
    %
    function H = FD_ds_hessian( self, nseg, tL, tR, XL, XR, UC, PARS, L )

      nx = self.nx;
      nu = self.nu;
      np = self.np;

      id1 = 1:nx;
      id2 = nx+id1;
      id3 = 2*nx+(1:nu);
      id4 = 2*nx+nu+(1:np);

      GRAD = @(W) self.ds_jacobian( ...
        nseg, tL, tR, XL+W(id1), XR+W(id2), UC+W(id3), PARS+W(id4) ...
      ).' * L;

      % finite difference approximation of the hessian
      % Based on a code by Brendan C. Wood
      % Copyright (c) 2011, Brendan C. Wood <b.wood@unb.ca>
      NN = 2*nx+nu+np;
      H  = zeros(NN,NN);
      h  = max(1,abs([XL;XR;UC;PARS]))*eps^(1/3); % ricetta di cucina
      for i=1:NN
        % derivative at first point (left)
        x1    = zeros(NN,1);
        x1(i) = - h(i);
        df1   = GRAD(x1);

        % derivative at second point (right)
        x2    = zeros(NN,1);
        x2(i) = h(i);
        df2   = GRAD(x2);

        % differentiate between the two derivatives
        d2f = (df2-df1) ./ (2*h(i));

        % assign as column i of Hessian
        H(:,i) = d2f;
      end
      H = 0.5*(H+H.');
    end

    %     _
    %  _ | |_  _ _ __  _ __
    % | || | || | '  \| '_ \
    %  \__/ \_,_|_|_|_| .__/
    %                 |_|
    %
    function ODE = jump_standard( ~, nsegL, t, XL, XR, P )
      ODE = XR - XL;
    end

    %
    function JAC = jump_standard_jacobian( self, nsegL, t, XL, XR, P )
      JAC = sparse(self.nx, 2*self.nx+self.np);
      idx = 1:self.nx;
      JAC(idx,idx)         = -speye(self.nx,self.nx);
      JAC(idx,idx+self.nx) =  speye(self.nx,self.nx);
    end

    %
    function JAC = jump_standard_jacobian_pattern( self )
      JAC = sparse(self.nx, 2*self.nx+self.np);
      idx  = 1:self.nx;
      JAC(idx,idx)         = speye(self.nx,self.nx);
      JAC(idx,idx+self.nx) = speye(self.nx,self.nx);
    end

    %
    function H = jump_standard_hessian( self, nsegL, t, XL, XR, P, L )
      dim = 2*self.nx+self.np;
      H   = zeros(dim,dim);
    end

    %  _
    % | |   __ _ __ _ _ _ __ _ _ _  __ _ ___
    % | |__/ _` / _` | '_/ _` | ' \/ _` / -_)
    % |____\__,_\__, |_| \__,_|_||_\__, \___|
    %           |___/              |___/
    %
    function L = lagrange_zero( self, mneg, tL, tR, XL, XR, UC, P )
      L = 0;
    end
    %
    function gradL = lagrange_zero_gradient( self, mneg, tL, tR, XL, XR, UC, P )
      gradL = zeros(1,2*self.nx+self.nu+self.np);
    end
    %
    function hessL = lagrange_zero_hessian( self, mneg, tL, tR, XL, XR, UC, P )
      hessL = sparse(2*self.nx+self.nu+self.np,2*self.nx+self.nu+self.np);
    end
    %
    function hessL = lagrange_zero_hessian_pattern( self )
      hessL = sparse(2*self.nx+self.nu+self.np,2*self.nx+self.nu+self.np);
    end

    %  __  __
    % |  \/  |__ _ _  _ ___ _ _
    % | |\/| / _` | || / -_) '_|
    % |_|  |_\__,_|\_, \___|_|
    %              |__/
    %
    function M = mayer_zero( self, tL, tR, XL, XR, P )
      M = 0;
    end
    %
    function gradM = mayer_zero_gradient( self, tL, tR, XL, XR, P )
      gradM = zeros(1,2*self.nx+self.np);
    end
    %
    function hessM = mayer_zero_hessian( self, tL, tR, XL, XR, P )
      hessM = sparse(2*self.nx+self.np,2*self.nx+self.np);
    end
    %
    function hessM = mayer_zero_hessian_pattern( self )
      hessM = sparse(2*self.nx+self.np,2*self.nx+self.np);
    end

    %  _______     __                  _             _
    % |_   _\ \   / /   ___ ___  _ __ | |_ _ __ ___ | |___
    %   | |  \ \ / /   / __/ _ \| '_ \| __| '__/ _ \| / __|
    %   | |   \ V /   | (_| (_) | | | | |_| | | (_) | \__ \
    %   |_|    \_/     \___\___/|_| |_|\__|_|  \___/|_|___/
    % tvU
    function tvU = TVU_zero( self, epsi, dt, UCL, UCR )
      tvU = 0;
    end

    function tvG = TVU_zero_gradient( self, epsi, dt, UCL, UCR )
      tvG = zeros(1,2*self.nu);
    end

    function tvH = TVU_zero_hessian( self, epsi, dt, UCL, UCR )
      tvH = sparse(2*self.nu,2*self.nu);
    end

    function tvH = TVU_zero_hessian_pattern( self )
      tvH = sparse(2*self.nu,2*self.nu);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function tvU = TVU_standard( self, epsi, dt, UCL, UCR )
      tvU = sum( epsi .* (UCR-UCL).^2)/dt;
    end

    function tvG = TVU_standard_gradient( self, epsi, dt, UCL, UCR )
      tmp = (2/dt) * reshape( epsi .* (UCR-UCL), 1, self.nu);
      tvG = [ -tmp, tmp ];
    end

    function tvH = TVU_standard_hessian( self, epsi, dt, UCL, UCR )
      tmp = spdiags( 2 * epsi / dt, 0, self.nu,self.nu);
      tvH = sparse([ tmp, -tmp; -tmp, tmp ]);
    end

    function tvH = TVU_standard_hessian_pattern( self )
      tmp = speye(self.nu);
      tvH = [ tmp, tmp; tmp, tmp ];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    function tvU = TVU_reg( self, epsi, dt, UCL, UCR )
      tvU = sum( epsi .* self.absreg(UCR-UCL) );
    end
    %
    function tvG = TVU_reg_gradient( self, epsi, dt, UCL, UCR )
      tmp = reshape( epsi .* self.absreg_D(UCR-UCL), 1, self.nu);
      tvG = [ -tmp, tmp ];
    end
    %
    function tvH = TVU_reg_hessian( self, epsi, dt, UCL, UCR )
      tmp  = reshape( epsi .* self.absreg_DD(UCR-UCL), 1, self.nu);
      tmp2 = spdiags( reshape(tmp,self.nu,1), 0, self.nu,self.nu);
      tvH  = sparse([ tmp2, -tmp2; -tmp2, tmp2 ]);
    end
    %
    function tvH = TVU_reg_hessian_pattern( self )
      tmp = speye(self.nu);
      tvH = [ tmp, tmp; tmp, tmp ];
    end

    %        _    _           _     _
    %  _ __ (_)__| |_ __  ___(_)_ _| |_
    % | '  \| / _` | '_ \/ _ \ | ' \  _|
    % |_|_|_|_\__,_| .__/\___/_|_||_\__|
    %              |_|
    %
    function C = midpoint_ds( ~, nseg, tL, tR, XL, XR, UC, PARS, RHS )
      tM = (tR+tL)/2;
      XM = (XR+XL)./2;
      C  = (XR-XL)/(tR-tL) - feval( RHS, nseg, tM, XM, UC, PARS );
    end
    %
    function CJ = midpoint_ds_jacobian( self, nseg, tL, tR, XL, XR, UC, PARS, JAC )
      nx = self.nx;
      tM = (tR+tL)/2;
      XM = (XR+XL)./2;
      JJ = feval( JAC, nseg, tM, XM, UC, PARS );
      B1 = (-0.5)*JJ(1:nx,1:nx);
      B2 = JJ(1:nx,nx+1:end);
      bf = 1/(tR - tL);
      CJ = [ B1-bf*speye(nx), B1+bf*speye(nx), -B2 ];
    end
    %
    function CJ = midpoint_ds_jacobian_pattern( self, JAC_pattern )
      nx  = self.nx;
      pat = feval( JAC_pattern );
      B1  = pat(1:end,1:nx);
      B2  = pat(1:end,nx+1:end);
      CJ  = [ B1+speye(nx), B1+speye(nx), B2 ];
    end
    %
    function CH = midpoint_ds_hessian( self, nseg, tL, tR, XL, XR, UC, PARS, L, HESS )
      nx = self.nx;
      tM = (tR+tL)/2;
      XM = (XR+XL)./2;
      HH = feval( HESS, nseg, tM, XM, UC, PARS, L );
      D1 = (-0.25)*HH(1:nx,1:nx);
      R1 = (-0.5)*HH(1:nx,nx+1:end);
      D2 = -HH(nx+1:end,nx+1:end);
      CH = [ D1,   D1,   R1; ...
             D1,   D1,   R1; ...
             R1.', R1.', D2 ];
    end
    %
    function CH = midpoint_ds_hessian_pattern( self, HESS_pattern )
      nx  = self.nx;
      pat = feval( HESS_pattern );
      D1  = pat(1:nx,1:nx);
      R1  = pat(1:nx,nx+1:end);
      D2  = pat(nx+1:end,nx+1:end);
      CH  = [ D1,   D1,   R1; ...
              D1,   D1,   R1; ...
              R1.', R1.', D2 ];
    end

    %        _    _           _     _
    %  _ __ (_)__| |_ __  ___(_)_ _| |_
    % | '  \| / _` | '_ \/ _ \ | ' \  _|
    % |_|_|_|_\__,_| .__/\___/_|_||_\__|
    %              |_|
    %
    function C = midpoint_lagrange( ~, nseg, tL, tR, XL, XR, UC, PARS, LAG )
      tM = (tR+tL)/2;
      XM = (XR+XL)./2;
      C  = (tR-tL) * feval( LAG, nseg, tM, XM, UC, PARS );
    end

    function CJ = midpoint_lagrange_gradient( self, nseg, tL, tR, XL, XR, UC, PARS, GRAD )
      nx = self.nx;
      tM = (tR+tL)/2;
      XM = (XR+XL)./2;
      GG = feval( GRAD, nseg, tM, XM, UC, PARS );
      B1 = 0.5*GG(1,1:nx);
      B2 = GG(1,nx+1:end);
      CJ = (tR-tL) * [ B1, B1, B2 ];
    end

    function CH = midpoint_lagrange_hessian( self, nseg, tL, tR, XL, XR, UC, PARS, HESS )
      nx = self.nx;
      tM = (tR+tL)/2;
      XM = (XR+XL)./2;
      HH = feval( HESS, nseg, tM, XM, UC, PARS );
      D1 = 0.25*HH(1:nx,1:nx);
      R1 = 0.5*HH(1:nx,nx+1:end);
      D2 = HH(nx+1:end,nx+1:end);
      CH = (tR-tL) * [ D1,   D1,   R1; ...
                       D1,   D1,   R1; ...
                       R1.', R1.', D2 ];
    end

    function CH = midpoint_lagrange_hessian_pattern( self, HESS_pattern )
      nx  = self.nx;
      pat = feval( HESS_pattern );
      D1  = pat(1:nx,1:nx);
      R1  = pat(1:nx,nx+1:end);
      D2  = pat(nx+1:end,nx+1:end);
      CH  = [ D1,   D1,   R1; ...
              D1,   D1,   R1; ...
              R1.', R1.', D2 ];
    end
  end
end
