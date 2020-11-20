%
% Matlab code for the Course:
%
%     Computational Methods for Mechatronix
%
% by
% Enrico Bertolazzi
% Dipartimento di Ingegneria Industriale
% Universita` degli Studi di Trento
% email: enrico.bertolazzi@unitn.it
%
classdef Problem3 < Problem_Base_1D

  properties (SetAccess = protected)
    a;
    b;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %     / a
    %    |     L( x, x', t ) dt
    %   / b
    %
    % evaluate the function L(x,x',t)
    function y = eval_Guess( self, x )
      y = -1;
    end

    % evaluate the function L(y,y',x)
    function L = eval_L( self, y, z, x )
      L = sqrt( max( (1+z*z)/(-y), 0 ) ); % avoid complex number
    end

    % evaluate the function DL(y,y',x) / Dy
    function L_D_1 = eval_L_D_1( self, y, z, t )
      tmp = 1+z^2;
      L_D_1 = -0.5*tmp/sqrt(-y*tmp);
    end

    % evaluate the function DL(y,y',t) / Dy'
    function L_D_2 = eval_L_D_2( self, y, z, x )
      tmp = 1+z^2;
      L_D_2 = z/sqrt(max(eps,-y*tmp)); % avoid complex number
    end

    % evaluate the function BC(y(a),y'(a),y(b),y'(b))
    function BC = eval_BC( self, ya, ya_dot, yb, yb_dot )
       BC = [ ya; yb+1];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy(a)
    function BC_D_1 = eval_BC_D_1( self, ya, ya_dot, yb, yb_dot )
      BC_D_1 = [ 1; 0];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy'(a)
    function BC_D_2 = eval_BC_D_2( self, ya, ya_dot, yb, yb_dot )
      BC_D_2 = [ 0; 0];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy(b)
    function BC_D_3 = eval_BC_D_3( self, ya, ya_dot, yb, yb_dot )
      BC_D_3 = [ 0; 1];
    end

    % evaluate the function DBC(y(a),y'(a),y(b),y'(b)) / Dy'(b)
    function BC_D_4 = eval_BC_D_4( self, ya, ya_dot, yb, yb_dot )
      BC_D_4 = [ 0; 0];
    end

    % evaluate the function IC(x,x',t)
    function IC = eval_IC( self, x, x_dot, t )
      IC = [];
    end

    % evaluate the function DIC(x,x',t)/Dx
    function IC_D_1 = eval_IC_D_1( self, x, x_dot, t )
      IC_D_1 = [];
    end

    % evaluate the function DIC(x,x',t)/Dx'
    function IC_D_2 = eval_IC_D_2( self, x, x_dot, t )
      IC_D_2 = [];
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Problem3( L )
      self@Problem_Base_1D( 'brachiostocrona', 2, 0 );
      self.a = 0;
      self.b = L;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function res = exact( self, t )
      % no exact solution for the moment
      res = 0;
    end
  end
end
