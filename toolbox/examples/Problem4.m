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
classdef Problem4 < Problem_Base_1D

  properties (SetAccess = protected)
    a;
    b;
    len;
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %     / a
    %    |     L( x, x', t ) dt
    %   / b
    %
    % evaluate the function L(x,x',t)
    function y = eval_Guess( self, x )
      y = (x-self.a)*(self.b-x);
    end

    % evaluate the function L(y,y',x)
    function L = eval_L( self, y, z, x )
      L = -y;
    end

    % evaluate the function DL(y,y',x) / Dy
    function L_D_1 = eval_L_D_1( self, y, z, t )
      L_D_1 = -1;
    end

    % evaluate the function DL(y,y',t) / Dy'
    function L_D_2 = eval_L_D_2( self, y, z, x )
      L_D_2 = 0;
    end

    % evaluate the function BC(y(a),y'(a),y(b),y'(b))
    function BC = eval_BC( self, ya, ya_dot, yb, yb_dot )
       BC = [ ya; yb];
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
    function IC = eval_IC( self, y, y_dot, t )
      IC = sqrt( 1 + y_dot^2 ) - self.len/(self.b-self.a);
    end

    % evaluate the function DIC(x,x',t)/Dx
    function IC_D_1 = eval_IC_D_1( self, y, y_dot, t )
      IC_D_1 = 0;
    end

    % evaluate the function DIC(x,x',t)/Dx'
    function IC_D_2 = eval_IC_D_2( self, y, y_dot, t )
      IC_D_2 = y_dot/sqrt( 1 + y_dot^2 );
    end
  end

  methods
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function self = Problem4( a, b, len )
      self@Problem_Base_1D( 'isoperimetric', 2, 1 );
      self.a   = a;
      self.b   = b;
      self.len = len;
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
