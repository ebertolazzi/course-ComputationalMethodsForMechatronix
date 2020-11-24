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
classdef Problem_Base_XY < handle
  %
  % base class for store 2D problem in the calculus of variations
  %
  properties (SetAccess = protected, Hidden = true)
    class_type;
    name;           % name of the problem
    number_of_bc;   % number of boundary conditions
    number_of_ic;   % number of integral conditions
  end

  methods (Abstract)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %     / a
    %    |     L( x, y, x', y', t ) dt
    %   / b
    %
    %   subject to  BC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) = 0
    %
    %               / a
    %   subject to  |   IC( x, y, x', y', t ) dt = 0
    %               / b
    %
    % evaluate the function L(x,x',t)
    eval_Guess( self, t )

    % evaluate the function L( x, y, x', y', t )
    eval_L( self, x, y, x_dot, y_dot, t )

    % evaluate the function DL(x, y, x', y', t) / Dx
    eval_L_D_1( self, x, y, x_dot, y_dot, t )

    % evaluate the function DL(x, y, x', y', t) / Dy
    eval_L_D_2( self, x, y, x_dot, y_dot, t )

    % evaluate the function DL(x, y, x', y', t) / Dx'
    eval_L_D_3( self, x, y, x_dot, y_dot, t )

    % evaluate the function DL(x, y, x', y', t) / Dy'
    eval_L_D_4( self, x, y, x_dot, y_dot, t )

    % evaluate the function BC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) )
    eval_BC( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dx(a)
    eval_BC_D_1( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dy(a)
    eval_BC_D_2( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dx'(a)
    eval_BC_D_3( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dy'(a)
    eval_BC_D_4( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dx(b)
    eval_BC_D_5( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dy(b)
    eval_BC_D_6( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dx'(b)
    eval_BC_D_7( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function DBC( x(a), y(a), x'(a), y'(a), x(b), y(b), x'(b), y'(b) ) / Dy'(b)
    eval_BC_D_8( self, xa, ya, xa_dot, ya_dot, xb, yb, xb_dot, yb_dot )

    % evaluate the function IC( x, y, x', y', t )
    eval_IC( self, x, y, x_dot, y_dot, t )

    % evaluate the function DIC( x, y, x', y', t )/Dx
    eval_IC_D_1( self, x, y, x_dot, y_dot, t )

    % evaluate the function DIC( x, y, x', y', t )/Dy
    eval_IC_D_2( self, x, y, x_dot, y_dot, t )

    % evaluate the function DIC( x, y, x', y', t )/Dx'
    eval_IC_D_3( self, x, y, x_dot, y_dot, t )

    % evaluate the function DIC( x, y, x', y', t )/Dy'
    eval_IC_D_4( self, x, y, x_dot, y_dot, t )
  end

  methods
    function self = Problem_Base_XY( name, number_of_bc, number_of_ic )
      %
      % initialize the class with the name of the problem
      % and the number of boundary condition and integral consitions
      self.class_type   = 'functional_1D';
      self.name         = name;
      self.number_of_bc = number_of_bc;
      self.number_of_ic = number_of_ic;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function delete( self )
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function name = get_name( self )
      name = self.name;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function t = get_type( self )
      t = self.class_type;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function nb = get_number_of_bc( self )
      nb = self.number_of_bc;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function ni = get_number_of_ic( self )
      ni = self.number_of_ic;
    end
  end
end
