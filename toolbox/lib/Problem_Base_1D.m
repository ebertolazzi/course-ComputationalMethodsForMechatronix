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
classdef Problem_Base_1D < handle
  %
  % base class for store 1D problem in the calculus of variations
  %
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
    %    |     L( x, x', t ) dt
    %   / b
    %
    %   subject to  BC( x(a), x'(a), x(b), x'(b) ) = 0
    %
    %               / a
    %   subject to  |   IC( x, x', t ) dt = 0
    %               / b
    %
    % evaluate the function L(x,x',t)
    eval_Guess( self, t )

    % evaluate the function L(x,x',t)
    eval_L( self, x, x_dot, t )

    % evaluate the function DL(x,x',t) / Dx
    eval_L_D_1( self, x, x_dot, t )

    % evaluate the function DL(x,x',t) / Dx'
    eval_L_D_2( self, x, x_dot, t )

    % evaluate the function BC(x(a),x'(a),x(b),x'(b))
    eval_BC( self, xa, xa_dot, xb, xb_dot )

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx(a)
    eval_BC_D_1( self, xa, xa_dot, xb, xb_dot )

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx'(a)
    eval_BC_D_2( self, xa, xa_dot, xb, xb_dot )

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx(b)
    eval_BC_D_3( self, xa, xa_dot, xb, xb_dot )

    % evaluate the function DBC(x(a),x'(a),x(b),x'(b)) / Dx'(b)
    eval_BC_D_4( self, xa, xa_dot, xb, xb_dot )

    % evaluate the function IC(x,x',t)
    eval_IC( self, x, x_dot, t )

    % evaluate the function DIC(x,x',t)/Dx
    eval_IC_D_1( self, x, x_dot, t )

    % evaluate the function DIC(x,x',t)/Dx'
    eval_IC_D_2( self, x, x_dot, t )
  end

  methods
    function self = Problem_Base_1D( name, number_of_bc, number_of_ic )
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
