classdef Beale
    % This class defines the following Beale function:
    % f(x1,x2) = (1.5-x1+x1*x2)^2+(2.25-x1+x1*x2^2)^2+(2.625-x1+x1*x2^3)^2
    
    % List of properties needed for the class.
    properties 
        % No properties needed for the Rosenbrock function.
    end
    
    methods
        
        % Evaluate objective function f at x.
        function f = func(self, x)
            if size(x,1) ~= 2
                error([inputname(2), ' must be a 2 dimensional vector!']);
            end
            x1 = x(1);
            x2 = x(2);
            f = (1.5-x1+x1*x2)^2+(2.25-x1+x1*x2^2)^2+(2.625-x1+x1*x2^3)^2;
        end
        
        % Evaluate the gradient of f at x.            
        function g = grad(self, x)
            if size(x,1) ~= 2
                error([inputname(2), ' must be a 2 dimensional vector!']);
            end
            x1 = x(1);
            x2 = x(2);
            g = [2*x2^3+2*x2^2+2*x2-6 ;
                 6*x1*x2^2+4*x1*x2+2*x1];
        end
        
        % Evaluate the Hessian of f at x.
        function h = hess(self, x)
            if size(x,1) ~= 2
                error([inputname(2), ' must be a 2 dimensional vector!']);
            end
            x1 = x(1);
            x2 = x(2);
            h = [       0     ,6*x2^2+4*x2+x;
                 6*x2^2+4*x2+x, 12*x1*x2+4*x1];
        end
        
        % Compute the product of the Hessian of f at x with vector x.
        function hv = hessvecprod(self, x, v)
            if size(x,1) ~= 2
                error([inputname(2), ' must be a 2 dimensional vector!']);
            end
            if size(v,1) ~= 2
                error([inputname(3), ' must be a 2 dimensional vector!']);
            end
            x1 = x(1);
            x2 = x(2);
            h = [       0     ,6*x2^2+4*x2+x;
                 6*x2^2+4*x2+x, 12*x1*x2+4*x1];
            hv = h*v;
        end
        
    end
end





