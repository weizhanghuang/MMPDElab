function [y,dt0,dt1] = MovFEM_radau45TS(fun_f,fun_jac,t,dt,y0,yt0,relTol, ...
                       absTol,fixed_step,direct_solver,ControlWeights,varargin)
%
% usage: [y,dt0,dt1] = MovFEM_radau45TS(fun_f,fun_jac,t,dt,y0,yt0,relTol, ...
%                      absTol,fixed_step,direct_solver,ControlWeights,varargin)
%
% this function implements the two-step fifth-order Radau IIA method with
% a two-step error estimator; for the latter, see
%
%   S. Gonzalez-Pinto, J. I. Montijano, and S. Perez-Rodrguez:
%   Two-step error estimators for implicit Runge-Kutta methods applied
%   to stiff systems,
%   ACM Trans. Math. Software, 30 (2004), 1-18.
%
% the system of ODEs is assumed to be in the implicit form, defined through
% functions
% 
%              f = fun_f(t,y,yt,varargin{:})
%   [dfdy,dfdyt] = fun_jac(t,y,yt,varargin{:})
%
% ControlWeights: (optional input) nonnegative vector of size length(y0)-by-1.
%       it defines the weights of the components of the solution used in
%       the error estimation for time step selection.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

preconditioning = true;
bicgstab_gmres = true; %% true => bicgstab                                                      

% define Butcher's Tableau for Radau5 (Radau IIA order 5)

   sqrt6 = sqrt(6.0);
   rkA = [(88-7*sqrt6)/360,(296-169*sqrt6)/1800,(-2+3*sqrt6)/225;
          (296+169*sqrt6)/1800,(88+7*sqrt6)/360,(-2-3*sqrt6)/225;
          (16-sqrt6)/36,(16+sqrt6)/36,1.0/9];
   rkC = sum(rkA,2);
   T = [9.1232394870892942792e-2,-0.14125529502095420843,-3.0029194105147424492e-2;
        0.24171793270710701896,0.20412935229379993199,0.38294211275726193779;
        0.96604818261509293619,1,0];
   gamma = 3.637834252744496;
   alpha = 2.681082873627752;
   beta = 3.050430199247410;
   Tinv = [gamma, 0, 0; 0, alpha, -beta; 0, beta, alpha]/T;
   u = 0.0000529585077373525889677785167637;
   %u = 1.0/2880;
   rkB = u*4/5*[19-14*sqrt(6),19+14*sqrt(6),52,-29-51*sqrt(6),-29+51*sqrt(6),-32];
      
% initialization

%   Mit = 200;
   Mit = 30;
   Mit_l = 1000;

   N = size(y0,1);
   Yt = zeros(N,3);
   Y = zeros(N,3);
   rhs = zeros(N,3);
   dt = 0.5*dt; %% for two steps
   
   TTinv = kron(Tinv,speye(N,N));
   TT = kron(T,speye(N,N));
   
   relTol_n = max(eps,relTol*0.01);
   relTol_n = min(sqrt(eps),relTol_n);
   absTol_n = max(eps,absTol*0.01);
   absTol_n = min(sqrt(eps),absTol_n);
      
   Tol_l = max(eps,min(relTol_n,absTol_n)*0.1);
   
   [JAC0,M0] = fun_jac(t,y0,yt0,varargin{:});
   I = sqrt(-1);
   
% begin integration

if (fixed_step)

      A = M0*(gamma/dt)+JAC0;
      B = M0*(alpha+I*beta)/dt+JAC0;
      if (preconditioning&&~direct_solver)
         [L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
         [LL,UU] = ilu(B,struct('type','ilutp','droptol',1e-6));
      end
      
    y = y0;
    for k=1:2 %% two steps
    
      if (k==1)
         Yt = repmat(yt0,1,3);
      else
         Yt = repmat(Yt(:,3),1,3);
         t = t + dt;
      end
      Y = repmat(y,1,3) + dt*Yt*rkA';
            
      for n=1:Mit  
         % compute rhs
         for i=1:3
            rhs(:,i)=fun_f(t+rkC(i)*dt,Y(:,i),Yt(:,i),varargin{:});            
         end
         if (norm(rhs)<=2*eps*sqrt(N*3)), break; end
         % solve the first linear system
         rhs(:) = TTinv*rhs(:)/dt;
         if (direct_solver)
            rhs(:,1) = A\rhs(:,1);
         else
            if (preconditioning)
               if (bicgstab_gmres)
                  [rhs(:,1),flag] = bicgstabl(A,rhs(:,1),Tol_l,Mit_l,L,U);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [rhs(:,1),flag] = gmres(A,rhs(:,1),20,Tol_l,Mit_l,L,U);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            else
               if (bicgstab_gmres)
                  [rhs(:,1),flag] = bicgstabl(A,rhs(:,1),Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [rhs(:,1),flag] = gmres(A,rhs(:,1),20,Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            end
            if (flag>0)
               fprintf('iteration for solving linear systems fails. flag = %d\n',flag); 
            end
         end
         % solve the second and third linear systems (complex system)
         R = rhs(:,2)+I*rhs(:,3);
         if (direct_solver)
            R = B\R;
         else
            if (preconditioning)
               if (bicgstab_gmres)
                  [R,flag] = bicgstabl(B,R,Tol_l,Mit_l,LL,UU);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [R,flag] = gmres(B,R,20,Tol_l,Mit_l,LL,UU);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            else
               if (bicgstab_gmres)
                  [R,flag] = bicgstabl(B,R,Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [R,flag] = gmres(B,R,20,Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            end
            if (flag>0)
               fprintf('iteration for solving linear systems fails. flag = %d\n',flag); 
            end
         end
         rhs(:,2) = real(R);
         rhs(:,3) = imag(R);
         % convert back to the original variable
         rhs(:) = TT*rhs(:);
         % update the solution
         Yt = Yt - rhs;
         Y = repmat(y,1,3) + dt*Yt*rkA';
         % test convergence
         if (norm(rhs)<=absTol_n*sqrt(N*3)), break; end
         if (norm(rhs)<=relTol_n*norm(Yt)), break; end
      end
      if (n == Mit) % Newton iteration fails
          error('integartion with fixed step fails (reach max iter number),%e\n',dt);
      end
      y = Y(:,end);
   end %% end of two steps
   t = t - dt;      
   dt0 = 2*dt;
   dt1 = 2*dt;
      
else % for variable stepping
   
   while (dt>2*eps)
   
      % Newton's iteration for first step
   
      A = M0*(gamma/dt)+JAC0;
      B = M0*(alpha+I*beta)/dt+JAC0;    
      if (preconditioning&~direct_solver)
         [L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
         [LL,UU] = ilu(B,struct('type','ilutp','droptol',1e-6));
      end

    flag_s = false;
    y = y0;
    for k=1:2 %% two steps
    
      if (k==1)
         Yt = repmat(yt0,1,3);
      else
         Yt = repmat(Yt(:,3),1,3);
         t = t + dt;
      end
      Y = repmat(y,1,3) + dt*Yt*rkA';
            
      for n=1:Mit
         % compute rhs
         for i=1:3
            rhs(:,i)=fun_f(t+rkC(i)*dt,Y(:,i),Yt(:,i),varargin{:});            
         end
         if (norm(rhs)<=2*eps*sqrt(N*3)), break; end
         % solve the first linear system
         rhs(:) = TTinv*rhs(:)/dt;
         if (direct_solver)
            rhs(:,1) = A\rhs(:,1);
         else
            if (preconditioning)
               if (bicgstab_gmres)
                  [rhs(:,1),flag] = bicgstabl(A,rhs(:,1),Tol_l,Mit_l,L,U);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [rhs(:,1),flag] = gmres(A,rhs(:,1),20,Tol_l,Mit_l,L,U);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            else
               if (bicgstab_gmres)
                  [rhs(:,1),flag] = bicgstabl(A,rhs(:,1),Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [rhs(:,1),flag] = gmres(A,rhs(:,1),20,Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            end
            if (flag>0)
               fprintf('iteration for solving linear systems fails. flag = %d\n',flag); 
            end
         end
         % solve the second and third linear systems (complex system)
         R = rhs(:,2)+I*rhs(:,3);
         if (direct_solver)
            R = B\R;
         else
            if (preconditioning)
               if (bicgstab_gmres)
                  [R,flag] = bicgstabl(B,R,Tol_l,Mit_l,LL,UU);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [R,flag] = gmres(B,R,20,Tol_l,Mit_l,LL,UU);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            else
               if (bicgstab_gmres)
                  [R,flag] = bicgstabl(B,R,Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('bicgstab fails. flag = %d\n',flag);
                  end
               else
                  [R,flag] = gmres(B,R,20,Tol_l,Mit_l);
                  if (flag > 0)
                     fprintf('gmres fails. flag = %d\n',flag);
                  end
               end
            end
            if (flag>0)
               fprintf('iteration for solving linear systems fails. flag = %d\n',flag); 
            end
         end
         rhs(:,2) = real(R);
         rhs(:,3) = imag(R);
         % convert back to the original variable
         rhs(:) = TT*rhs(:);
         % update the solution
         Yt = Yt - rhs;
         Y = repmat(y,1,3) + dt*Yt*rkA';
         % test convergence
         if (norm(rhs)<=absTol_n*sqrt(N*3)), break; end
         if (norm(rhs)<=relTol_n*norm(Yt)), break; end
      end
      if (n == Mit) % Newton iteration fails
          dt = dt*0.5;
          flag_s = true;
          break;
      end
      y = Y(:,end);
      if (k==1), Yt0 = Yt; end
   end %% end of two steps
   t = t - dt;
   if (flag_s), continue; end
      
      % estomate the error

      rhs(:,1) = dt*[Yt0,Yt]*rkB'; 
      rhs(:,1) = rhs(:,1).*ControlWeights./(absTol+relTol*max(abs(y),abs(y0)));
      delta = norm(rhs(:,1))/sqrt(2*N);
      if (delta<=1.0)
            delta1 = 0.9*delta^(-1.0/4);
            %delta1 = min(4.0, delta1);
            delta1 = min(1.5, delta1);
            dt1 = 2*delta1*dt;
            dt0 = 2*dt;
            break; % integration succeeds.
      else          
            delta1 = 0.9*delta^(-1.0/4);
            delta1 = min(0.1, delta1);
            delta1 = max(0.5, delta1);
            dt = delta1*dt;
      end
            
   end
   
   if (dt <= 2*eps)
      error('integartion fails with very small time stepsize, %e\n', dt);
   end
      
end

% end of MovFEM_radau45TS()
