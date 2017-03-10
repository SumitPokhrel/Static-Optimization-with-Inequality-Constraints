% Sumit Pokhrel                                3/10/2013
% Static Optimization problem with inequality contraints 
% Solved using two methos: 
% ====================================================== 
% 1 a) Using FMINCON
clc
[x,minval,exit,out,lambda] = fmincon(@(x) 100*(x(2)-x(1)^2)^2+(1-x(1))^2,[.5;.5],-[1 0; 0 -1; -1 4; -0.5 -1],-[-1;-1;1;-1])

% 1 b) Using Penalty Function 
% Guess x and Tau
x = [.5 .5]';
sigma_2 = 1; tau = zeros(4,1);
W = sigma_2; %Fixed weight function, should be reasonably large 
M = 10;
N = 4;

% Outer loop for updating Tau 
for ii = 1:M; 
     % Bracket Operator 
     A = [1 0; 0 -1; -1 4; -.5 -1];B = [1;1;-1;1];
     C = A*x+B+tau; 
     for i = 1:4
         if C(i) > 0
             A(i,:)=0;
             B(i)=0;tau(i)=0;
         end 
     end 
     
 %Multiplier penalty Function minimization loop 
 Jx= ones (2,1);
 %  While norm(Jx) >=1e-6
 %performing descent steps to sufficiently decrease the penalty function
 %for a fixed tau
 for jj = 1:N
     %Gradient 
      jx = [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1));200*(x(2)-x(1)^2)]+2*W*A'*(A*x+B+tau);
      %Hessian 
      jxx = [-400*(x(2)-x(1)^2)+800*x(1)^2+2    -400*x(1);-400*x(1) 200]+2*W*A'*A;
      %Search Direction 
      S = -inv(jxx)*jx;
      %Optimal step length 
      alpha_0 = 0.5; % for simplicity . May have to change if it doesn't works 
      
      %Update x
      x = x + alpha_0*S
      L = 100*(x(2)-x(1)^2)^2+(1-x(1))^2
      Lamda2 = 2*sigma_2*tau
 end %Penalty function minimization ends 
 tau = A*x+B+tau; % Update tau after the end of the minimization iterations
                  % Step size control parameter may be required 
end 

%%Plots
f = @(x,y) 100*(y-x.^2).^2+(1-x).^2;
g1 = @(x,y) x+1 ;
g2 = @(x,y) 1-y;
g3 = @(x,y) 4*y-x-1;
g4 = @(x,y) 1-0.5*x-y;

% Plotting the constraints 
ezplot(g1,[-10,10,-10,10])
hold on
ezplot(g2,[-10,10,-10,10])
hold on
ezplot(g3,[-10,10,-10,10])
hold on
ezplot(g4,[-10,10,-10,10])

hold on
ezcontour(f,[-10,10,-10,10])
plot(.7813,.6094,'ro'); % Plotting the minimum 
legend('constraint1','constraint2','constraint3','constraint4','f contours','minimum');
hold off
saveas(gcf,'figure1.png')
      
      
    