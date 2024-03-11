clc;
clear;
clear all;

%Ahmet Arda Zobar
%Newton Rhapson Method
%--------------------------STEP 1---------------------------------

% constants
epsilon= 10^-10; %stopping criteria
a=250; 
b=25;
c=20;

%function and its handle
syms x1 x2
f= @(x1,x2) c + ( x1.^2 - a*cos(((2*pi)/b).*x1)) + (x2.^2 - a*cos(((2*pi/b)).*x2));


%%%display derivative of the initial functions%%
ff= @(x1,x2) 2*x1 + 2*x2 + 20*pi*sin((2*pi*x1)/25) + 20*pi*sin((2*pi*x2))/25;
% fff= @(x1,x2) (8*pi^2*cos((2*pi*x1)/25))/5 + (8*pi^2*cos((2*pi*x2)/25))/5 + 4 ;
diff_f1= diff(f,x1)+diff(f,x2); %  first derivative of the function
diff_f2= diff(ff,x1)+diff(ff,x2); %  second derivative of the function
disp('First derivative of the index function')
disp(diff_f1)
disp('Second derivative of the index function')
disp(diff_f2)
%%%%

%variables
xk= [40;45]; %x initial
x1= xk(1,1); % initially 40
x2= xk(2,1); % initially 45

%derivative matrices
g= [ 2*x1 + 20*pi*sin((2*pi*x1)/25) ; 2*x2 + 20*pi*sin((2*pi*x2)/25)]; %first derivative
h= [ (8*pi^2*cos((2*pi*x1)/25))/5 + 2 , 0 ; 0 , ((8*(pi^2)*cos((2*pi*x2)/25))/5)+2 ]; %second derivative

% %%%%% 3d meshgrid plot%%%%%
% this part is mostly same from my steepest descent homework
figure(1)
x=linspace(-50,50,100); 
y=linspace(-50,50,100);
[X,Y] = meshgrid(x,y);
surf(X,Y,f(X,Y));
surfc(X,Y,f(X,Y)); % contour command
shading flat; %in order to get a better visual on the marked point, I shade them with light colors
hold on
plot3(40,45,f(40,45),'ko') %black circle on data point
legend('3d function','initial point') %legend

%------------------------STEP 2-------------------
delta_xk= -inv(h)*g; % when dividing two matrices, you take the demonitor's inverse
%------------------------STEP 3-------------------
xk_later= xk + delta_xk; % x2=x1 + deltax
k=1; 
%------------------------STEP 4-------------------
while norm(delta_xk)>epsilon

    xk=xk_later;
    x1=xk(1,1);
    x2=xk(2,1);
    
    g= [ 2*x1 + 20*pi*sin((2*pi*x1)/25) ; 2*x2 + 20*pi*sin((2*pi*x2)/25)];
    h= [ (8*pi^2*cos((2*pi*x1)/25))/5 + 2 , 0 ; 0 , ((8*(pi^2)*cos((2*pi*x2)/25))/5)+2 ]; 
    
    delta_xk= -inv(h)*g;
    xk_later= xk + delta_xk;
    

    
    %%data part
    
    %we take squareroots of their sums of squares because theyre matrices, 
    %but we need skalars to get a graph of them
    gin(k)= sqrt (g(1)^2 + g(2)^2); % g input data 
    hinput(k)= sqrt (h(1,1)^2 + h(1,2)^2 + h(2,1)^2 + h(2,2)^2); 
      %doesn't work
    delta_xkin(k)= abs(sqrt(delta_xk(1)^2 + delta_xk(2)^2)); %delta_xk input data
    x1in(k)= x1; % x1 input data
    x2in(k)= x2; % x2 input data
    xin (k)= f(x1,x2); % f(x) input data
    
    
    if norm(delta_xk) < epsilon;
        break
    else 
            k=k+1;
   
    end
    
    
end

%------------------------------STEP 5----------------------

%%%subplots%%%%
figure (2)

subplot (3,1,1);
plot( 1 : length(xin),xin);
legend('fx[k]');
grid
hold on

subplot (3,1,2);
plot( 1 : length(x1in),x1in);
legend('x1[k]');
grid
hold on

subplot (3,1,3);
plot( 1 : length(x2in),x2in);
legend('x2[k]');
grid



figure (3)

subplot (3,1,1);
plot( 1 : length(gin),gin);
legend('{\nabla} f[k]');
grid
hold on

% subplot (3,1,2);
% plot( 1 : length(hinput),hinput);
% legend('h[k]');
% grid
% hold on

subplot (3,1,3);
plot( 1 : length(delta_xkin),delta_xkin);
legend('delta_xk[k]');
grid

figure(4)
x=linspace(-50,50,100); 
y=linspace(-50,50,100);
[X,Y] = meshgrid(x,y);
surf(X,Y,f(X,Y));
surfc(X,Y,f(X,Y)); % contour command
shading flat; %in order to get a better visual on the marked point, I shade them with light colors
hold on

plot3(x1in(2),x2in(2),f(x1in(2),x2in(2)),'ko')
plot3(x1in(4),x2in(4),f(x1in(4),x2in(4)),'ko')
plot3(x1in(3),x2in(3),f(x1in(3),x2in(3)),'ko') %black circle on data point
plot3(x1in(5),x2in(5),f(x1in(5),x2in(5)),'ko')
plot3(x1in(6),x2in(6),f(x1in(6),x2in(6)),'ko')
plot3(x1in(7),x2in(7),f(x1in(7),x2in(7)),'ko')
plot3(x1in(8),x2in(8),f(x1in(8),x2in(8)),'ko')
plot3(x1in(9),x2in(9),f(x1in(9),x2in(9)),'ko')
plot3(x1in(10),x2in(10),f(x1in(10),x2in(10)),'ko')
% plot3(x1in(11),x2in(11),f(x1in(11),x2in(11)),'ko')
legend('3d function','every iteration until optimal point') %legend

user_input=6; %arbitrary init

while user_input~=0 
    
user_input= input("Do you want to see the collected variable data? If so, choose one of them. \n 1)delta_x[k] \n 2)fx[k] \n 3)first partial derivative of fx[k] \n 4)second partial derivative of fx[k] \n 5)all of them \n 0)none of them, cancel it \n");

switch user_input
case 1 
    disp('delta_x[k]')
    disp(delta_xkin)
    case 2 
        disp('fx[k]')
        disp(xin)
    case 3
        disp('first partial derivative of fx[k]')
        disp(gin)
        case 4
        disp('second partial derivative of fx[k]')
        disp(hinput)
    case 5
        disp('delta_x[k]')
        disp(delta_xkin)
        disp('fx[k]')
         disp(xin)
         disp('partial derivative of fx[k]')
          disp(gin)
          disp('second partial derivative of fx[k]')
        disp(hinput)
    case 0
        if 5
            disp('ending program...')
        end
       
    otherwise 
        disp("Please enter a valid value. " )
        
end
end

