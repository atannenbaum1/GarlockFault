% %%%%% Fitting Garlock Fault Magnetic Anomaly %%%%%

%clear, previous work
clear, clc ,close all

%load in data
load('lat11.mat')
load('long11.mat')
load('fld11.mat')
load('elv11.mat')
%break of survey line 
lat = a(1:111); 
lon = b(1:111);
B = c(1:111);

%setting up traverse line and lat long conversions
theta = 19;                                                             %Set rotation angle for data projection along normal to fault (0=E/W, 90=N/S).
fault_norm=[sin(theta/32),cos(theta/32)];                                   %Unit vector of normal to fault strike in x/y.
lat_ref=min(lat); lon_ref=min(lon);                                        %Set reference latitude and longitude(to use for flattening) as min lat and min lon.
xd=GCD(lat_ref,lon,lat_ref,lon_ref);                                       %Flatten geometry to where x is displacement due east.
yd=GCD(lat,lon_ref,lat_ref,lon_ref);                                       %Flatten geometry to where y is displacement due north.
dist=xd.*fault_norm(1)+yd.*fault_norm(2);                                  %Calculate each point's distance from reference projected along theta.

%plot raw data over distance
f1 = figure(1);
plot(dist,B,'b')
hold on, grid on
label('Horizontal Distance (m)','Field (nT)','Garlock Fault Line 1: Raw Data',14)

%%%%%% Line of Poles Model %%%%%%
%a=[X offset from reference(m), depth(m), pole intensity per unit length,magnetic shift (Y axis)(nT),fault rotation angle]  
a0=[130, 37, 1.8, 47400, 150, 4, 0.15, 47760, 3.8]; %guess vectors
f=MAG_MODEL2(a0,dist);  %plug into forward model                                                    
[a, R,J, SIG]=nlinfit(dist,B,@MAG_MODEL2,a0);  %obtain fitted parameters
devpole =sqrt(diag(SIG)); %standard deviations of fitted parameters
f_fit=MAG_MODEL2(a,dist);  

%plot raw data, forward model, and fitted model
f2 = figure(2); 
plot(dist,B,'b')
hold on, grid on
plot(dist,f,'g--')
plot(dist,f_fit,'--r');                 
legend('Raw Data','Forward Model','Inverted Model','location', 'northwest')
label('Horizontal Distance (m)','Field (nT)','Garlock Fault Line 1: Line of Poles Model',14)
ssq = (sum((B - f_fit).^2))./length(B);

%%%%%% dipping trapezoidal model %%%%%%%%
f3 = figure(3);
%beta0 = [location, sus, depth to top, yshift, location2, sus2, depth to top2, width2]
beta0 = [120, 2.8e-4, 5, 30, 47732, 145, 2.9e-4 , 0.872, 500]; %guess vectors
[F] = dipdikemodel2(beta0,dist); %forward model

%plot raw data and forward model
plot(dist, B, 'b')
hold on, grid on
plot(dist, F, '--g')

%nlinfit 
[F1guess,R, J, SIG] = nlinfit(dist,B,@dipdikemodel2,beta0);
devtrap =sqrt(diag(SIG)); %std dev
F1 = dipdikemodel2(F1guess,dist);

%plot overlay of inverted model
plot(dist, F1, '--r')
legend('Raw Data','Forward model','Inverted model','location', 'northwest')
label('Horizontal Distance (m)','Field (nT)','Garlock Fault Line 1: Dipping Trapezoid Model', 14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Second Survey line %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup locations for second survey line
load('lat11.mat')
load('long11.mat')
load('fld11.mat')
load('elv11.mat')
lat = a(169:220); 
lon = b(169:220);
B = c(169:220);

%setting up traverse line and lat long conversions
theta = 23;                                                             %Set rotation angle for data projection along normal to fault (0=E/W, 90=N/S).
fault_norm=[sin(theta/29),cos(theta/29)];                                   %Unit vector of normal to fault strike in x/y.
lat_ref=min(lat); lon_ref=min(lon);                                        %Set reference latitude and longitude(to use for flattening) as min lat and min lon.
xd=GCD(lat_ref,lon,lat_ref,lon_ref);                                       %Flatten geometry to where x is displacement due east.
yd=GCD(lat,lon_ref,lat_ref,lon_ref);                                       %Flatten geometry to where y is displacement due north.
dist=xd.*fault_norm(1)+yd.*fault_norm(2);                                  %Calculate each point's distance from reference projected along theta.

%plot raw data over distance
f4 = figure(4);
plot(dist,B,'b')
hold on, grid on
label('Horizontal Distance (m)','Field (nT)','Garlock Fault Line 2: Raw Data',14)

%%%%%% Line of Poles Model %%%%%%
%a=[X offset from reference(m), depth(m), pole intensity per unit length,magnetic shift (Y axis)(nT),fault rotation angle]  
a0=[40, 17, 0.4, 47797]; %guess vectors
f=MAG_MODEL1(a0,dist);  %plug into forward model                                                    
[a, R,J, SIG]=nlinfit(dist,B,@MAG_MODEL1,a0);  %obtain fitted parameters
devpole =sqrt(diag(SIG)); %standard deviations of fitted parameters
f_fit=MAG_MODEL1(a,dist);  

%plot raw data, forward model, and fitted model
f5 = figure(5); 
plot(dist,B,'b')
hold on, grid on
plot(dist,f,'g--')
plot(dist,f_fit,'--r');                 
legend('Raw Data','Forward Model','Inverted Model','location', 'northwest')
label('Horizontal Distance (m)','Field (nT)','Garlock Fault Line 2: Line of Poles Model',14)
ssq = (sum((B - f_fit).^2))./length(B);

%%%%%% dipping trapezoidal model %%%%%%%%
f6 = figure(6);
%beta0 = [location, sus, depth to top, yshift, location2, sus2, depth to top2, width2]
beta0 = [39, 3.4e-4, 1.9, 106, 47787, 140]; %guess vectors
[F] = dipdikemodel3(beta0,dist); %forward model

%plot raw data and forward model
plot(dist, B, 'b')
hold on, grid on
plot(dist, F, '--g')

%nlinfit 
[F1guess,R, J, SIG] = nlinfit(dist,B,@dipdikemodel3,beta0);
devtrap =sqrt(diag(SIG)); %std dev
F1 = dipdikemodel3(F1guess,dist);

%plot overlay of inverted model
plot(dist, F1, '--r')
legend('Raw Data','Forward model','Inverted model','location', 'northwest')
label('Horizontal Distance (m)','Field (nT)','Garlock Fault Line 2: Dipping Trapezoid Model', 14)

%save for image 2
f =[f1, f2, f3, f4, f5, f6];
for i = 1:length(f)
set(f(i), 'Position', [0 0 1400 900])
f(i).PaperPositionMode = 'auto';
print(f(i), '-r250','-dpng',['Tannenbaum-136Cfinalproject', num2str(i)])
end
