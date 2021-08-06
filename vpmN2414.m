c = 1;
% Defining Chord Length

Vmag = 100.0; 
% Velocity Magnitude

AoA = 2.75; 
% Angle of Attack (alpha) of airfoil


airfoil = load('n2414.dat');

NumPan=length(airfoil)-1;
for i=1:NumPan
    point1(i,1) = airfoil(i,1);
    point1(i,2) = airfoil(i,2);
    point2(i,1) = airfoil(i+1,1);
    point2(i,2) = airfoil(i+1,2);
end

for i=1:NumPan
    dx = point2(i,1)-point1(i,1);
    dy = point2(i,2)-point1(i,2);
    dl(i) = sqrt( (point2(i,1)-point1(i,1))^2 + (point2(i,2)-point1(i,2))^2);
    th(i) = atan2(dy,dx);
    tnx(i) = cos(th(i));
    tny(i) = sin(th(i));
    vnx(i) = -tny(i);
    vny(i) = tnx(i);
end

for i=1:NumPan
    midp(i,1) = 0.5*(point1(i,1)+point2(i,1));
    midp(i,2) = 0.5*(point1(i,2)+point2(i,2));
end

for i=1:NumPan
    
    for j=1:NumPan
        
        xt = midp(i,1) - point1(j,1);
        yt = midp(i,2) - point1(j,2);
        x = xt*tnx(j) + yt*tny(j);
        y = -xt*tny(j) + yt*tnx(j);
        x1 = 0.0; y1 = 0.0;
        x2t = point2(j,1) - point1(j,1);
        y2t = point2(j,2) - point1(j,2);
        x2 = x2t*tnx(j) + y2t*tny(j);
        r1 = sqrt(x^2+y^2);
        r2 = sqrt((x-x2)*(x-x2)+y^2);
        th1 = atan2(y,x);
        th2 = atan2(y,x-x2);
        if(i==j) % self-induced velocity
            ax1 = 0.5*(x/x2-1.0);
            ay1 = 1.0/(2*pi);
            ax2 =-0.5*x/x2;
            ay2 =-1.0/(2*pi);
        else
            dth = th2-th1;
            rrt = r2/r1;
            rrtl = log(rrt);
            fcc = 1/(2*pi*x2);
            ax1 = fcc*( y*rrtl + (x-x2)*dth );
            ay1 = fcc*((x-x2)*rrtl - y*dth + x2);
            ax2 = -fcc*(y*rrtl + x*dth );
            ay2 = -fcc*(x*rrtl - y*dth + x2);
        end
        
        ux1 = ax1*tnx(j) - ay1*tny(j);
        uy1 = ax1*tny(j) + ay1*tnx(j);
        ux2 = ax2*tnx(j) - ay2*tny(j);
        uy2 = ax2*tny(j) + ay2*tnx(j);
        
        if(j==1)
            a(i,1)= ux1*vnx(i) + uy1*vny(i);
            holda = ux2*vnx(i) + uy2*vny(i);
        elseif(j==NumPan)
            a(i,NumPan) = ux1*vnx(i) + uy1*vny(i) + holda;
            a(i,NumPan+1) = ux2*vnx(i) + uy2*vny(i);
        else
            a(i,j)= ux1*vnx(i) + uy1*vny(i) + holda;
            holda = ux2*vnx(i) + uy2*vny(i);
        end
        
        if(j==1)
            b(i,1)= ux1*tnx(i) + uy1*tny(i);
            holdb = ux2*tnx(i) + uy2*tny(i);
        elseif(j==NumPan)
            b(i,NumPan) = ux1*tnx(i) + uy1*tny(i) + holdb;
            b(i,NumPan+1) = ux2*tnx(i) + uy2*tny(i);
        else
            b(i,j)= ux1*tnx(i) + uy1*tny(i) + holdb;
            holdb = ux2*tnx(i) + uy2*tny(i);
        end
    end
end

a(NumPan+1,1) = 1.0;%the Kutta condition
a(NumPan+1,NumPan+1) = 1.0;%the Kutta condition
alfa = AoA*pi/180.0;
coalf = cos(alfa); sialf = sin(alfa);
Ux = Vmag*coalf; Uy = Vmag*sialf;
for i=1:NumPan
    rhs(i) = -Ux*vnx(i)-Uy*vny(i);
end
rhs(NumPan+1)=0.0;

gamma = rhs/a';

circ = 0.0; % circulation
circgam = 0.0; % circulation in terms of gamma
for i=1:NumPan
    tnvel = Ux*tnx(i)+Uy*tny(i); % tangential velocity
    for j=1:NumPan+1
        tnvel = tnvel + b(i,j)*gamma(j);
    end
    circ = circ - tnvel*dl(i); % circulation
    circgam = circgam +0.5*(gamma(i)+gamma(i+1))*dl(i); % circulation
    cp(i) = 1.0-tnvel*tnvel/(Vmag*Vmag); % pressure coefficient
end
cp(NumPan+1)=cp(1);

hold on
patch(airfoil(:,1),-airfoil(:,2),'y')
plot(airfoil(:,1),cp,'--')
axis ij
xlabel('x')
ylabel('Cp')


cy = 0.0;
cx = 0.0;

for i=1:NumPan
    
    cy = cy-cp(i)*dl(i)*vny(i);
    cx = cx-cp(i)*dl(i)*vnx(i);
    
    
end


cy = cy/c;
cx = cx/c;

cd =  cx*coalf+cy*sialf; %drag coefficient in wind axis
cl = -cx*sialf+cy*coalf; %lift coefficient in wind axis

fprintf('CL= %f \n',cl)
fprintf('CD= %f \n',cd)

