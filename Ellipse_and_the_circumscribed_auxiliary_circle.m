%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg

% Ellipse and the circumscribed auxiliary circle
% Get normal and tangential velocity around ellipse 
% Get normal and tangential acceleration around ellipse 
close all; clear all; clc;
%% inputs
e=0.8; % eccentricity
a=5;   % semi-major axis
theta=linspace(0,2*pi,200); % theta vctor
mue=398600; % mu of the earth
hz=1000; % simulation frequancy
%% calculations
E=2*atan(sqrt((1-e)/(1+e))*tan(theta/2)); % Eccentric anomaly
for mm=1:length(E)
    if E(mm) < 0
        E(mm)=E(mm)+2*pi;
    end
end
r=a*(1-e^2)./(1+e*cos(theta));  % ellipce equation
xe=r.*cos(theta);  % ellipce x-axis equation
ye=r.*sin(theta);  % ellipce y-axis equation
xc=a.*cos(E)-a*e;  % circle x-axis equation
yc=a.*sin(E);  % circle y-axis equation
Vc=sqrt(mue/a)*ones(1,length(theta)); % circle velocity
Ve=sqrt(2*(mue*r.^-1+mue/2/a)); % ellipce velocity
%% plotting
figure(1);
nn=1;
% while figure(1) && nn==1
set(gcf,'color','w');
for n=1:length(E)
    if figure(1)==false && nn==1
        nn=0;
    elseif figure(1)==true && nn==1
    subplot(2,2,1);
    % ellipce
    plot(xe,ye,'color','black','linewidth',2);
    % circle
    hold;
    plot(xc,yc,'color','red','linewidth',2);
    % ra
    plot([0,a*(1-e)],[0,0],'color',[0.7,0.2,0.5],'linewidth',2);
    % rp
    plot([0,-a*(1+e)],[0,0],'color',[0.5,0.25,0.7],'linewidth',2);
    hold on;
    grid on;
    vectarrow([0,0],[xe(n),ye(n)],2,[0.1,0.5,0.4]);
    vectarrow([-a*e,0],[xc(n),yc(n)],2,[0.1,0.5,0.9]);
    vectarrow([xe(n),ye(n)],[xc(n),yc(n)],2,[0.9,0.5,0.4]);
    hold off;
    ylabel('Y-axis','Fontsize',18);
    xlabel('X-axis','Fontsize',18);
    axis equal;
    legend('Ellipse','Circle','ra','rp','r_e','r_c','\Delta r');
    ax1=subplot(2,2,2);
    hold all;
    gscatter(Ve(n),E(n));
    grid on;
    axis([ax1],[min(Ve) max(Ve) 0 2*pi]);
    ylabel('E','Fontsize',18);
    xlabel('V_e_l_l_i_p_s_e','Fontsize',18);
    ax2=subplot(2,2,3);
    hold all;
    gscatter(norm([xe(n),ye(n)]-[xc(n),yc(n)]),theta(n));
    grid on;
    axis([ax2],[0 a*e 0 2*pi]);
    ylabel('\theta','Fontsize',18);
    xlabel('\Delta r','Fontsize',18);
    ax3=subplot(2,2,4);
    hold all;
    gscatter(theta(n),E(n));
    grid on;
    axis([ax3],[0 2*pi 0 2*pi])
    ylabel('E','Fontsize',18);
    xlabel('\theta','Fontsize',18);
    pause(1/hz);
    end
end
% clf;
% end




