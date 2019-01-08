figure(1)
subplot(1,2,1)
%thick line
plot(Tsim,x_upd(:,1),Tsim,x_k(:,1),'--g','LineWidth',2)
set(gca, 'FontSize', 12)
%y & x lable
% xlabel('Time(min)','FontSize',20,'FontWeight','bold','Color','b')
% ylabel('C_A(Kmol/m^3)','FontSize',20,'FontWeight','bold','Color','b')
xlabel('Time(min)','FontSize',15,'FontWeight','bold')
ylabel('C_A(Kmol/m^3)','FontSize',15,'FontWeight','bold')
legend('Estimated State','Plant State')
%
% figure(2)
%thick line
subplot(1,2,2)
plot(Tsim,x_upd(:,2),Tsim,x_k(:,2),'--g','LineWidth',2)
set(gca, 'FontSize', 12)
%y & x lable
% xlabel('Time(min)','FontSize',20,'FontWeight','bold','Color','b')
% ylabel('C_A(Kmol/m^3)','FontSize',20,'FontWeight','bold','Color','b')
xlabel('Time(min)','FontSize',15,'FontWeight','bold')
ylabel('T(Kmol/m^3)','FontSize',15,'FontWeight','bold')
legend('Estimated State','Plant State')
%