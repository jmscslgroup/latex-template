% (C) RAHUL BHADANI
% Problem 5

% Kennedy Receiver - BPSK

N = 0.0000001:0.0001:10;
P = zeros(size(N));
b_root = zeros(size(N));

for j = 1:length(N)
    a = sqrt(N(j));
    
    %Newton Raphson Method, %Find the roots of dPE/dbeta = 0
    bn = 0.5; %Initial Guess of the roopts of dPE/dbeta
    %We are doing 10000 iterations
    for i  = 1:10000
        b_next = bn - (f(bn,a)./fdash(bn,a));
        bn = b_next;
    end
    b_root(j) = b_next;
    P(j) = PE(b_next, a);

end

fig= figure;
fig.Position = [652   490   593   479];
PEDD = 0.5*exp(-4.*N);

P_Dolinar = 0.5.*(1- sqrt(1 - exp(-4.*N)));
subplot(2,1,1);
semilogx(N, PEDD,'r','LineWidth',1);
hold on;
semilogx(N, P,'b--','LineWidth',1);
semilogx(N, P_Dolinar,'k-','LineWidth',1);
xlabel('N','Interpreter','latex');
ylabel('$P_e$ Average probability of the error','Interpreter','latex');
set(gca,'FontSize',16);
set(gca,'FontName','Times');
title('$P(error)$ for the Kennedy receiver','Interpreter','latex');
grid on;
grid minor;
legend({'$\textbf{Direct Detection-Equal Prior}$',...
    '$\textbf{Kennedy Receiver-Equal Prior}$',...
    '$\textbf{Dolinar Receiver}$'...
    },'Interpreter','latex');

subplot(2,1,2);
plot(N, b_root, 'g--','LineWidth',2);
xlabel('N','Interpreter','latex');
ylabel('$\beta$','Interpreter','latex');
set(gca,'FontSize',16);
set(gca,'FontName','Times');
title('$\textbf{Optimal Beta as a function of N}$','Interpreter','latex');
grid on;
grid minor;

% savefig(fig, 'figures/rahulbhadani_OPTI_595B_Q5.fig');
% saveas(gcf,'figures/rahulbhadani_OPTI_595B_Q5.pdf');

function PRET = PE(b, a)
    PRET = 0.5*(1 - exp(-b.*b)) + 0.5*exp(-( 2.*a + b).^2);
end

function fRET = f(b, a)
    fRET = b*exp(-b*b) - (b+2.*a)*exp(-(b+2.*a)^2);
end

function fdashRET = fdash(b, a)
    fdashRET = 2*((b+2.*a)^2)*exp(-(b+2.*a)^2) - exp(-(b+2.*a)^2) - 2*b*b*exp(-b*b) + exp(-b*b);
end
