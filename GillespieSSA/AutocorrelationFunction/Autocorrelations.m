clear
Ns = 1e5; % Number of samples

%------ SYSTEM SIZE ------%
W = 120;

%------ PARAMETERS ------%
p = [0.5, 3.0, 9.0, 3.0];
a0 = p(1);      a = p(2);
K = p(3);       n = p(4);

%------ INPUT ------%
Nu = 12;
Nm = 300;
AC = zeros(Nm,Nu);

figure(6); clf; hold on; box off;
set(gca,'visible','off')
for ui=1:Nu
    tic
    u = 0.08 + 0.0125*ui;

    %------ JUMP TIME ------%
    tau_jump = 1e3;

    %------ ENSAMBLE OF SIMULATIONS ------%
    %------ INITIALIZE ------%
    t = 0; X = 0.55*W;
    %------ EQUILIBRATION ------%
    while t < 4.0*tau_jump % first jump
        %------ PROPENSITY ------%
        r = W*(a0 + a*(u+X/W)^n/(K + (u+X/W)^n));
        %------ TIME UPDATE ------%
        t = t + log(1/rand)/(r+X);
        %------ STATE UPDATE ------%
        if rand < r/(r+X)
            X = X+1;
        else
            X = X-1;
        end
    end
    %------ MEASUREMENTS ------%
    x_ave = 0;
    x2_ave = 0;
    Nm = 300;
    Xt = zeros(1,Nm);
    ts = 5;
    for i=1:Ns
        t=0;
        while t < ts
            %------ PROPENSITY ------%
            r = W*(a0 + a*(u+X/W)^n/(K + (u+X/W)^n));
            %------ TIME UPDATE ------%
            t = t + log(1/rand)/(r+X);
            %------ STATE UPDATE ------%
            if rand < r/(r+X)
                X = X+1;
            else
                X = X-1;
            end
        end
        x_ave = x_ave + X/W;
        x2_ave = x2_ave + X/W*X/W;
        Xt = [Xt(2:end),X/W];
        AC(:,ui) = AC(:,ui) + Xt.'*X/W;
    end

    toc
    AC(:,ui) = (AC(:,ui) - x_ave*x_ave/Ns)/(x2_ave - x_ave*x_ave/Ns);

    ax(ui) = axes('Position',[0+0.03*ui,0.0+0.07*ui,0.25,0.1]);
    plot(ax(ui),ts*(1:Nm),fliplr(AC(:,ui).'),'o','color','k','markersize',2.5);
    set(ax(ui),'Box','off');
    set(ax(ui),'Color','none');
    axis([0 ts*Nm -0.1 1])

    ax(Nu+ui) = axes('Position',[0.35+0.03*ui,0.0+0.07*ui,0.25,0.1]);
    plot(ax(Nu+ui),ts*(1:Nm),fliplr(AC(:,ui).'),'o','color','k','markersize',2.5);
    set(ax(Nu+ui),'Box','off');
    set(ax(Nu+ui),'Color','none');
    set(ax(Nu+ui),'xscale','log');
    axis([0 ts*Nm -0.1 1])

    drawnow
end