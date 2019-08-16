clear

Nm = 20; % Number of measurements
W = 90; % System size

%------ PARAMETERS ------%
a0 = 0.5;      a = 3.0;
K = 9.0;       n = 3.0;

%------ INPUT ------%
u = 0.15;
   
%------ INITIALIZATION ------%
X = 2.55*W; % choose IC
t = 0;

%------ HISTORY MEASUREMENTS ------%
Xt = zeros(Nm);
T = zeros(Nm);

dt = 50; % history size
dx = 50; % time displacement

%------ PLOTTING ------%
figure(1); clf; hold on; box on;
hand = plot(1:Nm,Xt,"o-","color", 'k');
xlabel("Time")
ylabel("Concentration")
set(gca,"xlim",[-dx 0]);
set(gca,"ylim",[0 4]);

%------ SIMULATION LOOP ------%
for i=1:inf
    for j=1:dt
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
    %------ STATE HISTORY UPDATE ------%
    Xt = [Xt(2:end),X/W];
    T = [T(2:end),t];
    %------ PLOTTING ------%
    set(hand,"xdata", T);
    set(hand,"ydata", Xt);
    set(gca,"xlim",[t-dx t]);
    %------ DISPLAY PLOT ------%
    %if mod(i,5) == 0
        drawnow
    %end
end