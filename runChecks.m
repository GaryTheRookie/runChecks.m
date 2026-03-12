function checksTbl = runChecks(params, sweepTbl)
%RUNCHECKS (Project C) Learner Template
%
% Signature must match solution:
%   checksTbl = runChecks(params, sweepTbl)
%
% REQUIRED checksTbl columns (names + ORDER must match exactly):
%   check
%   error
%   pass
%
% Notes:
% - Exactly 3 checks are required (as defined in the handout).
% - Do NOT write any files here.
% - sweepTbl is expected to contain (at minimum) t and ymax.

checks = strings(3,1);
err    = zeros(3,1);
pass   = false(3,1);

%% Check 1 (handout): Zero load => zero stress and zero deflection
checks(1) = "Zero load: w=0 gives max|sigma|=0 and max|y|=0";
% TODO: form a temporary parameter struct with w set to zero

p = params;
p.w = 0;  % Set the load parameter to zero

% TODO: compute max stress and max deflection at the baseline thickness
% [smax, ymax] = beamMax(p, t_baseline);

[smax, ymax] = beamMax(p, t_baseline);

% TODO: set numeric residual (a single scalar) and pass/fail

err(1)  = max(abs([smax ymax]));
pass(1) = err(1) < 1exp(-12)

%% Check 2 (handout): Monotone stiffness over sweep
checks(2) = "Monotone stiffness: ymax decreases as thickness increases";
% TODO: compute numeric residual that quantifies violation of monotone decrease
%       (assume sweepTbl.t is increasing)

y = sweepTbl.ymax;
viol = max(diff(y));
err(2)  = max(0, viol);
pass(2) = err <= 0;

%% Check 3 (handout): Scaling (y ~ 1/I)
checks(3) = "Scaling: ymax ratio matches inverse I ratio (y~1/I)";
% TODO: choose two thickness values t1, t2 from the sweep (e.g., endpoints)

t1 = sweepTbl.t(1);
t2 = sweepTbl.t(end);

% TODO: compute ymax(t1), ymax(t2) using a helper (or by re-running a baseline calc)

[~, y1] = beamMax(params, t1);
[~, y2] = beamMax(params, t2);

% TODO: compute I(t) using the handout relationship and form the comparison ratio

I1 = params.b*t1^3/12;
I2 = params.b*t2^3/12;

ratio_y = y2 / y1;
ratio_I = I1 / I2;

% TODO: set a relative-error residual and pass/fail
err(3)  = abs(ratio_y - ratio_I)/abs(ratio_I);
pass(3) = err(3) < 1exp(-3);

%% Assemble table (EXACT names + ORDER)
checksTbl = table(checks, err, pass, 'VariableNames', {'check','error','pass'});

end

% --- Optional local helper stub ---
function [sigmax, ymax] = beamMax(params, t)
%BEAMMAX (Project C) Learner helper 
% TODO: Using the Project C handout formulas, compute:
%   sigmax : maximum bending stress magnitude
%   ymax   : maximum deflection magnitude
%
% Inputs:
%   params : struct containing at minimum L, b, E, w (and any needed constants)

L = params.L;
b = params.b;
E = params.E;
w = params.w;

I = b*t^3/12;
Mmax = w*L^2/8;
c = t/2;

%   t      : thickness (m)

sigmax = abs(Mmax*c/I);
ymax   = abs((5*w*L^4)/(384*E*I));

% TODO: implement per handout
end