function [dydt] = ODE_STAGE1(t, y, params)
% This ODE set simulates through aggregation of antibodies around antigen
% Assumptions
%   2 binding spots per antigen
%   All IgG subclasses bind with equal affinity
%   Once at steady state, the dissociation is limited enough to be
%   neglected before progression to the FcR engagement

%% Parameter assignments
f1 = params(1);%env-IgG1 kon
r1 = params(2);%koff
f2 = params(3);%env-IgG2
r2 = params(4);
f3 = params(5);%env-IgG3
r3 = params(6);
f4 = params(7);%env-IgG4
r4 = params(8);

%% Complex assignments
% Single bound antigen
e1 = y(1); 
e2 = y(2); 
e3 = y(3); 
e4 = y(4);

% Doubly bound antigen
e11 = y(5); 
e12 = y(6);
e13 = y(7); 
e14 = y(8); 
e22 = y(9); 
e23 = y(10);   
e24 = y(11);
e33 = y(12);
e34 = y(13);
e44 = y(14);

%% Parameter assignemnt (by reaction - might want to edit later)
k1f = 2*f1;
k1r = r1;
k2f = 2*f2;
k2r = r2;
k3f = 2*f3;
k3r = r3;
k4f = 2*f4;
k4r = r4;
k5f = f1;
k5r = 2*r1;
k6f = f2;
k6r = r2;
k7f = f1;
k7r = r1;
k8f = f3;
k8r = r3;
k9f = f1;
k9r = r1;
k10f = f4;
k10r = r4;
k11f = f1;
k11r = r1;
k12f = f2;
k12r = 2*r2;
k13f = f3;
k13r = r3;
k14f = f2;
k14r = r2;
k15f = f4;
k15r = r4;
k16f = f2;
k16r = r2;
k17f = f3;
k17r = 2*r3;
k18f = f4;
k18r = r4;
k19f = f3;
k19r = r3;
k20f = f4;
k20r = 2*r4;

% Conservation equations
g1 = g1tot-(e1+2*e11+e12+e13+e14);
g2 = g2tot-(e2+2*e22+e12+e23+e24);
g3 = g3tot-(e3+2*e33+e23+e13+e34);
g4 = g4tot-(e4+2*e44+e24+e34+e14);
e = etot-(e1+e2+e3+e4+e11+e12+e13+e14+e22+e23+e24+e33+e34+e44);

% Reaction rates
react1  = k1f*g1*e-k1r*e1;
react2  = k2f*g2*e-k2r*e2;
react3  = k3f*g3*e-k3r*e3;
react4  = k4f*g4*e-k4r*e4; 
react5  = k5f*g1*e1-k5r*e11;
react6  = k6f*g2*e1-k6r*e12;
react7 = k7f*g1*e2-k7r*e12;
react8 = k8f*g3*e1-k8r*e13;
react9 = k9f*g1*e3-k9r*e13;
react10 = k10f*g4*e1-k10r*e14;
react11 = k11f*g1*e4-k11r*e14;
react12 = k12f*g2*e2-k12r*e22; 
react13 = k13f*g3*e2-k13r*e23;
react14 = k14f*g2*e3-k14r*e23;
react15 = k15f*g4*e2-k15r*e24;
react16 = k16f*g2*e4-k16r*e24;
react17 = k17f*g3*e3-k17r*e33;
react18 = k18f*g4*e3-k18r*e34;
react19 = k19f*g3*e4-k19r*e34;
react20 = k20f*g4*e4-k20r*e44;

de1 = react1-react5-react6-react8-react10;
de2 = react2-react7-react12-react13-react15;
de3 = react3-react9-react14-react17-react18;
de4 = react4-react11-react16-react19-react20;
de11 = react5;
de12 = react6+react7;
de13 = react8+react9;
de14 = react10+react11;
de22 = react12;
de23 = react13+react14;
de24 = react15+react16;
de33 = react17;
de34 = react18+react19;
de44 = react20;

% [nM/s] product
dydt=[de1;de2;de3;de4;de11;de12;de13;de14;de22;de23;de24;de33;de34;de44];

% Manipulation of global variables for storing fluxes, intermediate variables
global tStep tArray varArray;
if t > tArray(tStep)            % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;  
end
tArray(tStep) = t;              
varArray(tStep,:) = [g1,g2,g3,g4,e,f];

end