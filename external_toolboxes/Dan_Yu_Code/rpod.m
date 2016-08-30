function [Ar, Br, Cr, VF, UF] = rpod(A, B, C, NUM_NON, NUM_PRIMAL, NUM_STEP)

[NUM_SYS, NUM_IN] = size(B);
[NUM_OUT, NUM_SYS] = size(C);



%NUM_PRIMAL          = 200;               % no. of Primal snapshots

NUM_ADJOINT = NUM_PRIMAL;


%NUM_STEP = 10;                              % time step

NUM_STEP_A = NUM_STEP;
TIME = 4001;

Xt = zeros(NUM_SYS, TIME);
Yt = zeros(NUM_SYS, TIME);

w = 1* randn(NUM_IN,TIME);
for i = 1:TIME
    Xt(:,i+1) = A * Xt(:,i)+ B * w(:,i);
end

v = 1* randn(NUM_OUT, TIME);
for i = 1:TIME
    Yt(:,i+1) = A' * Yt(:,i) + C' * v(:,i);
end

XST = Xt(:, 2000 : NUM_STEP:2000+ NUM_STEP*NUM_PRIMAL);
YST =  Yt(:,2000: NUM_STEP_A: 2000+ NUM_STEP_A*NUM_ADJOINT);
%%
Ht =YST' * XST;

%%
[Ut, St, Vt] = svd(Ht);

Unt = Ut(:, 1:NUM_NON);
Snt = St(1:NUM_NON, 1:NUM_NON);
Vnt = Vt(:,1:NUM_NON);

Trt= XST* Vnt* Snt^(-1/2);
Tlt = Snt^(-1/2) * Unt' * YST';


At = Tlt * A * Trt;
Bt = Tlt * B;
Ct = C *Trt;


VF = Trt*ct;
UF = inv(ct)*Tlt;

Cr = C * VF;
Br = UF * B;
Ar = UF *A * VF;