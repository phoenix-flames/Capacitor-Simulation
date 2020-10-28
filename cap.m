%all of the constants
Lx = .5;
Ly = .5;
n = 50;
m = 50;
d = 0.1;
delx = 0.01;
dely = 0.01;
delV = 3;
epsilon = 8.854 * 10^(-12);
N = n * m;

%Voltage Column Vector
V1 = zeros(1,N);
for i = 1:N
V1(i) = 1.5;
end
V2 = zeros(1,N);
for i = 1:N
V2(i) = -1.5;
end
V = [V1 V2];
V = V';

Z = zeros(2*N);
xi = zeros(2*N,1);
xj = zeros(2*N,1);
%Z
for i = 1:2*N
    for j = 1:2*N
        if i == j
            Z(i,j) = (delx*asinh(1))/(pi*epsilon);
        else
            if ((i > N) && (j > N)) %both on bottom
                xi = [(mod((i-N-1),50)+1)*delx (floor((i-N-1)/50)+1)*dely];
                xj = [(mod((j-N-1),50)+1)*delx (floor((j-N-1)/50)+1)*dely];
            elseif ((i > N) && (j <= N)) %xi on bottom and xj on top
                xi = [(mod((i-N-1),50)+1)*delx (floor((i-N-1)/50)+1)*dely];
                xj = [(mod((j-1),50)+1)*delx (floor((j-1)/50)+1)*dely];
            elseif ((i <= N) && (j > N))%xi on top and xj on bottom
                xi = [(mod((i-1),50)+1)*delx (floor((i-1)/50)+1)*dely];
                xj = [(mod((j-N-1),50)+1)*delx (floor((j-N-1)/50)+1)*dely];
            elseif ((i <= N) && (j <= N))%both are on top plate
                xi = [(mod((i-1),50)+1)*delx (floor((i-1)/50)+1)*dely];
                xj = [(mod((j-1),50)+1)*delx (floor((j-1)/50)+1)*dely];
            end
       
            %if both vectors are on the same plate
            if (i <= N && j <= N)||(i > N && j > N)
                Z(i,j) = (1/(4*pi*epsilon))*((delx*dely)/norm(xi-xj));
            % if the vectors are on different plates
            else
                Z(i,j) = (1/(4*pi*epsilon))*((delx*dely)/sqrt((norm(xi-xj))^2+d^2));
            end
        end
    end
end

%ZV = Q
Q=linsolve(Z,V);

%changing Q from a vector to a matrix
%positive plate
k = 1;
Qpos = zeros(n);
for i = 1:n
    for j = 1:n
        Qpos(i,j) = Q(k);
        k = k+1;
    end
end
%negative plate
k = 2501;
Qneg = zeros(n);
for i = 1:n
    for j = 1:n
        Qneg(i,j) = Q(k);
        k = k + 1;
    end
end

%plot Q
X = 1:n;
Y = 1:m;
surf(X, Y, Qpos);
figure;
surf(X, Y, Qneg);
