function [] = truss_solver(m)%m is number of members
% member and joint related inputs are taken through the excel sheet
inp1 = xlsread('t_i_j.xlsx');
inp2 = xlsread('t_i_m.xlsx');
[a,b] = size(inp1);
for i = 1:a % for inputting joints data to a structured matrix
    for j=1:b
        if isnan(inp1(i,j))
            inp1(i,j)=0;
        end
    end
end
[a,b] = size(inp2);
for i = 1:a % for inputting members data to a structured matrix
    for j=1:b
        if isnan(inp2(i,j))
            inp2(i,j)=0;
        end
    end
end
dof = 2*length(inp1(1,:));% Total no of dof's
L = zeros(length(inp2(1,:)),1);
a = zeros(length(inp2(1,:)),1);
for i = 1:length(inp2(1,:)) % Finding length and alpha of memebers
    L(i) = sqrt((inp1(2,inp2(3,i))-inp1(2,(inp2(4,i))))^2+(inp1(3,inp2(3,i))-inp1(3,(inp2(4,i))))^2); %length
    y2 = inp1(3,inp2(4,i));
    y1 = inp1(3,inp2(3,i)); 
    x2 = inp1(2,inp2(4,i));
    x1 = inp1(2,inp2(3,i));
    a(i) = atan((y2-y1)/(x2-x1));
    if(a(i)<0)
        a(i)=pi+a(i); %angle
    end
end
f1 = @(a) [cos(a),sin(a),0,0;-1*sin(a),cos(a),0,0;0,0,cos(a),sin(a);0,0,-1*sin(a),cos(a)];
f2 = @(AE,L) (AE/L)*[1,0,-1,0;0,0,0,0;-1,0,1,0;0,0,0,0];
Sm = @(AE,L,a) (AE/L)*[(cos(a))^2, cos(a)*sin(a), -1*(cos(a))^2, -1*cos(a)*sin(a);...
                       cos(a)*sin(a), (sin(a))^2, -1*cos(a)*sin(a), -1*(sin(a))^2;...
                       -1*(cos(a))^2, -1*cos(a)*sin(a), (cos(a))^2, cos(a)*sin(a);...
                       -1*cos(a)*sin(a), -1*(sin(a))^2, cos(a)*sin(a), (sin(a))^2];
Ss = zeros(dof);
fea1 = zeros(m);
fea2 = zeros(m);% fixed end actions due to strain or heat
for i = 1:length(inp2(1,:)) % Finding Global stiffness matrix
    Ssi = zeros(dof);
    s = Sm(inp2(2,i),L(i),a(i));
    Ssi((2*inp2(3,i)-1):2*inp2(3,i),(2*inp2(3,i)-1):2*inp2(3,i)) = s(1:2,1:2);
    Ssi((2*inp2(4,i)-1):2*inp2(4,i),(2*inp2(4,i)-1):2*inp2(4,i)) = s(3:4,3:4);
    Ssi((2*inp2(3,i)-1):2*inp2(3,i),(2*inp2(4,i)-1):2*inp2(4,i)) = s(1:2,3:4);
    Ssi((2*inp2(4,i)-1):2*inp2(4,i),(2*inp2(3,i)-1):2*inp2(3,i)) = s(3:4,1:2);
    Ss = Ss + Ssi;
    
        fea1(i) = inp2(2,i)*inp2(7,i)*inp2(5,i); 
    
    
        fea2(i) = inp2(2,i)*inp2(6,i)/L(i); 
    
end
fea = fea1+fea2; %member end actions due to heat and prestrain
Rt = zeros(4,4,length(inp2(1,:)));
for i=1:m     % changing the member end actions to global structrual system
    Rt(:,:,i) = f1(a(i));
    ds = Rt(:,:,i)' * [fea(i);0;-fea(i);0];
    inp1(7,inp2(3,i)) = inp1(7,inp2(3,i))+ds(1);
    inp1(8,inp2(3,i)) = inp1(8,inp2(3,i))+ds(2);
    inp1(7,inp2(4,i)) = inp1(7,inp2(4,i))+ds(3);
    inp1(8,inp2(4,i)) = inp1(8,inp2(4,i))+ds(4);
end
j =1;l=1; % r for restrained, k for free
for i = 1:length(inp1(1,:)) % Finding Free dof's and restrained dof's
    if (inp1(6,i)==0)
        k(j) = 2*inp1(1,i)-1;
        j=j+1;
        k(j) = 2*inp1(1,i);
        j=j+1;
    elseif (inp1(6,i)==10)
        k(j) = 2*inp1(1,i);
        j = j+1;
        r(l) = 2*inp1(1,i)-1;
        l=l+1;
    elseif (inp1(6,i)==1)
        k(j) = 2*inp1(1,i)-1;
        j = j+1;
        r(l) = 2*inp1(1,i);
        l=l+1;
    elseif (inp1(6,i)==11)
        r(l) = 2*inp1(1,i)-1;
        l=l+1;
        r(l) = 2*inp1(1,i);
        l=l+1;
    end
end
Af = zeros(length(k),1);
Dr = zeros(dof-length(k),1);
Df = zeros(length(k),1);
Ar = zeros(dof-length(k),1);
Ssn1 = zeros(dof);l1 = 1;
for i=1:length(k) %row interchange
    if rem(k(i),2)~= 0
        Af(i) = inp1(7,((k(i)+1)/2));
    elseif rem(k(i),2) == 0
        Af(i) = inp1(8,(k(i)/2));
    end
    Ssn1(l1,:) = Ss(k(i),:);
    l1 = l1+1;
end
for i=1:dof-length(k)%row interchange
    if rem(r(i),2)~= 0
        Dr(i) = inp1(4,((r(i)+1)/2));
    elseif rem(r(i),2) == 0
        Dr(i) = inp1(5,(r(i)/2));
    end
    Ssn1(l1,:) = Ss(r(i),:);
    l1 = l1+1;
end
l1 = 1;Ssx = Ssn1;
for i=1:length(k)%column interchange
    Ssn1(:,l1) = Ssx(:,k(i));
    l1 = l1+1;
end
for i=1:dof-length(k)%column interchange
    Ssn1(:,l1) = Ssx(:,r(i));
    l1 = l1+1;
end
Ss11 = Ssn1(1:length(k),1:length(k));
Ss12 = Ssn1(1:length(k),(length(k)+1):dof);
Ss21 = Ssn1((length(k)+1):dof,1:length(k));
Ss22 = Ssn1((length(k)+1):dof,(length(k)+1):dof);
Df = inv(Ss11)*(Af-Ss12*Dr);
Ar = Ss21*Df + Ss22*(Dr);
for i=1:length(k) %printing displacements results
    if rem(k(i),2) == 0
       fprintf("The y displacement at node %d is %f\n",(k(i)/2),Df(i));
    elseif rem(k(i),2)~= 0
       fprintf("The x displacement at node %d is %f\n",((k(i)+1)/2),Df(i));
    end
end
for i=1:dof-length(k) %printing force results
    if rem(r(i),2) == 0
       fprintf("The y force at node %d is %f\n",(r(i)/2),Ar(i));
    elseif rem(r(i),2)~= 0
       fprintf("The x force at node %d is %f\n",((r(i)+1)/2),Ar(i));
    end
end
dsr = zeros(dof,1);
dsf = zeros(dof,1);
for i=1:length(k) %placing the free cordinates back in respective structural positions
    dsf(k(i)) = Df(i);
end
for i=1:length(r) %placing the restrained cordinates back in respective structural positions
    dsr(r(i)) = Dr(i);
end
Ds = dsf+dsr;
Dm = zeros(4,length(inp2(1,:)));
Am = zeros(4,length(inp2(1,:)));
Rt = zeros(4,4,length(inp2(1,:)));
Smi = zeros(4,4,length(inp2(1,:)));
for i = 1:length(inp2(1,:)) %printing the member forces and their nature(compression and tension)
    Rt(:,:,i) = f1(a(i));
    Smi(:,:,i) = f2(inp2(2,i),L(i)); 
    Dm(:,i) = Rt(:,:,i)*[Ds(2*inp2(3,i)-1);Ds(2*inp2(3,i));Ds(2*inp2(4,i)-1);Ds(2*inp2(4,i))];
    Am(:,i) = Smi(:,:,i)*Dm(:,i);
    if sign(Am(1,i)+fea(i)) == 1
        fprintf("The force and the nature of member %d is %f, %s respectively.\n",i,Am(1,i)+fea(i),'Compression');
    elseif sign(Am(1,i)+fea(i)) == -1
        fprintf("The force and the nature of member %d is %f, %s respectively.\n",i,-1*Am(1,i)+fea(i),'Tension');
    elseif Am(1,i)==0
         fprintf("The member %d is a zero force member.\n",i);
    end
end
end