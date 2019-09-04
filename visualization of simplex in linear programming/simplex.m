function [dcxb,x,fval,exitflag,flag]=simplex(f,A,b,Aeq,beq)
%????????
%?1???????????MATLAB??????????
%?2?dcxb????????????flag=1???????????flag=-1?????????????????????????????
%?3???????dcxb?inf??????????inf???????????????????????1??
%?4?dcxb???????????????inf???inf????????
%?5?dcxb???????inf???????????

%?6?exitflag?????????x,fval??????????????????????????????????????
%?7?????????simplex?simplex???????????????????????????
%?8?flag=-1?dcxb????????inf????????????????????????????????(???????2?)
%?9??????linprog?????????
%?10?exitflag=1?????????exitflag=2???????exitflag=-1??????????exitflag=-2???????
%?11?flag=2???????????????????????????????????????????????????????????????????(3)
%?????1???????????????????????P44?23??1???
%??????????????
%f=[1 -2 1 -3]';
%A=[1 1 3 1;0 -2 1 1;0 -1 6 -1];
%b=[6 3 4]';
%[dcxb,x,fval,exitflag,flag]=simplex(f,A,b)
%?????2?:?????????????????????P23?1.8
%???????????
%f=[-0.75 20 -0.5 6 0 0 0]';
%Aeq=[0.25 -8 -1 9 1 0 0;0.5 -12 -0.5 3 0 1 0;0 0 1 0 0 0 1];
%A=[];
%b=[];
%beq=[0 0 1]';
%[dcxb,x,fval,exitflag,flag]=simplex(f,A,b,Aeq,beq)
%????(3):??????????????????P29?1.9
%???????????
%f=[1 -2]';
%A=[-1 -1;1 -1;0 1];
%b=[-2 -1 3]';
%[dcxb,x,fval,exitflag,flag]=simplex(f,A,b)
flag=1;
if nargin < 5, beq = [];
    if nargin < 4, Aeq = [];
    end
end
[A_m,A_n]=size(A);
[Aeq_m,Aeq_n]=size(Aeq);
bb=[b;beq];
m=A_m+Aeq_m;
n = max([length(f),A_n,Aeq_n]); % In case A is empty
if isempty(f), f=zeros(n,1); end%4
if isempty(A), A=zeros(0,n); end
if isempty(b), b=zeros(0,1); end
if isempty(Aeq), Aeq=zeros(0,n); end
if isempty(beq), beq=zeros(0,1); end
if(all(b>=0)==0)
    flag=3;
    b_len=length(b);
    k=1;
    for r=1:inf
        if(k<=b_len)
            if b(k)<0
                Aeq=[Aeq,zeros(Aeq_m,1);A(k,:),1];%1%5
                beq=[beq;b(k)];%2
                A=[A(1:k-1,:);A(k+1:A_m,:)];%3
                A=[A,zeros(A_m-1,1)];
                b=[b(1:k-1);b(k+1:A_m)];
                k=k-1;
                b_len=b_len-1;
                A_m=A_m-1;
                Aeq_m=Aeq_m+1;
                Aeq_n=Aeq_n+1;
                f=[f;0];
            else
                Aeq=[Aeq,zeros(Aeq_m,1);A(k,:),1];
                beq=[beq;b(k)];
                A=[A(1:k-1,:);A(k+1:A_m,:)];
                A=[A,zeros(A_m-1,1)];
                b=[b(1:k-1);b(k+1:A_m)];
                k=k-1;
                b_len=b_len-1;
                A_m=A_m-1;
                Aeq_m=Aeq_m+1;
                Aeq_n=Aeq_n+1;
                f=[f;0];
            end
            k=k+1;
        else break;
        end
    end
    A=[];
    b=[];
end
for k=1:Aeq_m
    if beq(k)<0
        Aeq(k,:)=-Aeq(k,:);
        beq(k)=-beq(k);
    end
end
[A_m,A_n]=size(A);
[Aeq_m,Aeq_n]=size(Aeq);
bb=[b;beq];
m=A_m+Aeq_m;
n = max([length(f),A_n,Aeq_n]); 
if isempty(f), f=zeros(n,1); end
if isempty(A), A=zeros(0,n); end
if isempty(b), b=zeros(0,1); end
if isempty(Aeq), Aeq=zeros(0,n); end
if isempty(beq), beq=zeros(0,1); end
AA=[A;Aeq];
pq=0;
PQ=zeros(n+m,1);
QP=zeros(n,1);
for k=1:n
    D=AA(:,k);
    [maxk,maxK]=max(D);
    [mink,minK]=min(D);
    if(maxk~=0)%?????????
        if(mink==0)
            D(maxK)=0;
            maxk=max(D);
            if(maxk==0)
                pq=pq+1;
                PQ(maxK)=k;
                QP(pq)=maxK;
                for s=1:pq-1
                    if(A(:,s)==A(:,pq))
                        pq=pq-1;
                    end
                end
            end
        end
    end
end
PQ=PQ(1:m);
QP=QP(1:pq);
W=eye(m);
for k=1:pq
    W(:,QP(k))=zeros(m,1);
end
WW=zeros(m,0);
kk=1;
for k=1:m
    if(W(:,k)==zeros(m,1))
    else
        for r=1:m
            if PQ(r)==0
                PQ(r)=n+kk;
                break;
            end
        end
        WW=[WW,W(:,k)];
        kk=kk+1;
    end
end
AA=[AA,WW];
Ab=[AA,bb];
FF=zeros(1,m+n-pq);
for r=n+1:m+n-pq
    for s=1:m
        if AA(s,r)~=0
            FF=FF+AA(s,:)/AA(s,r);
        end
    end
end
TT=f;
FF=[FF(1:n),zeros(1,m-pq)];
if flag==3
    ff=-FF';
    T=[f;zeros(m-pq,1)];
else
    ff=[f;zeros(m-pq,1)];
end
F=ff;
B=PQ;
BB=B;
dcxb=[inf,1:m+n-pq,inf;B,Ab;inf,ff',inf];
for loop=1:inf
    [minff,i]=min(ff);
    if(flag~=1)
        for k=1:i
            if(ff(k)<0)
                i=k;
                break;
            end
        end
    end
    C=inf*ones(m,1);
    if(minff>=0)
        exitflag=1;
        break;
    else
        if(AA(:,i)+abs(AA(:,i))==0)%????
            exitflag=-1;%exitflag??????
            dcxb=[dcxb;inf(1,m+n-pq)];
            break;
        else
            for h=1:m
                if(AA(h,i)>0)
                    C(h)=bb(h)/AA(h,i);
                end
            end
            [minC,j]=min(C);
            B(j)=i;
            Ab(j,:)=Ab(j,:)/Ab(j,i);
            ff=ff-AA(j,:)'*ff(i)/AA(j,i);
            if flag>1
                T=T-AA(j,:)'*T(i)/AA(j,i);
            end
            for k=1:m
                if(k~=j)%???<>????
                    Ab(k,:)=Ab(k,:)-Ab(j,:)*Ab(k,i)/Ab(j,i);
                end
            end
            bb=Ab(:,n+m-pq+1);
            AA=Ab(:,1:n+m-pq);
            dcxb=[dcxb;B,Ab;inf,ff',inf];
            BB=[BB,B];
        end
    end
    if flag==1
        for r=1:(loop)
            if(BB(:,r)==BB(:,loop+1))
                flag=-1;
                dcxb=[dcxb;inf*ones(1,m+n+2-pq)];
                break;
            end
        end
    end
end
if flag==3
    x=zeros(length(ff),1);
    for k=1:m
        x(B(k))=bb(k);
    end
    if x'*ff~=0
        exitflag=-2
    end
    fprintf('??????????\n??????????\n')
    dcxb
    flag=2;
    BB=[];
    ff=T(1:n);
end
if flag==2
    F=ff;
    AA=AA(:,1:n);
    Ab=[Ab(:,1:n),bb];
    BB=B;%%%%%
    dcxb=[inf,1:n,inf;B,Ab;inf,ff',inf];
    for loop=1:inf
        [minff,i]=min(ff);
        if(flag==-1)
            for k=1:i
                if(ff(k)<0)
                    i=k;
                    break;
                end
            end
        end
        C=inf*ones(m,1);
        if(minff>=0)
            exitflag=1;
            break;
        else
            if(AA(:,i)+abs(AA(:,i))==0)%????
                exitflag=-1;%exitflag??????
                dcxb=[dcxb;inf(1,m+n-pq)];
                break;
            else
                for h=1:m
                    if(AA(h,i)>0)
                        C(h)=bb(h)/AA(h,i);
                    end
                end
                [minC,j]=min(C);
                B(j)=i;
                Ab(j,:)=Ab(j,:)/Ab(j,i);
                ff=ff-AA(j,:)'*ff(i)/AA(j,i);
                for k=1:m
                    if(k~=j)%???<>????
                        Ab(k,:)=Ab(k,:)-Ab(j,:)*Ab(k,i)/Ab(j,i);
                    end
                end
                bb=Ab(:,n+1);
                AA=Ab(:,1:n);
                dcxb=[dcxb;B,Ab;inf,ff',inf];
                BB=[BB,B];
            end
        end
        if flag==1
            for r=1:(loop)
                if(BB(:,r)==BB(:,loop+1))
                    flag=-1;
                    dcxb=[dcxb;inf*ones(1,n+2)];
                    break;
                end
            end
        end
    end
end
if flag==2
    x=zeros(n,1);
    F=TT;
    fprintf('??????????\n')
else
    x=zeros(m+n-pq,1);
end
for k=1:m
    x(B(k))=bb(k);
end
if sum(ff==0)~=m
    exitflag=2;
end
fval=x'*F;
