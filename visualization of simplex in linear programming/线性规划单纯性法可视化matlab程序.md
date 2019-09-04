```matlab
function [dcxb,x,fval,exitflag,flag]=simplex(f,A,b,Aeq,beq)
%本程序相关说明：
%（1）本函数基本使用方法与MATLAB自身函数使用方法类似
%（2）dcxb表示整个过程的单纯形表，flag=1表示没有出现退化情况，flag=-1表示出现退化情况，改为依照勃兰德法则确定进基变量和出基变量
%（3）程序运行结果dcxb中inf无意义，表的第一行两inf中间的数字表示变量的下标，本程序默认变量下标从1开始
%（4）dcxb其余各行如果仅仅行首和行尾出现inf，则两inf中间内容为检验数
%（5）dcxb第一列列相邻两inf间的内容为基变量的下标
%（6）exitflag为负表示无解，此时x,fval并不是相应的基可行解和最有值，仅仅代表程序中止时的基解和该基解所对应的函数值
%（7）函数名暂时命名为simplex，simplex是谷歌翻译给出的单纯形法的英文翻译，不知是否为专业术语
%（8）flag=-1时dcxb中会出现一行全为inf，无实际意义，仅表示该行一下依照勃兰德法则确定进基变量和出基变量(参见用法示例（2）)
%（9）其他参数与linprog中相应参数意义相似
%（10）exitflag=1表示有唯一最优解，exitflag=2表示有无穷解，exitflag=-1表示目标函数无下界，exitflag=-2表示可行域为空
%（11）flag=2表示用到了两阶段法，此时程序运行结果将出现“第一阶段单纯形表为”“第二阶段单纯形表为”字样，并显示两个阶段各自的单纯形表，详见用法示例(3)
%用法示例（1）：题目为束金龙编的《线性规划理论与模型应用》P44第23题（1）小题
%在命令窗口依次输入以下语句：
%f=[1 -2 1 -3]';
%A=[1 1 3 1;0 -2 1 1;0 -1 6 -1];
%b=[6 3 4]';
%[dcxb,x,fval,exitflag,flag]=simplex(f,A,b)
%用法示例（2）:题目为束金龙编的《线性规划理论与模型应用》P23例1.8
%命令窗口输入的命令为：
%f=[-0.75 20 -0.5 6 0 0 0]';
%Aeq=[0.25 -8 -1 9 1 0 0;0.5 -12 -0.5 3 0 1 0;0 0 1 0 0 0 1];
%A=[];
%b=[];
%beq=[0 0 1]';
%[dcxb,x,fval,exitflag,flag]=simplex(f,A,b,Aeq,beq)
%用法示例(3):题目为课本《线性规划理论与模型应用》P29例1.9
%命令窗口输入的命令为：
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
    if(maxk~=0)%这种情况是不正确的
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
        if(AA(:,i)+abs(AA(:,i))==0)%注意等号
            exitflag=-1;%exitflag为负表示无解
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
                if(k~=j)%不等号<>不对吗？
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
    fprintf('本题采用两阶段法求解\n第一阶段单纯形表为：\n')
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
            if(AA(:,i)+abs(AA(:,i))==0)%注意等号
                exitflag=-1;%exitflag为负表示无解
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
                    if(k~=j)%不等号<>不对吗？
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
    fprintf('第二阶段单纯形表为：\n')
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

```

```
题目为课本《线性规划理论与模型应用》P29例1.9
```

![simplex_run中的problem_1](/Users/mac/Desktop/研一课程/最优化理论基础/课程相关matlab程序/线性规划的单纯性法的可视化/simplex_run中的problem_1.jpg)

```
<线性规划理论与模型应用>书籍在线阅读网站
http://img.sslibrary.com/n/slib/book/slib/11162405/8b529a0100524429a5183d5fb36217b6/945c0c653998c17724eef6a2be46b9f1.shtml?dxbaoku=false&deptid=1071&fav=http%3A%2F%2Fwww.sslibrary.com%2Freader%2Fpdg%2Fpdgreader%3Fd%3D501c6455d1552b2543a2036e7e089d2c%26ssid%3D11162405&fenlei=1301140101&spage=1&t=5&username=202.120.19.10&view=-1
```

