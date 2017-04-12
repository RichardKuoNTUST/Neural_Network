clear all;clc;close all;
read
readtrainoutput
readtrain
readtestdata
readtem
for i=1:21
    te(:,i)=te(i)/40;
end
k=0;
for i=1:21
    for j=1:6
        k=k+1;
        inputv(j,i)=volume3134(k)/500;
        outputv(j,i)=output3156(k)/500;
        inputww1(j,i)=week3(k);
        inputww2(j,i)=week2(k);
        inputww3(j,i)=week1(k);
    end
end
k=0;
for i=1:7
    for j=1:6
        k=k+1;
        testdata(j,i)=test3134(k)/500;
    end
end
ww=[0 1 0 1 0 1 1;
    1 1 0 0 1 1 0;
    0 0 1 1 1 1 0];
test=[testdata;ww]';
k=0;
% for i=1:21
%     for j=1:6
%         k=k+1;
%         inputw(j,i)=week(k)/10;
%         %outputw(j,i)=week1(k)/10;
%     end
% end

input=[inputv;inputww1(1,:);inputww2(1,:);inputww3(1,:)]';
output=outputv';
numberoutput=size(output);
ls = [9, 55,55, 6];
n =size(ls);
a=0.5;
lr=1.1;
weight1=randn(ls(1),ls(2))*0.1;
weight2=randn(ls(2),ls(3))*0.1;
weight3=randn(ls(3),ls(4))*0.1;

bais1=randn(1,ls(2))*0.1;
bais2=randn(1,ls(3))*0.1;
bais3=randn(1,ls(4))*0.1;

output1=zeros(1,ls(1));
output2=zeros(1,ls(2));
output3=zeros(1,ls(3));
output4=zeros(1,ls(4));

delta1=zeros(1,ls(2));
delta2=zeros(1,ls(3));
delta3=zeros(1,ls(4));

ni=size(input);

oldweight1=weight1;%zeros(ls(1),ls(2));
oldweight2=weight2;%zeros(ls(2),ls(3));
oldweight3=weight3;%zeros(ls(3),ls(4));

oldbais1=bais1;%zeros(1,ls(2));
oldbais2=bais2;%zeros(1,ls(3));
oldbais3=bais3;%zeros(1,ls(4));
for c=1:50
    for i=1:ni(1)
        t=output(i,:);
        output1=input(i,:);
        z2=output1*weight1-bais1;
        for j=1:ls(2)
            %output2(j)=1/(1+exp(-z2(j)));
            output2(j)=tanh(z2(j));
        end
        z3=output2*weight2-bais2;
        for j=1:ls(3)
            %output3(j)=1/(1+exp(-z3(j)));
            output3(j)=tanh(z3(j));
        end
        z4=output3*weight3-bais3;
        for j=1:ls(4)
            if z4(j)>1
                z4(j)=1;
            end
            if z4(j)<-1
                z4(j)=-1;
            end
            output4(j)=z4(j);
            %output4(j)=-1+2/(1+exp(-z4(j)));
        end
        delta3=(t-output4);%((1-output4.*output4))*(t-output4);
        delta2=(1-output3.*output3).*((weight3*delta3')');%output3.*(1-output3).*((weight3*delta3')');
        delta1=(1-output2.*output2).*((weight2*delta2')');%output2.*(1-output2).*((weight2*delta2')');
        error(i).ee(c,:)=abs(t-output4)/t;
        for j=1:ls(1)
            for k=1:ls(2)
                weight1(j,k)=weight1(j,k)+a*output1(j)*delta1(k);%+0.5*(weight1(j,k)-oldweight1(j,k));
            end
        end
        oldweight1=weight1;
        for j=1:ls(2)
            for k=1:ls(3)
                weight2(j,k)=weight2(j,k)+a*output2(j)*delta2(k);%+0.6*(weight2(j,k)-oldweight2(j,k));
            end
        end
        oldweight2=weight2;
        for j=1:ls(3)
            for k=1:ls(4)
                weight3(j,k)=weight3(j,k)+a*output3(j)*delta3(k);%+0.8*(weight3(j,k)-oldweight3(j,k));
            end
        end
        oldweight3=weight3;
        bais1=bais1-a*delta1;%-0.5*(bais1-oldbais1);
        oldbais1=bais1;
        bais2=bais2-a*delta2;%-0.6*(bais2-oldbais2);
        oldbais2=bais2;
        bais3=bais3-a*delta3;%-0.8*(bais3-oldbais3);
        oldbais3=bais3;
    end
    figure(1)
    hold on;
    q=1:c;
    for i=1:numberoutput(1)
        e(i)=mean(error(i).ee(c,:));
    end
    msm(c)=mean(e);
    plot(q,msm);
    a=a/(lr*c);
end

numbertest=size(test);
for i=1:numbertest(1)
    output1=test(i,:);
    testz2=output1*weight1-bais1;
    for j=1:ls(2)
        %output2(j)=1/(1+exp(-z2(j)));
        output2(j)=tanh(testz2(j));
    end
    testz3=output2*weight2-bais2;
    for j=1:ls(3)
        %output3(j)=1/(1+exp(-z3(j)));
        output3(j)=tanh(testz3(j));
    end
    testz4=output3*weight3-bais3;
    for j=1:ls(4)
        if testz4(j)>1
            testz4(j)=1;
        end
        if testz4(j)<-1
            testz4(j)=-1;
        end
        output4(j)=testz4(j);
        %output4(j)=-1+2/(1+exp(-z4(j)));
    end
    testoutput(i,:)=output4*500;
end
