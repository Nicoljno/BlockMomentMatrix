clear
clc

d=3;
for i=1:2*d
    rho(:,:,i)=RandomDensityMatrix(2*d);
end
id=eye(d);
for i=1:d
    target(:,:,i)=id(i,:)'*id(i,:);
end
for i=1:d
    target(:,:,d+i)=FourierMatrix(d)*target(:,:,i)*FourierMatrix(d)';
end
for i=1:2*d
    tmp=zeros(2*d,2*d);
    tmp(1:d,1:d)=target(:,:,i);
    target_(:,:,i)=tmp;
end

for k=1:d^2
    f(k,:)=toSeveralBases(k-1,[d d])+1;
end


arr=[];
for n_eps=31:51
    eps=0.001*(n_eps-1);
    for n_steps=1:20
        yalmip('clear')
        M=sdpvar(2*d,2*d,d,2,'hermitian','complex');
        G=sdpvar(2*d,2*d,d^2,'hermitian','complex');
        
        constr=[];
        constr=[constr, sum(G,3)==eye(2*d)]
        for i=1:d^2
            constr=[constr, G(:,:,i)>=0];
        end
        
        for a=1:d
            for x=1:2
                r=0;
                for lambda=1:d^2
                    pdet=[f(lambda,x)==a];
                    r=r+G(:,:,lambda)*pdet;
                end
                constr=[constr, M(:,:,a,x)==r];
            end
        end
        obj=0;
        for a=1:d
            for x=1:2
                obj=obj+trace(rho(:,:,a+(x-1)*d)*M(:,:,a,x));
            end
        end
%         obj=trace(rho(:,:,1)*M(:,:,1,1))+trace(rho(:,:,2)*M(:,:,2,1))+trace(rho(:,:,3)*M(:,:,3,1))+ ...
%             trace(rho(:,:,4)*M(:,:,1,2))+trace(rho(:,:,5)*M(:,:,2,2))+trace(rho(:,:,6)*M(:,:,3,2));
        
        ops=sdpsettings('solver','mosek');
        diagnostic=solvesdp(constr, -obj, ops);
        
        M=value(M);
        
        yalmip('clear')
        
        rho=sdpvar(2*d,2*d,d*2,'hermitian','complex');
        
        constr=[];
        for i=1:d*2
            constr=[constr, trace(rho(:,:,i))==1];
            constr=[constr, trace(rho(:,:,i)*target_(:,:,i))>=1-eps];
            constr=[constr, rho(:,:,i)>=0];
        end
        obj=0;
        for a=1:d
            for x=1:2
                obj=obj+trace(rho(:,:,a+(x-1)*d)*M(:,:,a,x));
            end
        end
%         obj=trace(rho(:,:,1)*M(:,:,1,1))+trace(rho(:,:,2)*M(:,:,2,1))+trace(rho(:,:,3)*M(:,:,3,1))+ ...
%             trace(rho(:,:,4)*M(:,:,1,2))+trace(rho(:,:,5)*M(:,:,2,2))+trace(rho(:,:,6)*M(:,:,3,2));
        
        ops=sdpsettings('solver','mosek');
        diagnostic=solvesdp(constr, -obj, ops);
        
        rho=value(rho);
    end
    
    arr=[arr; eps value(obj)/d];
end


