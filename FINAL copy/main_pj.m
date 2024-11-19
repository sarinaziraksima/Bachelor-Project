%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auhter: Sarina Ziraksima %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [errorr,final_results]=main_pj(gain,r_c_val,fn,Fpass,Apass,Fstop,Astop,methodSolution,ftype,N,topology_type,n,FitnessLimit_Data,MaxGenerations_Data,opamps)
%format shortG
errorr = "There is no error";

MaxStallGenerations_Data = MaxGenerations_Data;


if (r_c_val==1) % E24 standard
    x4 = [1 1e1 1e2 1e3 1e4];
    x7 = [1 1e1 1e2 1e3 1e4 1e5 1e6 1e7];
    x8 = [1 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11];
    x11 = [1 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11];
    x10 = [1 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10];


    c = sort((1e-12)*[1.0*x11 1.1*x4 1.2*x4 1.3*x4 1.5*x10 1.6*x4 1.8*x4 2.0*x4 2.2*x10 2.4*x4 2.7*x4 3.0*x4 3.3*x10 3.6*x4 3.9*x4 4.3*x4 4.7*x10 5.1*x4 5.6*x4 6.2*x4 6.8*x10 7.5*x4 8.2*x4 9.1*x4]);
    r = sort([1.0*x8 1.1*x7 1.2*x7 1.3*x7 1.5*x7 1.6*x7 1.8*x7 2.0*x7 2.2*x7 2.4*x7 2.7*x7 3.0*x7 3.3*x7 3.6*x7 3.9*x7 4.3*x7 4.7*x7 5.1*x7 5.6*x7 6.2*x7 6.8*x7 7.5*x7 8.2*x7 9.1*x7]);
    
elseif (r_c_val==0) % customize
    v1 = [1 1e2 1e4 1e5 1e6 1e7 1e8 1e9];
    v2 = [1e1 1e3 1e4 1e5 1e6 1e7 1e8 1e9];
    v3 = [1e1 1e2 1e4 1e5 1e6 1e7 1e8 1e9];
    v4 = [1 1e2 1e4 1e5 1e6 1e7 1e8 1e9];%v1
    v5 = [1e1 1e2 1e4 1e5 1e6 1e7 1e8];
    v6 = [1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9];
    v7 = [1 1e1 1e2 1e3 1e4 1e5 1e7 1e8];
    v8 = [1 1e1 1e2 1e3 1e4 1e5 1e7 1e8];

    v=[1 1e1 1e2 1e3 1e4 1e5 1e6 1e7];

    c = sort((1e-12)*   [1*v1 1.5*v2 1.8*v3 2.2*v4 3.3*v5 4.7*v6 6.8*v7 8.2*v8]   );
    r = sort(   [1*v 1.2*v 1.5*v 1.8*v 2.2*v 2.7*v 3.3*v 4.7*v 5.6*v 6.8*v 8.2*v]   );
end

Ri = opamps(1);
Ro = opamps(2);
op_wu = opamps(3);


%select C & R

c=c(c>4.7e-10 & c<470e-6);
r=r(r>Ro*10 & r<Ri/10);
rm=r;




if (1==1)
  
    % Start with the default options
    options = optimoptions('ga');
    options = optimoptions(options,'MutationFcn', {  @mutationuniform 0.42308 });
    %options = optimoptions(options,'MutationFcn',{@mutationgaussian, scale, shrink});
    %options = optimoptions(options,'MutationFcn', @mutationadaptfeasible);

    options = optimoptions(options,'SelectionFcn',{@selectionroulette});
    options = optimoptions(options,'MaxGenerations', MaxGenerations_Data);
    options = optimoptions(options,'FitnessLimit', FitnessLimit_Data);
    options = optimoptions(options,'MaxStallGenerations', MaxStallGenerations_Data);
    options = optimoptions(options,'PlotFcn', @gaplotbestf);

    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (methodSolution(1)=='e') % set elliptic and chebyshev2 as zeross topology
    topology_type = 'zeross';
elseif (methodSolution(1)=='c')
    if (methodSolution(6)=='2')
        topology_type = 'zeross';      
    end
end

%%%%%%%%%%%%%%%%%%%%% main

[b,a,Q,w0,N,Gb,p,z,k] = make_Q_w0(Fpass,Fstop,Apass,Astop,ftype,methodSolution,N,fn)




if(ftype(1)=='l'||ftype(1)=='h')
    final_results=zeros(round(N/2),15)
    if mod(N,2)==1

        [x_final_f,exitflag_f,fval_f,final_results] = order1design(final_results,options,r,c,w0,op_wu,n)

    end


    for(i=2:2:N)


        j=i/2;
        if(mod(N,2)==1)
            j=j+1;  
        end
          
        [x_final,x,fval,final_results] = order2design(ftype,topology_type,Q,w0,n,r,c,options,i,j,op_wu,N,final_results,a,b)

    end

    
    


end






if(ftype(1)=='b')

    final_results=zeros(round(N),15)
    
    if(mod(N,2)==1)
        w0(N)=sqrt(abs(p(N))*abs(p(N+1)));
        Q(N)=w0(N)/(abs(p(N))+abs(p(N+1)));
        Q(N+1) = Q(N);
        w0(N+1) = w0(N);

    end

    w0=w0(1:2:end)
    Q=Q(1:2:end)


    w0l=w0(round(N/2)+1:N)
    Ql=Q(round(N/2)+1:N)
    if(mod(N,2)==1)
        w0b=w0(round(N/2))
        Qb=Q(round(N/2))
    end
    w0h=w0(1:floor(N/2))
    Qh=Q(1:floor(N/2))





    for(i=1:length(Qh))
        ftype='high'

        j=i;
        [x_final,x,fval,final_results] = order2design(ftype,topology_type,Qh,w0h,n,r,c,options,i,j,op_wu,N,final_results,a,b)

    end

    if (mod(N,2)==1)
        
        j=1+length(Qh);

        [x_final,x,fval,final_results] = order2design(ftype,topology_type,Qb,w0b,n,r,c,options,i,j,op_wu,N,final_results,a,b)


    end

    for(i=1:length(Ql))
        ftype='low'
        j=i+length(Qh);
        if (mod(N,2)==1)
            j=j+1;
        end

        [x_final,x,fval,final_results] = order2design(ftype,topology_type,Ql,w0l,n,r,c,options,i,j,op_wu,N,final_results,a,b)

    end



end







%if (topology_type(1)=='z')
%    gain(1,:)=final_results(:,11)./final_results(:,9)
%end























if(ftype(1)~='s')
    final_results
    errorr
   
  
    
end


    


%%%%%%%%%%%%%%%%%%%%% first order filter designer
function     [x_final_f,exitflag_f,fval_f,final_results] = order1design(final_results,options,r,c,w0,op_wu,n)




     if(w0(1)>op_wu/1)    
        errorr = 'An error apeared, selected opamp frequency response is invalide for practical cases!';
    end
    
    nvars=2; % number of undifiend R and C
    design = make_integrator_RC_ff(w0(1),r,c)
    ll = ones(1,nvars);
    ml = [ones(1,nvars-1)*length(r) ones(1,1)*length(c)];

    for i=1:n

        [x_f(i,1:nvars),fval_f(i),exitflag_f(i),output_f,population_f,scores_f] = ga(design,nvars,[],[],[],[],ll,ml,[],[],options);

    end

    
    

    x_final_f=(x_f(fval_f==min(fval_f),1:nvars))
    x_final_f=x_final_f(1,:);
    final_results(1,1) = r(round(x_final_f(1)));
    final_results(1,2) = c(round(x_final_f(2)));
    
  
    
    
end
%%%%%%%%%%%%%%%%%%%%% second order filter designer
                                         
function [x_final,x,fval,final_results] = order2design(ftype,topology_type,Q,w0,n,r,c,options,i,j,op_wu,N,final_results,a,b)
    if(w0(i)>op_wu/gain)    
        errorr = 'An error apeared, selected opamp frequency response is invalide for practical cases!';
    end
    [x_final,x,fval] = main_optimizer_for_senthisizing_circuit(ftype,topology_type,Q(i),w0(i),n,r,c,options,j,a,b);
   
    
    
   
    %final_results(j,1:length(x_final)) = x_final;
    index_c=2;
    if (ftype(1)=='h')
        if (topology_type(1)=='M'||topology_type(1)=='a')
            index_c=3;
        elseif (topology_type(1)=='t')               
                if(topology_type(8)=='2')
                    index_c=3;
                end
        end
    end  
    if (topology_type(1)=='z')
        index_c=2;
    end 
    index_c=index_c-1;
    for(i=1:length(x_final))
        
        if(i<length(x_final)-index_c)
            final_results(j,i) = r(round(x_final(i)));
        else
            final_results(j,i) = c(round(x_final(i)));
            
        end
    end
    
end
%%%%%%%%%%%%%%%%%%%%% generating circuit functions
function [x_final,x,fval] = main_optimizer_for_senthisizing_circuit(ftype,topology_type,Q,w0,n,r,c,options,j,a,b)


if (topology_type(1)=='s')
    if gain ~= 1
        errorr = "An error apeared, this topology only supports unity gain!";
    end
    nvars=4; % number of undifiend R and C
    if(ftype(1)=='b')
        nvars=5;
    end
    design = make_sallenkey_stn_ff(ftype,Q,w0,r,c);
    ll = ones(1,nvars);
    ml = [ones(1,nvars-2)*length(r) ones(1,2)*length(c)];
end 
if (topology_type(1)=='M')
    if gain ~= 1
        errorr = "An error apeared, this topology only supports unity gain!";
    end
    nvars=5; % number of undifiend R and C
    design = make_MFB_stn_ff(ftype,Q,w0,r,c);
    c_index=2;
    if(ftype(1)=='h')
        c_index=3;
    end
    ll = ones(1,nvars);
    ml = [ones(1,nvars-c_index)*length(r) ones(1,c_index)*length(c)];
end 
if (topology_type(1)=='a')
    nvars=8; % number of undifiend R and C
    design = make_akerberg_stn_ff(ftype,Q,w0,r,c,gain);
    c_index=2;
    if(ftype(1)=='h')
        c_index=3;
    end
    ll = ones(1,nvars);
    ml = [ones(1,nvars-c_index)*length(r) ones(1,c_index)*length(c)];
end 
if (topology_type(1)=='t')
    if(topology_type(8)=='1')
        nvars=8; % number of undifiend R and C
        if(ftype(1)=='h')
            nvars=9;
        end
        design = make_thomas1_stn_ff(ftype,Q,w0,r,c,gain);
        ll = ones(1,nvars);
        ml = [ones(1,nvars-2)*length(r) ones(1,2)*length(c)];
    end
end 
if (topology_type(1)=='t')
    if(topology_type(8)=='2')
        nvars=8; % number of undifiend R and C
        c_index=2;
        if(ftype(1)=='h')
            c_index=3;
        end
        design = make_thomas2_stn_ff(ftype,Q,w0,r,c,gain);
        ll = ones(1,nvars);
        ml = [ones(1,nvars-c_index)*length(r) ones(1,c_index)*length(c)];
    end
end 



if (topology_type(1)=='z')
    if gain ~= 1
        errorr = "An error apeared, this topology only supports unity gain!";
    end
    nvars=15; % number of undifiend R and C
    design = zeros_ff(ftype,r,rm,c,a,b);
    ll = ones(1,nvars);
    ml = [ones(1,nvars-2)*length(r) ones(1,2)*length(c)];
    
end 

x=zeros(n,nvars); % output array
%fval=zeros(n,1);
%exitflag=zeros(n,1);
%output=zeros(n,1);

for i=1:n
    
    [x(i,1:nvars),fval(i),exitflag(i),output,population,scores] = ga(design,nvars,[],[],[],[],ll,ml,[],[],options);
    
   

end


%%% documantation
% R is k-ohm    c is n-F





x_final=(x(fval==min(fval),1:nvars));
    
 



x_final=x_final(1,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%% calculate fc and get Q and w0
function [b,a,Q,w0,N,Gb,p,z,k] = make_Q_w0(Fpass,Fstop,Apass,Astop,ftype,methodSolution,N,fn)


z=0;
k=0;
p=0;


    

if (methodSolution(1) == 'b')

    if(N==0)
        % Design Butterworth low-pass filter
        [N,wn] = buttord(2*pi*Fpass,2*pi*Fstop,Apass,Astop,'s')
        fn=wn/(2*pi);
    
    end

               
        
    [z,p,k] = butter(N,2*pi*fn,ftype,'s');
    [b,a] = zp2tf(z,p,k);
    
    

    % Get frequency response using freqs
    [H, W] = freqs(b,a,4096);
    
  
    figure;  
    plot(W/(2*pi),20*log10(gain*abs(H)),'Color',[0,0.7,0.9])
    title('Amplitude of frequency response')
    xlabel('f')
    ylabel('|H(j2.pi.f)|')

elseif(methodSolution(1)=='e' )
    
    % Define desired filter characteristics

    
    if(N==0)
        % Design Butterworth low-pass filter
        [N,~] = ellipord(2*pi*Fpass,2*pi*Fstop,Apass,Astop,'s')
    end
    
    % Design ellip low-pass filter
    [b, a] = ellip(N, Apass, Astop, 2*pi*Fpass, ftype, 's');

   
    
   

    % Get frequency response using freqs
    [H, W] = freqs(b,a,4096);
    
  
    figure;
    plot(W/(2*pi),20*log10(gain*abs(H)))
    title('Amplitude of frequency response')
    xlabel('f')
    ylabel('|H(j2.pi.f)|')

elseif((methodSolution(1)=='c') & (methodSolution(6)=='1'))

    % Define desired filter characteristics
    if(N==0)
        % Design Butterworth low-pass filter
        [N,~] = cheb1ord(2*pi*Fpass,2*pi*Fstop,Apass,Astop,'s')
    end
    

    % Design Chebyshev Type I low-pass filter
    [b, a] = cheby1(N, Apass, 2*pi*Fpass, ftype,'s');

    % Get frequency response using freqs
    [H, W] = freqs(b, a, 4096);
 
 
  
    
    figure;
    plot(W/(2*pi),20*log10(gain*abs(H)))
    title('Amplitude of frequency response')
    xlabel('f')
    ylabel('|H(j2.pi.f)|')

    
elseif((methodSolution(1)=='c') & (methodSolution(6)=='2'))
    % Define desired filter characteristics
    
    if(N==0)
        % Design Butterworth low-pass filter
        [N,~] = cheb2ord(2*pi*Fpass,2*pi*Fstop,Apass,Astop,'s')
    end
   
    

    % Design Chebyshev Type I low-pass filter
    [b, a] = cheby2(N, Astop, 2*pi*Fstop, ftype,'s');

    % Get frequency response using freqs
    [H, W] = freqs(b, a, 4096);
 
 
  
    
    figure;
    plot(W/(2*pi),20*log10(gain*abs(H)))
    title('Amplitude of frequency response')
    xlabel('f')
    ylabel('|H(j2.pi.f)|')


end






%wn = fn*2*pi;

%[n,d] = butter(N,wn,'s')

Gb = tf(b,a)

[w0, zeta, p] = damp(Gb);

Q = 1./(2*zeta);
filt_params = table(Q, w0);



display(filt_params);


if (ftype(1)=='b')
    gain = gain^(1/N)
else
    gain = gain^floor(N/2)
end



end

%%%%%%%%%%%%%%%%%%%%% different filters fitness function
function p = make_sallenkey_stn_ff(ftype,Q,w0,r,c)
p = @sallenkey_stn_ff;

   




    function output = sallenkey_stn_ff(input)


    





 
    %%% documantation
    % R is k-ohm    c is n-F


    %%%%%%%%%%%%%%%%%%%%% formula functions
    %Rmax = r(end);
    %Rmin = r(1);
    %Cmax = c(end);
    %Cmin = c(1);
    if(ftype(1)=='l')
        r1 = r(round(input(1)));%r1
        r2 = r(round(input(2)));%r2
        c1 = c(round(input(3)));%c1
        c2 = c(round(input(4)));%c2
        f1 = abs((abs(r1)+abs(r2))/(abs(c1)*abs(r1)*abs(r2))-(w0/Q))/(w0/Q);
        f2 = abs((((abs(r1))*(abs(r2))*(abs(c1))*(abs(c2)))^-1)-(w0)^2)/((w0)^2);
        f3 = (f1+f2)*0.5; 
        %f3 = f3^2;
        output=f3;
    end
    
    if(ftype(1)=='b')
        r1 = r(round(input(1)));%r1
        r2 = r(round(input(2)));%r2
        r5 = r(round(input(3)));%r5
        c1 = c(round(input(4)));%c1
        c2 = c(round(input(5)));%c2
        k=1;
        tt=(   ((abs(c1)*abs(r1))^-1)  +   ((abs(c1)*abs(r2))^-1)   +   (((abs(c1)*abs(r5))^-1))   +   (((abs(c2)*abs(r2))^-1))    );
        f1 = abs(tt-(w0/Q))/(w0/Q);
        f2 = abs((abs(r1)+abs(r5))*((((abs(r1))*(abs(r2))*(abs(c1))*(abs(c2)))^-1))-(w0)^2)/((w0)^2);
        f3= abs(  (1/(tt*0.5*abs(r1)*abs(c1)))  -  k)/k;
        f4 = (f1+f2)*0.45+f3*0.1; 
        %f3 = f3^2;
        output=f4;
        
    end
    
    if(ftype(1)=='h')
        r1 = r(round(input(1)));%r1
        r2 = r(round(input(2)));%r2
        c1 = c(round(input(3)));%c1
        c2 = c(round(input(4)));%c2
        f1 = abs((abs(c1)+abs(c2))/(abs(r2)*abs(c1)*abs(c2))-(w0/Q))/(w0/Q);
        f2 = abs((((abs(r1))*(abs(r2))*(abs(c1))*(abs(c2)))^-1)-(w0)^2)/((w0)^2);
        f3 = (f1+f2)*0.5; 
        %f3 = f3^2;
        output=f3;
    end
    
    
    
    
 

 
    end
end
function p = make_MFB_stn_ff(ftype,Q,w0,r,c)
p = @MFB_stn_ff;

   


    function output = MFB_stn_ff(input)

    %input=x_final
    
    %[r1 r2 r3 c1 c2]



    %%%%%%%%%%%%%%%%%%%%% formula functions
    
    
    
    
    if(ftype(1)=='l')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        c1 = c(round(input(:,4)));
        c2 = c(round(input(:,5)));
        
        %k=r1/r3
        k=gain;
        f1 = abs((abs(r1)*abs(r2)*abs(r3)*abs(c1))/(abs(r1)*abs(r2)+abs(r2)*abs(r3)+abs(r3)*abs(r1))-(Q/w0))/(Q/w0);
        f2 = abs((abs(r1))*(abs(r2))*(abs(c1))*(abs(c2))-(w0)^-2)/(w0)^-2;
        f3 = abs(((abs(r1))/(abs(r3))) - k)/k;
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
    if(ftype(1)=='b')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        c1 = c(round(input(:,4)));
        c2 = c(round(input(:,5)));
        
        k=gain;
        f1 = abs((abs(r2)*abs(c1)*abs(c2))/(abs(c1)+abs(c2))-(Q/w0))/(Q/w0);
        f2 = abs((abs(c1)^-1)*(abs(c2)^-1)*(abs(r2)^-1)*((abs(r1)^-1)+(abs(r3)^-1))-(w0)^2)/(w0)^2;
        f3 = abs(      (abs(r2)*abs(c2)*(abs(r3)^-1))/(abs(c1)+abs(c2))     - k)/k;
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
   
    if(ftype(1)=='h')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        c1 = c(round(input(:,3)));
        c2 = c(round(input(:,4)));
        c3 = c(round(input(:,5)));
        
        %k=c3/c1
        k=gain;

        f1 = abs((abs(r2)*abs(c1)*abs(c2))/(abs(c1)+abs(c2)+abs(c3))-(Q/w0))/(Q/w0);
        f2 = abs((abs(c1))*(abs(c2))*(abs(r2))*(abs(r1))-(w0)^-2)/(w0)^-2;
        f3 = abs(   (abs(c3)/abs(c1))   - k)/k;
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end

    output=f4;
    end
end
function p = make_akerberg_stn_ff(ftype,Q,w0,r,c,gain)
p = @akerberg_stn_ff;

   


    function output = akerberg_stn_ff(input)
  


    if(ftype(1)=='l')

        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        c1 = c(round(input(:,7)));
        c2 = c(round(input(:,8)));
        %%%%%%%%%%%%%%%%%%%%% formula functions
        %k=r6r3/r1r5
        k=gain;
        f1 = abs(((abs(c1)*abs(r2))^-1)-(w0/Q))/(w0/Q);
        f2 = abs(abs(r6)/(abs(r3)*abs(r4)*abs(r5)*abs(c1)*abs(c2))-(w0)^2)/(w0)^2;
        f3 = abs((abs(r3)*abs(r5))/(abs(r1)*abs(r6))-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end

    if(ftype(1)=='b')

        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        c1 = c(round(input(:,7)));
        c2 = c(round(input(:,8)));
        %%%%%%%%%%%%%%%%%%%%% formula functions
        
        k=gain;
        f1 = abs(((abs(c1)*abs(r2))^-1)-(w0/Q))/(w0/Q);
        f2 = abs(abs(r6)/(abs(r3)*abs(r4)*abs(r5)*abs(c1)*abs(c2))-(w0)^2)/(w0)^2;
        f3 = abs(abs(r2)/abs(r1)-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
    if(ftype(1)=='h')

        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        c1 = c(round(input(:,6)));
        c2 = c(round(input(:,7)));
        c3 = c(round(input(:,8)));
        %%%%%%%%%%%%%%%%%%%%% formula functions
        
        k=gain;
        f1 = abs(((abs(c1)*abs(r2))^-1)-(w0/Q))/(w0/Q);
        f2 = abs(abs(r1)/(abs(r3)*abs(r4)*abs(r5)*abs(c1)*abs(c2))-(w0)^2)/(w0)^2;
        f3 = abs(abs(c3)*abs(r2)-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
   

    output=f4;
    end
end
function p = make_thomas1_stn_ff(ftype,Q,w0,r,c,gain)
p = @thomas1_stn_ff;

   


    function output = thomas1_stn_ff(input)

    if(ftype(1)=='l')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        c1 = c(round(input(:,7)));
        c2 = c(round(input(:,8)));

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=gain;
        f1 = abs(abs(c1)*abs(r3)-(Q/w0))/(Q/w0);
        f2 = abs(abs(r2)*abs(r4)*abs(r6)*abs(c1)*abs(c2)/abs(r5)-((w0)^-2))/((w0)^-2);
        f3 = abs(((abs(r6)*abs(r4))/(abs(r5)*abs(r1)))-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
        output=f4;
    end
    
    if(ftype(1)=='b')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        c1 = c(round(input(:,7)));
        c2 = c(round(input(:,8)));

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=gain;
        f1 = abs(abs(c1)*abs(r3)-(Q/w0))/(Q/w0);
        f2 = abs(abs(r2)*abs(r4)*abs(r6)*abs(c1)*abs(c2)/abs(r5)-((w0)^-2))/((w0)^-2);
        f3 = abs((abs(r3)*abs(r5))/(abs(r1)*abs(r6))-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
        output=f4;
    end
    if(ftype(1)=='h')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        r7 = r(round(input(:,7)));
        c1 = c(round(input(:,8)));
        c2 = c(round(input(:,9)));

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=gain;
        f1 = abs(((abs(c1)*abs(r3))^-1)-(w0/Q))/(w0/Q);
        f2 = abs(abs(r2)*abs(r4)*abs(r6)*abs(c1)*abs(c2)/abs(r5)-((w0)^-2))/((w0)^-2);
        f3 = abs((abs(r5)/abs(r1))-k);
        f4 = abs(   (abs(r7)^-1)   -   (abs(r6)/(abs(r1)*abs(r3)))   )*r5/(r6*c1);%%%%%%%%%%%%%%
        f5 = 0.3*f1+0.3*f2+f3*0.1+0.3*f4; 
        output=f5;
    end
 



    
    end
end
function p = make_thomas2_stn_ff(ftype,Q,w0,r,c,gain)
p = @thomas2_stn_ff;

   


    function output = thomas2_stn_ff(input)


    if(ftype(1)=='l')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        c1 = c(round(input(:,7)));
        c2 = c(round(input(:,8)));

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=gain;
        f1 = abs(abs(c1)*abs(r3)-(Q/w0))/(Q/w0);
        f2 = abs(abs(r2)*abs(r4)*abs(r5)*abs(c1)*abs(c2)/abs(r6)-((w0)^-2))/((w0)^-2);
        f3 = abs(((abs(r2)*abs(r5))/(abs(r6)*abs(r1)))-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
    
    if(ftype(1)=='b')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        c1 = c(round(input(:,7)));
        c2 = c(round(input(:,8)));

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=gain;
        f1 = abs(abs(c1)*abs(r3)-(Q/w0))/(Q/w0);
        f2 = abs(abs(r2)*abs(r4)*abs(r5)*abs(c1)*abs(c2)/abs(r6)-((w0)^-2))/((w0)^-2);
        f3 = abs(abs(r3)/abs(r1)-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
    if(ftype(1)=='h')
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        c1 = c(round(input(:,6)));
        c2 = c(round(input(:,7)));
        c3 = c(round(input(:,8)));

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=gain;
        f1 = abs(abs(c1)*abs(r3)-(Q/w0))/(Q/w0);
        f2 = abs(abs(r2)*abs(r4)*abs(r5)*abs(c1)*abs(c2)/abs(r1)-((w0)^-2))/((w0)^-2);
        f3 = abs((abs(c3)/abs(c1))-k);
        f4 = 0.45*f1+0.45*f2+0.1*f3; 
    end
    
    
    
    output=f4;
    end
end

function p = make_integrator_RC_ff(w0,r,c)
p = @integrator_RC_ff;



    function output = integrator_RC_ff(input)

    
    r1 = r(round(input(:,1)));
    c1 = c(round(input(:,2)));


    %%%%%%%%%%%%%%%%%%%%% formula functions

  
    f1 = abs(abs(c1)*abs(r1)-(w0)^-1)/(w0)^-1;
    output = f1;
   
    
    
    end
end




function p = zeros1_ff(ftype,r,c,a,b)
p = @zeros1_ff;

   


    function output = zeros1_ff(input)


    
        r1 = 82000;%r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        r7 = r(round(input(:,7)));
        r8 = r(round(input(:,8)));
        r9 = r(round(input(:,9)));
        r10 = r(round(input(:,10)));
        c1 = c(round(input(:,11)));
        c2 = c(round(input(:,12)));
        
        if (a(2)<0.001)
            a(2)=0.001;
        end
        if (a(3)<0.001)
            a(3)=0.001;
        end
        if (b(2)<0.001)
            b(2)=0.001;
        end
        if (b(3)<0.001)
            b(3)=0.001;
        end
    
            

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=1;
        f1 = abs(    (    abs(r9)*(abs(r7)+abs(r8))    )/(    abs(r1)*abs(r8)*abs(c1)*(abs(r9)+abs(r10))      )  -b(2)  )/(b(2));
        f2 = abs(    abs(r7)/(abs(r1)*abs(r2)*abs(r8)*abs(c1)*abs(c2))   -  b(3)  )/(b(3));
        
        A = abs(   abs(r4)*(abs(r5)+abs(r6))   /    (abs(r5)*(abs(r3)+abs(r4)) )  );
        
        B = abs(   A*abs(r3)/(abs(r1)*abs(r4)*abs(c1))    );
        D = abs(   abs(r6)/(abs(r5)*abs(r1)*abs(r2)*abs(c1)*abs(c2))    );
        
        f3 = abs(B-a(2))/a(2);
        f4 = abs(D-a(3))/a(3);
        
        %k1 = abs(r4)*abs(r8)*(  abs(r9)+abs(r10)  )*(  abs(r5)+abs(r6)   );
        %k2 = abs(r5)*abs(r10)*(  abs(r7)+abs(r8)  )*(  abs(r3)+abs(r4)   );
        %f5 = abs((k1/k2)-k);
        f7=abs(r1*c1-0.82)/0.82;
        f6 = f2*0.13+f3*0.33+f4*0.43+f7*0.10;
        %f6 = 0.22*f2+0.22*f3+0.22*f4+0.12*f1f7*10;

    
    
    output=f6;
    end
end



function p = zeros2_ff(ftype,r,c,a,b)
p = @zeros2_ff;

   


    function output = zeros2_ff(input)


    
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        r7 = r(round(input(:,7)));
        r8 = r(round(input(:,8)));       
        c1 = c(round(input(:,9)));
        c2 = c(round(input(:,10)));
        c3 = c(round(input(:,11)));
        if(b(2)<0.0001)
            b(2)=0.001;
        end
            

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=1;
        f1 = abs(         (1/abs(c3))*(    (1/abs(r4))   -    (  abs(r6)/(abs(r3)*abs(r8))  )    )              -b(2)  )/(b(2));
        f2 = abs(       abs(r6)/(abs(r3)*abs(r5)*abs(r7)*abs(c2)*abs(c3))       -b(3)  )/(b(3));
        f3 = abs(      1/(abs(r1)*abs(c1))                  -a(2)  )/(a(2));
        f4 = abs(           abs(r6)/(abs(r2)*abs(r3)*abs(r5)*abs(c1)*abs(c2))                -a(3)  )/(a(3));
        f5 = abs(  abs(c3)/abs(c1) - k );
        
        f6 = 0.04*f1+0.24*f2+0.24*f3+0.47*f4+0.01*f5;
        

    
    
    output=f6;
    end
end

function p = zeros3_ff(ftype,r,c,a,b)
p = @zeros3_ff;

   


    function output = zeros3_ff(input)


    
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        r7 = r(round(input(:,7)));
        r8 = r(round(input(:,8)));  
        r9 = r(round(input(:,9)));
        r10 = r(round(input(:,10))); 
        c1 = c(round(input(:,11)));
        c2 = c(round(input(:,12)));
        c3 = c(round(input(:,13)));
        if(b(2)<0.0001)
            b(2)=0.001;
        end
            

        %%%%%%%%%%%%%%%%%%%%% formula functions
        
        k=1000;
        B1 = abs(c3);
        B2 = abs(r6)/(abs(r3)*abs(r8))   +   1/abs(r4);
        B3 = abs(r6)/(abs(r3)*abs(r5)*abs(r7)*abs(c2));
        A1 = (   abs(c1)+abs(c3)   )*(   abs(r9)/(abs(r9)*abs(r10))   )+abs(c1);
        A2 = (   1+1/(abs(r1))+1/(abs(r4))    )*(    abs(r9)/(abs(r9)*abs(r10))    )+   1/abs(r1);
        A3 = abs(r6)/(abs(r2)*abs(r3)*abs(r5)*abs(c2));
        
        f1 = abs(  B2/B1    -b(2)  )/(b(2));
        f2 = abs(  B3/B1    -b(3)  )/(b(3));
        f3 = abs(  A2/A1    -a(2)  )/(a(2));
        f4 = abs(  A3/A1    -a(3)  )/(a(3));
        f5 = abs(  B1/A1  - k );
        
        f6 = 0.04*f1+0.24*f2+0.24*f3+0.47*f4;
        

    
    
    output=f6;
    end
end

function p = zeros_ff(ftype,r,rm,c,a,b)
p = @zeros_ff;

   


    function output = zeros_ff(input)


    
        r1 = r(round(input(:,1)));
        r2 = r(round(input(:,2)));
        r3 = r(round(input(:,3)));
        r4 = r(round(input(:,4)));
        r5 = r(round(input(:,5)));
        r6 = r(round(input(:,6)));
        r7 = r(round(input(:,7)));
        r8 = r(round(input(:,8)));
        r9 = r(round(input(:,9)));
        r10 = r(round(input(:,10)));
        r11 = r(round(input(:,11)));
        r12 = r(round(input(:,12)));
        r13 = r(round(input(:,13)));
        c1 = c(round(input(:,14)));
        c2 = c(round(input(:,15)));
        
        
          

        %%%%%%%%%%%%%%%%%%%%% formula functions
        k=1;
        
        A = abs(   abs(r4)*(abs(r5)+abs(r6))   /    (abs(r5)*(abs(r3)+abs(r4)) )  );     
        B = abs(   A*abs(r3)/(abs(r1)*abs(r4)*abs(c1))    );
        D = abs(   abs(r6)/(abs(r5)*abs(r1)*abs(r2)*abs(c1)*abs(c2))    );
        
        k1 = abs(r13)*abs(r8)/(abs(r9)*abs(r7));
        k2 = abs(r13)*abs(r11)/(abs(r12)*abs(r10));
        
        f1 = abs(      A*k1    -b(1)  )/(b(1));
        f2 = abs(   A*k2/(abs(r1)*abs(r2)*abs(c1)*abs(c2))   -  b(3)  )/(b(3));
        
        f3 = abs(B-a(2))/a(2);
        f4 = abs(D-a(3))/a(3);
        
        
       
        f5 = 0.1*f3+0.1*f1+0.1*f2+0.8*f4;
        

    
    
    output=f5;
    end
end



end

