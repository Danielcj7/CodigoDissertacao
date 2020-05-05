clc; clear all; close all;

H = tf([1000.57653929],[1 0.8934 0]); %SISTEMA IDENTIFICADO

format long;  
t=3; %tempo de simulaçao do sinal.
Ts = 0:0.01:t; %discretização
ref=gensig('sin',2,3,0.01); %geração de sinal de referência
ref=ref+0.01*randn(size(ref));    %Ruído 
A=1;   %aleatoriedade dos indivíduos
ordem = 2; %ordem do sistema
popini = 50;  %população inicial
sigma = 50;   %tamanho do kernel
precmax = 20;   % Faixa de entropia.
precmin = -10;               
normmax = 8;  %Norma H2
normmin = 1;
pm = 0.25;  %porcentagem de mutação
maxepoca = 100; %iterações

Result(200,1)=0; %inicialização de matrizes
Resultfuncoes(200,:)=tf(1,1);
G(popini,:)=tf(1,1);
Hmfg(popini,:)=tf(1,1);
xx=1;
Exx=zeros(1,2);

K=2;   %quantidade de resets

for cont=1:K;
tic
randn('state',sum(A*randn*clock)); %reseta o gerador aleatório
indiv = A*rand(popini,ordem*2+2); %gera individuos
epoca=0;
E=zeros(popini,1);

    % Funcao fitness
    for i=1:size(indiv,1)
        [E(i,:),G(i,:),Hmfg(i,:)] = fitness_n(indiv,H,Ts,ref,sigma,i);
    end
    i=1; j=1;
    while i ~= (length(E(:,1))+1)  %substituindo os NaN e inf
        B=isnan(E(i,:));  C=isinf(E(i,:));
        if B==1 || C==1
            indiv(i,:)=A*rand(1,ordem*2+2);    %substitui por outro
            [E(i,:),G(i,:),Hmfg(i,:)] = fitness_n(indiv,H,Ts,ref,sigma,i);
        else i=i+1;
        end
    end
    Ex(:,:)=50;   %coloca-se um valor fora das precisoes para que se passe nos "if's"
    for i=1:size(indiv,1)
        if ( (E(i,:)~=0) && (E(i,:) < precmax) && (E(i,:) > precmin) )
            Ex(j,:)= E(i,:);   %joga em novas variaveis pra ver a minima depois
            Ex(j,2)= i;
            j=j+1;
        end    
    end
    j=1; Exx(:,:)=50;
    for i=1:size(Ex,1)
        [n]=(Ex(i,2));
        if (((max(real(roots(G(n,:).den{1})))) < -0.002) && (norm(Hmfg(n,:)) < normmax) && (norm(Hmfg(n,:)) > normmin)) %verifica polos parta real negativa
            Exx(j,:)= Ex(i,:);   
            j=j+1;
        end    
    end
    if (  (min(Exx(:,1)) < precmax) && ((min(Exx(:,1))) > precmin) )
        [m,h] = min(Exx(:,1))
        [k]=(Exx(h,2));
        G(k,:)
        epoca=epoca+1; disp(['Epoca= ' num2str(epoca)]); disp(' ');
        Result(xx,:)= min(Exx(:,1));
        Resultfuncoes(xx,:)=G(k,:);
        salvar = ['Results sig ',num2str(sigma),'.mat'];
        save (salvar, 'Result', 'Resultfuncoes', 'H', 'Ts', 'A', 'popini', 'sigma', 'precmax', 'pm', 'maxepoca', 'cont', 'ref', 'epoca')
        xx=xx+1;
        precmax=min(Exx(:,1));
    end
    
while epoca<maxepoca
    % Seleção
    disp('== Seleção ==')    
    fit=1./E;    %Para inverter, o menor número ser o melhor    
    aux = cumsum(fit);   %Soma cumulativa do numero + anteriores
    prob = aux/sum(fit);  %Estabelece o intervalo de probabilidade.
    for i=1:size(indiv,1)   %Roleta no mesmo número de indivíduos.
        x = rand(1,1);      %Valor aleatório de 0 a 1.
        P = find(prob>x);  %Procura qual indivíduo possuí o valor sorteado.
        nindiv(i,:)=indiv(P(1,1),:);  %individuos selecionados.
    end
    indiv=nindiv;    %Substitui a população pela nova.
   
    % Crossover
    disp('== Crossover ==')
    for i=1:(size(indiv,1)/2)   %Seleção de cruzamento
        x1=randi([1 size(indiv,1)],1,1);  %escolha aleatória de pais
        pai1=indiv(x1,:);
        x2=randi([1 size(indiv,1)],1,1); 
        pai2=indiv(x2,:);
        alpha=rand;      %gerando filhos
        filho1=alpha*pai1+(1-alpha)*pai2;  %equação utilizada para cruzar os genes
        filho2=(1-alpha)*pai1+alpha*pai2;
        indivaux=nindiv;
        nindiv=[indivaux; filho1; filho2];  %adicionando filhos criados a população
    end 
    indiv=nindiv;
    for i=1:size(indiv,1)  %testando novamente os individuos
       [E(i,:),G(i,:),Hmfg(i,:)] = fitness_n(indiv,H,Ts,ref,sigma,i);
    end
       %Deixando apenas os melhores individuos
    for i = 1:popini    %exclui os com entropia de maior valor
        [m,n]=max(E);      
            G(n,:) = [];
            nindiv(n,:) = [];
            E(n,:) = [];
            Hmfg(n,:)=[];
    end
    indiv=nindiv;
    i=1; j=1;
    while i ~= (length(E(:,1))+1)    %substituindo os NaN e inf
        B=isnan(E(i,:));  C=isinf(E(i,:));
        if B==1 || C==1
            indiv(i,:)=A*rand(1,ordem*2+2);    %substitui por outro
         [E(i,:),G(i,:),Hmfg(i,:)] = fitness_n(indiv,H,Ts,ref,sigma,i);
        else i=i+1;
        end
    end
    Ex(:,:)=50;
    for i=1:size(indiv,1)
    if ( (E(i,:)~=0) && (E(i,:) < precmax) && (E(i,:) > precmin) )
        Ex(j,:)= E(i,:);   %joga em novas variaveis pra ver a minima depois
        Ex(j,2)= i;
        j=j+1;
    end    
    end
    j=1; Exx(:,:)=50;
    for i=1:size(Ex,1)
        [n]=(Ex(i,2));
        if (((max(real(roots(G(n,:).den{1})))) < -0.002) && (norm(Hmfg(n,:)) < normmax) && (norm(Hmfg(n,:)) > normmin)) %verifica polos parta real negativa
            Exx(j,:)= Ex(i,:);   
            j=j+1;
        end    
    end
        if (  (min(Exx(:,1)) < precmax) && ((min(Exx(:,1))) > precmin) )
        [m,h] = min(Exx(:,1))
        [k]=(Exx(h,2));
        G(k,:)
        epoca=epoca+1; disp(['Epoca= ' num2str(epoca)]); disp(' ');
        Result(xx,:)=min(Exx(:,1));
        Resultfuncoes(xx,:)=G(k,:);
        salvar = ['Results sig ',num2str(sigma),'.mat'];
        save (salvar, 'Result', 'Resultfuncoes', 'H', 'Ts', 'A', 'popini', 'sigma', 'precmax', 'pm', 'maxepoca', 'cont', 'ref', 'epoca')
        precmax=min(Exx(:,1));
        xx=xx+1;
    end
    
    % Mutacao
    disp('== Mutação ==')
    i=1;
    for i = 1:size(indiv,1)  %percorre todos os elementos de todos os individuos
        for j = 1:length(indiv(i,:))
           if pm>rand    %sofre alteração se o numero sorteado for menor que a porcentagem escolhida
            indiv(i,j) = indiv(i,j)*(1 + A*randn(1,1));
           end
        end
          [E(i,:),G(i,:),Hmfg(i,:)] = fitness_n(indiv,H,Ts,ref,sigma,i);
    end
    i=1; j=1;
    while i ~= (length(E(:,1))+1)  %substituindo os NaN e inf
        B=isnan(E(i,:));  C=isinf(E(i,:));
        if B==1 || C==1
            indiv(i,:)=A*rand(1,ordem*2+2); %substitui por outro
            [E(i,:),G(i,:),Hmfg(i,:)] = fitness_n(indiv,H,Ts,ref,sigma,i);
        else i=i+1;
        end
    end
    Ex(:,:)=50;
    for i=1:size(indiv,1)
    if ( (E(i,:)~=0) && (E(i,:) < precmax) && (E(i,:) > precmin) )
        Ex(j,:)= E(i,:);   %joga em novas variaveis pra ver a min depois
        Ex(j,2)= i;
        j=j+1;
    end    
    end
    j=1; Exx(:,:)=50;
    for i=1:size(Ex,1)
        [n]=(Ex(i,2));
        if (((max(real(roots(G(n,:).den{1})))) < -0.002) && (norm(Hmfg(n,:)) < normmax) && (norm(Hmfg(n,:)) > normmin)) %verifica polos parta real negativa
            Exx(j,:)= Ex(i,:);   
            j=j+1;
        end    
    end
        if (  (min(Exx(:,1)) < precmax) && ((min(Exx(:,1))) > precmin) )
        [m,h] = min(Exx(:,1))
        [k]=(Exx(h,2));
        G(k,:)
        epoca=epoca+1; disp(['Epoca= ' num2str(epoca)]); disp(' ');
        Result(xx,:)=min(Exx(:,1));
        Resultfuncoes(xx,:)=G(k,:);
        salvar = ['Results sig ',num2str(sigma),'.mat'];
        save (salvar, 'Result', 'Resultfuncoes', 'H', 'Ts', 'A', 'popini', 'sigma', 'precmax', 'pm', 'maxepoca', 'cont', 'ref', 'epoca')
        precmax=min(Exx(:,1));
        xx=xx+1;
        end
    epoca=epoca+1; disp(['Epoca= ' num2str(epoca)]); disp(' ');
    disp(['Minimo atual=' num2str(min(Exx))]);
    Exx=zeros(1,2);
    disp(' ');
end
toc
cont
end

    %% -------------- Função Fitness -----------------
function [E, G, Hmfg] = fitness_n(indiv,H,Ts,ref,sigma,i)

        G=tf([indiv(i,1:length(indiv(i,:))/2)],[indiv(i,((length(indiv(i,:))/2)+1):(length(indiv(i,:))))]);
        Hmfg = feedback(G*H,1);
        [y]=lsim(Hmfg,ref,Ts);
        e = ref - y;                
        E = renyi2(abs(e),sigma);
        
       
    %% --------------- Função MEE ------------------        
% Renyi's quadratic entropy estimation using Parzen windowing
% and Gaussian kernels
%
% H2 = renyi2(x,s)
%
%   H2 - Renyi's quadratic entropy (bits)
%   x - signal
%   s - Gaussian kernel size
%
% From 
% Weifeng Liu, Puskal P. Pokharel, and Jose C. Principe. 
% Correntropy: Properties and Applications in Non-Gaussian Signal Processing. 
% IEEE Transactions on Signal Processing, v. 55, n. 11, nov. 2007, 5286-5298
%
% PRS, 07/jun/2017
%
function H2 = renyi2(x,s)

if nargin ~= 2
    error('Wrong number of parameters')
end
x = x(:);
N = length(x);
X = repmat(x',N,1)-repmat(x,1,N);
IP = kernel_s2(X,sqrt(2)*s);
IP = sum(IP(:))/(N^2);
H2 = -log2(IP);
clear x N X IP 

    %% -------------- Função Kernel --------------------- 
% Gaussian Kernel estimation
%
% PRS, fev/08/2016
% modified in 07/jun/2017 -> x can be a matrix
%
function k = kernel_s2(x,s)

k = (1/(sqrt(2*pi)*s)) * exp((-x.^2)/(2*s^2));
