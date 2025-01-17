%% #############################################################################################
% Adaptive-Crossover-based-Marine-Predators-Algorithm-AC-MPA-
% This is the original code of Adaptive Crossover-based Marine Predators Algorithm (AC-MPA). 
% The code was developed by Assist. Prof. Dr. Shaymah Akram Yasear.
% 
% The mathematical details of AC-MPA can be found in: Shaymah Akram Yasear, 
% Adaptive crossover-based marine predators algorithm for global optimization problems, Journal of 
% Computational Design and Engineering, Volume 11, Issue 4, August 2024, Pages 124–150, 
% https://doi.org/10.1093/jcde/qwae060
% 
% Contact details of author: shayma.akram.yasear@gmail.com
% https://www.researchgate.net/profile/Shaymah-Yasear
% https://www.scopus.com/authid/detail.uri?authorId=57194554063
% https://www.webofscience.com/wos/author/record/AAN-5105-2021
%% #############################################################################################

function [Top_predator_fit, Top_predator_pos , Convergence_curve] = ACMPA(options)
    % Extract options
    fno = options.CEC_fun_no;
    objf = options.objf;
    dim = options.nVar;
    lb = options.lb;
    ub = options.ub;
    N = options.popSize;
    Max_iter = options.MaxFES;

    % Preallocate arrays
    Top_predator_pos = zeros(1, dim);
    Top_predator_fit = inf;
    fitness = inf(N, 1);
    MShcFitness = zeros(N, 1);
    MSvcFitness = zeros(N, 1);
    Convergence_curve = zeros(1, Max_iter);

    % Initialize Prey using adaptiveSamplingMaximinDistance function
    Prey = adaptiveSamplingMaximinDistance(lb, ub, dim, N, Max_iter);

    % Initialize variables
    Xmin=repmat(ones(1,dim).*lb,N,1);
    Xmax=repmat(ones(1,dim).*ub,N,1);
    fit_old = fitness;
    Prey_old = Prey;
    stepsize = zeros(N,dim);
    P =0.5;
    Pv = rand;
    FADs=0.2;
    % Main loop
    for Iter = 0:Max_iter-1
        %------------------- Detecting top predator -----------------
        % Apply bounds to Prey matrix
        Prey = max(min(Prey, ub), lb);
        
        for i = 1:N
            fitness(i) = objf(Prey(i, :)', fno);
            if fitness(i) < Top_predator_fit
                Top_predator_fit = fitness(i);
                Top_predator_pos = Prey(i, :);
            end
        end
        %------------------- Marine Memory saving ------------------- 
         if Iter==0
           fit_old=fitness;    Prey_old=Prey;
         end
         
         Inx=(fit_old<fitness);
         Indx=repmat(Inx,1,dim);
         Prey=Indx.*Prey_old+~Indx.*Prey;
         fitness=Inx.*fit_old+~Inx.*fitness;
         fit_old=fitness;    Prey_old=Prey;
        %------------------------------------------------------------  
        % Update Elite
        Elite = repmat(Top_predator_pos, N, 1);

        % Update CF
        CF = (1 - Iter / Max_iter) ^ (2 * Iter / Max_iter);

        % Levy and Brownian random numbers
        RL = 0.05 * levy(N, dim, 1.5);
        RB = randn(N, dim);

        % Update Prey positions
        for i = 1:N
            for j = 1:dim
                R=rand();
                  %------------------ Phase 1 (Eq. 4) ------------------- 
                if Iter<Max_iter/3 
                  stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*Prey(i,j));                    
                  Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j); 

                  %--------------- Phase 2 (Eqs. 5)----------------
                elseif Iter>Max_iter/3 && Iter<2*Max_iter/3 

                 if i>size(Prey,1)/2
                    stepsize(i,j)=RB(i,j)*(RB(i,j)*Elite(i,j)-Prey(i,j));
                    Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j); 
                 else
                    stepsize(i,j)=RL(i,j)*(Elite(i,j)-RL(i,j)*Prey(i,j));                     
                    Prey(i,j)=Prey(i,j)+P*R*stepsize(i,j);  
                 end  

                 %----------------- Phase 3 (Eq. 7)-------------------
                else 
                   stepsize(i,j)=RL(i,j)*(RL(i,j)*Elite(i,j)-Prey(i,j)); 
                   Prey(i,j)=Elite(i,j)+P*CF*stepsize(i,j);  
                end
            end
        end
        
      %------------------ Detecting top predator ------------------   
      Prey = max(min(Prey, ub), lb);
      for i=1:size(Prey,1)
         fitness(i,1) = objf(Prey(i,:)',fno);  
         if fitness(i,1)<Top_predator_fit 
            Top_predator_fit=fitness(i,1);
            Top_predator_pos=Prey(i,:);
         end     
      end

        %---------------------- Marine Memory saving ----------------

     if Iter==0
        fit_old=fitness;    Prey_old=Prey;
     end

        Inx=(fit_old<fitness);
        Indx=repmat(Inx,1,dim);
        Prey=Indx.*Prey_old+~Indx.*Prey;
        fitness=Inx.*fit_old+~Inx.*fitness;

        fit_old=fitness;    Prey_old=Prey;

        %---------- Eddy formation and FADs. effect (Eq. 8) -----------          
        if rand() < FADs
            U = rand(N, dim) < FADs;
            Prey = Prey + CF * ((Xmin + rand(N, dim) .* (Xmax - Xmin)) .* U);
        else
            r = rand();
            Rs = size(Prey, 1);
            stepsize = (FADs * (1 - r) + r) * (Prey(randperm(Rs), :) - Prey(randperm(Rs), :));
            Prey = Prey + stepsize;
        end

        %------------ HORIZONTAL CROSSOVER -------------------------------
        MShc = Prey;
        for i = 1:N
            L = randperm(N, 1);
            % Create array of possible values for M (all numbers except L)
            possible_M = setdiff(1:N, L);
            % Randomly select one value from possible_M
            M = possible_M(randi(N-1));
            for j = 1:dim
                c1 = -1 + 2 * rand;
                r1 = rand;
                MShc(L, j) = r1 * Prey(L, j) + (1 - r1) * Prey(M, j) + c1 * (Prey(L, j) - Prey(M, j));
                c2 = -1 + 2 * rand;
                r2 = rand;
                MShc(M, j) = r2 * Prey(M, j) + (1 - r2) * Prey(L, j) + c2 * (Prey(M, j) - Prey(L, j));
            end
        end
        
        MShc = max(min(MShc, ub), lb);
        
        % Update fitness and Prey positions based on new positions
        for i = 1:N
            MShcFitness(i) = objf(MShc(i, :)', fno);
            if MShcFitness(i) < fitness(i)
                fitness(i) = MShcFitness(i);
                Prey(i, :) = MShc(i, :);
            end
        end
        

        %------------------- VERTICAL CROSSOVER ---------------------------
        MSvc = Prey;
        B = randperm(dim);
        for j = 1:floor(dim / 2)
            r3 = rand;
            if r3 < Pv
                d1 = B(2 * j - 1);
                d2 = B(2 * j);
                for i = 1:N
                    r4 = rand;
                    MSvc(i, d1) = r4 * MSvc(i, d1) + (1 - r4) * MSvc(i, d2);
                end
            end
        end
        MSvc(i, :) = max(min(MSvc(i, :), ub), lb);

        % Update fitness and Prey positions based on new positions
        for i = 1:N
            MSvcFitness(i) = objf(MSvc(i, :)', fno);
            if MSvcFitness(i) < fitness(i)
                fitness(i) = MSvcFitness(i);
                Prey(i, :) = MSvc(i, :);
            end
        end
     % Update convergence curve
      Convergence_curve(Iter+1) = Top_predator_fit;
    end

end
