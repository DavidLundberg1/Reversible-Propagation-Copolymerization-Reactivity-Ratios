function[conversion_vector] = BinaryReversibleStochasticSimulation(rates,chain_length)
%Stochastic simulation script for binary copolymerization with reversible
%propagation
%Inputs: rates, a vector containing [f1_0, r1, r2, k2, k6, k8] 
%f1_0 = initial monomer 1 mole fraction
%M0 = initial monomer concentration (M)
%k1, k2, k3, k4, k5, k6, k7, k8 are rate constants as defined in Scheme 1
%of the manuscript. 

%Output is 'conversion_vector' which contains total monomer conversion, and
%individual monomer conversion values at evert 100th simulation step. 

k1 = 1; %Monomer 1 homo-propagation
k2 = rates(4); %Monomer 1 homo-depropagation *zero for now
k3 = 1/rates(2); %Monomer 1 cross-propagation 
k4 = rates(7); %Monomer 2 cross-depropagation 
k5 = 1/rates(3); %Monomer 2 cross-propagation
k6 = rates(5); %Monomer 1 cross-depropagation *zero for now
k7 = 1; %Monomer 2 homo-propagation
k8 = rates(6); %monomer 2 homo-depropagation 

M0 = 1; %M, Total initial monomer concentration (defined here)
M1_0 = rates(1)*M0; %Monomer 1 initial concentration
M2_0 = M0-M1_0; %Monomer 2 initial concentration

f1_0 = rates(1);
f2_0 = 1-f1_0;

%Define size of simulation (how many monomers to include)
num_mon = 100*chain_length;
num_chains = num_mon/chain_length;

prop_prob_tally = zeros(8,0); %Track the number of times each reaction occurs
frac_terminal_counter=zeros(2,0);
frac_homo_counter=zeros(2,0);

%Store chains in polymer_array and initiate chains randomly
polymer_array = zeros(num_chains,chain_length*3);
polymer_array(1,1)=1;
polymer_array(2,1)=2;
for i = 3:num_chains
    X = rand;
    if X <= f1_0
    polymer_array(i,1)= 1;
    else
        polymer_array(i,1)=2;
    end     
end

%Tally amount of monomers in monomer_amounts
[~,~,ix]=unique(polymer_array);
content = accumarray(ix,1);
monomer_amounts = zeros(2,1);
monomer_amounts(1) = f1_0*num_mon - content(2);
monomer_amounts(2) = f2_0*num_mon - content(3);

%Define terminal_index of each chain (which index is propagating) and its
%identity in terminal_monomer_list
terminal_index = ones(num_chains,1);
terminal_monomer_list = polymer_array(:,1);

%Identity of penultimate monomer
penultimate_monomer_list = terminal_monomer_list;

%Indices of chains with terminal monomer 1 or 2 units
terminal_one_index = find(terminal_monomer_list<2);
terminal_two_index = find(terminal_monomer_list>1);

%Fraction of chains with a given terminal unit
frac_terminal_1 = length(terminal_one_index)/(num_mon/chain_length);
frac_terminal_2 = 1-frac_terminal_1;

%Track identity of terminal dyad for each chain (1 or 0)
homo_dyad_list1 = zeros(1,length(terminal_monomer_list));
hetero_dyad_list1 = zeros(1,length(terminal_monomer_list));
homo_dyad_list2 = zeros(1,length(terminal_monomer_list));
hetero_dyad_list2 = zeros(1,length(terminal_monomer_list));

%Indeces of each type of terminal dyad
terminal_homo_index1 = find(homo_dyad_list1>0);
terminal_hetero_index1 = find(hetero_dyad_list1>0);
terminal_homo_index2 = find(homo_dyad_list2>0);
terminal_hetero_index2 = find(hetero_dyad_list2>0);

%Fraction of chains with a given type of terminal dyad
homo_dyad_frac2 = length(terminal_homo_index2)/length(terminal_two_index);%Fraction of chain ends 2 having a homo dyad
homo_dyad_frac1 = length(terminal_homo_index1)/length(terminal_one_index);%Fraction of chain ends 1 having a homo dyad

%Define initial fractions and chain/monomer amounts
active_chains = 1:num_chains;
mon_left = sum(monomer_amounts);
f1=monomer_amounts(1)*M0/num_mon; 
f2=1-f1;
mon_converted = num_mon - mon_left;

%Define rates in terms of global variables 
k11 = k1;
k12 = k3;
k22 = k7;
k21 = k5;

%Reaction rate (constants)
rate_11 = k11*f1;
rate_12 = k12*f2;
rate_21 = k21*f1;
rate_22 = k22*f2;
rate_1homodeprop = k2;
rate_1heterodeprop = k6;
rate_2homodeprop = k8;
rate_2heterodeprop = k4;

%Probability of each possible reaction
prop_prob_vec = [
rate_11*frac_terminal_1;
rate_12*frac_terminal_1; 
rate_22*frac_terminal_2;
rate_21*frac_terminal_2;
rate_1homodeprop*frac_terminal_1*homo_dyad_frac1;
rate_1heterodeprop*frac_terminal_1*(1-homo_dyad_frac1);
rate_2homodeprop*frac_terminal_2*homo_dyad_frac2;
rate_2heterodeprop*frac_terminal_2*(1-homo_dyad_frac2)];
norm_prop_prob_vec = prop_prob_vec/sum(prop_prob_vec);

%Define vector to store conversion versus time data in
conversion_vector = zeros(3,num_mon/100);

conversion_vector(1,1) = 1 - mon_left/num_mon; %Fraction of all monomers converted
conversion_vector(2,1) = 1 - (monomer_amounts(1)/(num_mon*f1_0));%Fraction of monomer type 1 converted
conversion_vector(3,1) = 1 - (monomer_amounts(2)/(num_mon*f2_0));%Fraction of monomer type 2 converted

%Track each loop in the simulation
loop_count = 0;

%Tally each reaction type
event_list = [0,0,0,0,0,0,0,0];

equilibrium_condition_count = 0; 

%Define vectors to count terminal run lengths
terminal_runs_lengths = zeros(2,chain_length*3,num_mon/1000);

%Propagation
while sum(monomer_amounts) > 0 
    loop_count = loop_count+1;
   
    %Randomly pick one of the reactions to occur by rxn_number
    x = rand;
    prop_cum_sum = cumsum(norm_prop_prob_vec)-x; 
    rxn_number = find(prop_cum_sum>0,1); 
    
    if rxn_number ==1 %M1 Homopropagation
       %Pick random chain
       chain = terminal_one_index(randperm(length(terminal_one_index),1)); 
       
       %Add unit into chain
       polymer_array(chain,terminal_index(chain)+1)=1; 
       
       %Subtract propagated monomer from pool
       monomer_amounts(1) = monomer_amounts(1) -1;
       
       %Update chain terminal index
       terminal_index(chain)=terminal_index(chain)+1;
       
       %Subtract from mon_converted
       mon_converted=mon_converted+1;
       
       %Tally event_list vector
       event_list(1)=event_list(1)+1;
       
       %Update homodyad list and remove from heterodyad list if it were. 
       homo_dyad_list1(chain)=1;
       hetero_dyad_list1(chain)=0;

    elseif rxn_number == 2 %M1 Crosspropagation
           chain = terminal_one_index(randperm(length(terminal_one_index),1)); 
           polymer_array(chain,terminal_index(chain)+1)=2;
           monomer_amounts(2) = monomer_amounts(2)-1;
          
           %Chainge terminal monomer identity
           terminal_monomer_list(chain)=2; 
           terminal_index(chain)=terminal_index(chain)+1; 
           mon_converted=mon_converted+1;
           
           %Update penultimate identity
           %Specify chain now as heterodyad
           hetero_dyad_list2(chain)=1;
           homo_dyad_list1(chain)=0;
           hetero_dyad_list1(chain)=0;
           
           %Tally event_list vector
           event_list(2)=event_list(2)+1;
   
    elseif rxn_number == 3 %M2 Homopropagation
        %Pick a random chain with M2 terminal index. 
            chain = terminal_two_index(randperm(length(terminal_two_index),1));
            polymer_array(chain,terminal_index(chain)+1)=2; 
            monomer_amounts(2) = monomer_amounts(2) - 1;
            terminal_index(chain)=terminal_index(chain)+1;
            mon_converted=mon_converted+1;
          
           %Specify as homodyad
           homo_dyad_list2(chain)=1;
           hetero_dyad_list2(chain)=0;
            event_list(3)=event_list(3)+1;
    
    elseif rxn_number == 4 %M2 Crosspropagation
            chain = terminal_two_index(randperm(length(terminal_two_index),1)); 
            polymer_array(chain,terminal_index(chain)+1)=1; 
            monomer_amounts(1) = monomer_amounts(1) - 1;
            terminal_index(chain)=terminal_index(chain)+1;
           mon_converted=mon_converted+1;
           terminal_monomer_list(chain)=1; 
           
           %Erase terminal dyad status since we irreversibly prop'd M1
           homo_dyad_list2(chain)=0;
           hetero_dyad_list2(chain)=0;
           hetero_dyad_list1(chain)=1;
           event_list(4)=event_list(4)+1;
   
    elseif rxn_number == 5; %M1 Homodyad depropagation
        if length(terminal_homo_index1)>1 %There are some to depropagate
            chain = terminal_homo_index1(randperm(length(terminal_homo_index1),1)); 
            if terminal_index(chain)>1
            polymer_array(chain,terminal_index(chain))= 0; 
            terminal_index(chain)=terminal_index(chain)-1;
            monomer_amounts(1)=monomer_amounts(1)+1;
            mon_converted = mon_converted-1;
            terminal_monomer_list(chain)=1;
            if terminal_index(chain) == 1
                    homo_dyad_list1(chain)=0;
                    hetero_dyad_list1(chain)=0;
                    homo_dyad_list2(chain)=0;
                    hetero_dyad_list2(chain)=0;
                else
            if  polymer_array(chain,terminal_index(chain)-1) == 1;
                homo_dyad_list1(chain)=1;
            else
                hetero_dyad_list1(chain)=1;
                homo_dyad_list1(chain)=0;
            end
            event_list(5)=event_list(5)+1;
            end
            end
        end

    elseif rxn_number == 6 %M1 Heterodyad depropagation
            if length(terminal_hetero_index1)>1
            chain = terminal_hetero_index1(randperm(length(terminal_hetero_index1),1)); 
               
            if terminal_index(chain)>1 %Don't deprop @ beginning
                polymer_array(chain,terminal_index(chain))= 0;
                 terminal_index(chain)=terminal_index(chain)-1;
                 monomer_amounts(1)=monomer_amounts(1)+1;
                    mon_converted = mon_converted-1;
                 terminal_monomer_list(chain)=2;
                if terminal_index(chain) == 1
                    homo_dyad_list1(chain)=0;
                    hetero_dyad_list1(chain)=0;
                    homo_dyad_list2(chain)=0;
                    hetero_dyad_list2(chain)=0;
                else
                if polymer_array(chain,terminal_index(chain)-1) == 2
                    homo_dyad_list2(chain)=1;
                    hetero_dyad_list1(chain)=0;
                else
                    hetero_dyad_list2(chain)=1;
                    hetero_dyad_list1(chain)=0;   
                end
            end
            %tally event list 
                event_list(6)=event_list(6)+1;
            end
            end
    elseif rxn_number == 7 %M2 Homodyad depropagation
        %Pick a random chain with a M2 homodyad at the end
        if length(terminal_homo_index2)>1
        chain = terminal_homo_index2(randperm(length(terminal_homo_index2),1)); 
        if terminal_index(chain)>1 %Don't deprop @ beginning
        polymer_array(chain,terminal_index(chain))= 0;
        %Change terminal index now. 
        terminal_index(chain)=terminal_index(chain)-1;
        monomer_amounts(2)=monomer_amounts(2)+1;
        mon_converted = mon_converted-1;
        terminal_monomer_list(chain)=2;
        
        %Update dyad status. Only need to change if it is a hetero dyad. 
        if terminal_index(chain) == 1
                    homo_dyad_list1(chain)=0;
                    hetero_dyad_list1(chain)=0;
                    homo_dyad_list2(chain)=0;
                    hetero_dyad_list2(chain)=0;
                else
        if  polymer_array(chain,terminal_index(chain)-1) ==1 
            hetero_dyad_list2(chain)=1;
            homo_dyad_list2(chain)=0;
        else
            homo_dyad_list2(chain)=1;
        end
        %Only tally if we actually do a step
        event_list(7)=event_list(7)+1;
        end
        end
        end 
        elseif rxn_number == 8 %M2 Heterodyad depropagation
             if length(terminal_hetero_index2)>1
            chain = terminal_hetero_index2(randperm(length(terminal_hetero_index2),1)); 
               
            if terminal_index(chain)>1 %Don't deprop @ beginning
                polymer_array(chain,terminal_index(chain))= 0;
                terminal_index(chain)=terminal_index(chain)-1; %Decrease terminal index
                monomer_amounts(2)=monomer_amounts(2)+1;
                mon_converted = mon_converted-1;
                terminal_monomer_list(chain)=1;
                %No longer a hetero dyad
                hetero_dyad_list2(chain)=0;
                %Do conditionl for previous unit change to no heterodyad lists when there
                %is some difference. Again check this out. 
                if terminal_index(chain) == 1
                    homo_dyad_list1(chain)=0;
                    hetero_dyad_list1(chain)=0;
                    homo_dyad_list2(chain)=0;
                    hetero_dyad_list2(chain)=0;
                else
                if polymer_array(chain,terminal_index(chain)-1) == 1
                    homo_dyad_list1(chain)=1;
                    hetero_dyad_list2(chain)=0;
                else
                    hetero_dyad_list1(chain)=1;
                    hetero_dyad_list2(chain)=0;
                end
            end
            %tally event list
                event_list(8)=event_list(8)+1;
            end
             end
    end

%Recalculate propagation probabilities
mon_left = monomer_amounts(1) + monomer_amounts(2);

%This means f1/f2 accounts for concentration. 
f1=monomer_amounts(1)*M0/num_mon; 
f2=monomer_amounts(2)*M0/num_mon;

%Reaction rate (constants)
rate_11 = k11*f1;
rate_12 = k12*f2;
rate_21 = k21*f1;
rate_22 = k22*f2;
rate_1homodeprop = k2;
rate_1heterodeprop = k6;
rate_2homodeprop = k8;
rate_2heterodeprop = k4;

%Need to update probability vector. 
terminal_one_index = find(terminal_monomer_list<2);
terminal_two_index = find(terminal_monomer_list>1);

terminal_one_index(isnan(terminal_one_index))=0;
terminal_two_index(isnan(terminal_two_index))=0;

%Fraction of chains with a given terminal unit
frac_terminal_1 = length(terminal_one_index)/(num_mon/chain_length);
frac_terminal_2 = 1-frac_terminal_1;


%Indeces of each type of terminal dyad
terminal_homo_index1 = find(homo_dyad_list1>0);
terminal_hetero_index1 = find(hetero_dyad_list1>0);
terminal_homo_index2 = find(homo_dyad_list2>0);
terminal_hetero_index2 = find(hetero_dyad_list2>0);

%Fraction of chains with a given type of terminal dyad
homo_dyad_frac2 = length(terminal_homo_index2)/length(terminal_two_index);%Fraction of chain ends 2 having a homo dyad
homo_dyad_frac1 = length(terminal_homo_index1)/length(terminal_one_index);%Fraction of chain ends 1 having a homo dyad

%Probability of each possible reaction
prop_prob_vec = [
rate_11*frac_terminal_1;
rate_12*frac_terminal_1; 
rate_22*frac_terminal_2;
rate_21*frac_terminal_2;
rate_1homodeprop*frac_terminal_1*homo_dyad_frac1;
rate_1heterodeprop*frac_terminal_1*(1-homo_dyad_frac1);
rate_2homodeprop*frac_terminal_2*homo_dyad_frac2;
rate_2heterodeprop*frac_terminal_2*(1-homo_dyad_frac2)];
prop_prob_vec(isnan(prop_prob_vec))=0;
norm_prop_prob_vec = prop_prob_vec/sum(prop_prob_vec);

%Define criteria for stopping simulation after equilibrium is reached
%(i.e., rate of propagation and depropagation become equal)
if (prop_prob_vec(5)+prop_prob_vec(6)) > (prop_prob_vec(1)+prop_prob_vec(4)) & (prop_prob_vec(7)+prop_prob_vec(8)) > (prop_prob_vec(3)+prop_prob_vec(2));
    equilibrium_condition_count = equilibrium_condition_count + 1;
end

if   equilibrium_condition_count > 0.001*num_mon %Define stopping condition for number of steps to stay at equilibrium
    break
end

if mod(loop_count,100) == 0 %Only write to 'converssion_vector' ever 100 loops to decrease the size of the output vector
conversion_vector(1,loop_count/100+1) = 1 - mon_left/num_mon; %Fraction of mon converted
conversion_vector(2,loop_count/100+1) = 1 - (monomer_amounts(1)/(num_mon*f1_0));
conversion_vector(3,loop_count/100+1) = 1 - (monomer_amounts(2)/(num_mon*f2_0));
prop_prob_tally(:,end+1) = norm_prop_prob_vec;
frac_terminal_counter(1,end+1) = frac_terminal_1;
frac_terminal_counter(2,end) = frac_terminal_2;
frac_homo_counter(1,end+1)=homo_dyad_frac2;
frac_homo_counter(2,end)=1-homo_dyad_frac2;
end

end

end