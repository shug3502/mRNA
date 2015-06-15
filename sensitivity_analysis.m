function sensitivity_analysis(k)
p_indices = [1,4,6];
par_params = [1.16,0.8,0.11,0.42,0.84,0.58,0.5,0];

M=5;
for i=1:length(p_indices)
    sens_params(i,1:M) = [10.^linspace(-2,2,M)];
    if p_indices(i)==4
        sens_params(i,1:M) = 10.^linspace(-4,0,M);
    end
end
q_size = k*3 + (1-k)*51;
%q_store = zeros(,10,q_size);

my_index = 0;
while my_index<length(p_indices)
    my_index=my_index+1;
    for i1=1:M
        for i2=1:M
            for i3=1:M
                [i1,i2,i3]
                par_params(p_indices) = [sens_params(1,i1),sens_params(2,i2),sens_params(3,i3)];
                
                [q_estimate] = summary_statistic_calculator(par_params,20,1,k);
                
                q_store((M-1)^2*i1+(M-1)*i2+i3,:) = q_estimate;
            end
        end
    end
end
q_store
end