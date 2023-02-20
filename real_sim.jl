#import Pkg
#Pkg.add("Plots")

# initializing the global conditions 
num_freqs = Int64(100)
N = Int64(10000)
two_dim_map = zeros(ComplexF64, num_freqs, N);
smarter_result_matrix = zeros(Complex{Float64},2,N);
end_V_I = [1 1/50]';


function initializing_omega(omega) 
    # this creates the ABCD matrix for a given omega
    #

    C_0 = 0.17e-12;
    L_0 = 0.238e-9; 
    d = 9;
    Z_0 = (L_0/C_0)^-0.5;

    L = 2.8e-9
    C = 2.7e-12

    
    
    series_matrix = [1 im*L_0*d*omega - im/(omega*C) ; 0 1];
    parallel_matrix = [1 0 ; im*C_0*d*omega - im/(omega*L)  1];

    cell_matrix = series_matrix*parallel_matrix
    dx_matrix = cell_matrix^(1/(N/20))
    
    return dx_matrix
end;


# calculating a dx matrix



# simulating the voltage across distance for 100 frequencies
f_min = Int64(1e9); f_max = Int64(8e9);
freq_space = LinRange(2*pi*f_min, 2*pi*f_max, num_freqs);

for j = 1:num_freqs
    dx_matrix = initializing_omega( freq_space[j] );
    
    smarter_result_matrix[:,end] = end_V_I;
    for i = 1:N-1
        smarter_result_matrix[:, N-i] = (dx_matrix)*smarter_result_matrix[:, N-i+1];
    end
    
    two_dim_map[Int64(j),:] = smarter_result_matrix[1,:];
end


# filtering extreme amplitudes

for m = 1:size(two_dim_map,1)
    for  n = 1:size(two_dim_map, 2)
        if real(two_dim_map[m,n]) > 5
            max_val = maximum(real.(two_dim_map[m,:]))
            two_dim_map[m,:] = two_dim_map[m,:] * (5 / (max_val))
        end
    end
end


# potting the results 

using Plots
f_min_MHz = Int64(f_min*1e-9); f_max_MHz = Int64(f_max*1e-9);
jumps = Int64((f_max - f_min)*1e-9) + 1;
heatmap(1:num_freqs, 1:N, real.(two_dim_map'), c=:viridis, 
        xticks = (  1:num_freqs/(jumps -1):num_freqs, [i for i in LinRange(f_min_MHz, f_max_MHz, jumps)] ),
        xlabel = "frequency [MHz]",
        yticks = (0:N/10:N, [j for j in LinRange(0,9,10)]),
        ylabel = "distance [cm]",
        title = "Simualtion 2.8 [nh] 2.7 [pf]" 
        
        )



