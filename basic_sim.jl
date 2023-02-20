#import Pkg
#Pkg.add("Plots")

# initializing the global conditions 
num_freqs = Int64(200)
N = Int64(1e3)
two_dim_map = zeros(ComplexF64, num_freqs, N);
smarter_result_matrix = zeros(Complex{Float64},2,N);
end_V_I = [1 1/50]'

function initializing_omega(omega) 
    # this creates the ABCD matrix for a given omega
    #

    Z_0 = 50;
    C_0 = 5e-3;
    L_0 = (Z_0^2)*C_0
    d = 1;
    
    L = 100
    C = 0.1
    
    series_matrix = [1 im*L_0*omega - im/(omega*C) ; 0 1];
    parallel_matrix = [1 0 ; im*C_0*omega - im/(omega*L)  1];
    
    return series_matrix*parallel_matrix
end;


# simulating the voltage across distance for 200 frequencies

freq_space = LinRange(0.75, 1.6, num_freqs);
for j = 1:200
    cell_mat = initializing_omega( freq_space[j] );
    
    for i = 1:N
        smarter_result_matrix[:, N+1-i] = (cell_mat^i)*end_V_I;
    end
    
    two_dim_map[Int64(j),:] = smarter_result_matrix[1,:];
end


# filtering extreme amplitudes

for m = 1:size(two_dim_map,1),  n = 1:size(two_dim_map, 2)
    if real(two_dim_map[m,n]) > 1
        two_dim_map[m,n] = 2/maximum(abs.(two_dim_map[m,:]))
    end
end


# potting the results 

using Plots
heatmap(1:N, 1:num_freqs, real(two_dim_map), c=:thermal)

