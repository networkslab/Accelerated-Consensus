%% Instructions
%
% 1) Run scripts "generate_opt_matrices_RG.m" and
% "generate_opt_matrices_path.m" for every desired value of network size, 
% which is set in variable "Nvec". These scripts create the
% .mat files in the current folder. The files contain the FDLA optimal
% matrices and polynomial weights necessary to perform comparison with the
% proposed algorithm. Move all these files into folder named "RGG_matrices".
%
% 2) Run scripts "main_consensus_circ_lin_saved_mtx.m" and 
% "main_consensus_RG_lin_saved_mtx.m" for desired network sizes.
% These scripts simulate the proposed algorithm and the benchmark
% algorithms from the literature. "main_consensus_circ_lin_saved_mtx.m"
% implemets the simulation on the path graph with "Slope" initialization. 
% "main_consensus_circ_lin_saved_mtx.m" implemets the simulation on the random
% geometric graph with "Slope" initialization. Changing "Slope" to "Spike"
% in "main_consensus_circ_lin_saved_mtx.m" results in the Spike
% initialization. The scripts save the simulation output in the current
% directory. Move the output to the directory named "sim_results_2".
%
% 3) Run the plotting scripts. Their name starts with "plot_results_".
% These scripts generate all the figures used in the paper based on the
% data collected during steps 1 and 2. The figures are saved in the .eps
% format in the folder named "fig3".


%% File list
% create_E.m                    Create the connectivity matrix for the 2D random 
%                               geometric graph
% create_E_chain                Create the connectivity matrix for the path of N nodes
% create_MH                     Create the Metropolis-Hastings weight matrix based 
%                               on the adjacency matrix
% create_OPT                    Create the optimal symmetric weight matrix based on
%                               the adjacency matrix. The optimal weight matrix is 
%                               the solution of the symmetric FDLA problem from L. 
%                               Xiao and S. Boyd, “Fast linear iterations for distributed 
%                               averaging.” Sys. and Control Letters, vol. 53, no. 1, 
%                               pp. 65–78, Sep. 2004.
% create_POLY                   Create the optimal weights for the polynomial filter
% get_hermite                   Create the sub-optimal weights for the polynomial 
%                               filter based on the Hermite interpolating polynomial
% get_sund_haj_time             Estimate the time that the algorithm proposed by 
%                               Sundaram and Hajicostis requires to achieve perfect 
%                               convergence.
% get_alpha                     Calculate the optimal value of alpha, according to
%                               Theorem 1
% initialize                    Initialize the state vector with a realization of 
%                               one of possible field types
% do_consensus                  Perform consensus iterations, output the MSE curve 
%                               and the final value
% do_consensus_acc              Estimate the eigenvalue of the foundational weight 
%                               matrix (if instructed) calculate the optimal value 
%                               of the mixing parameter, perform accelerated consensus 
%                               iterations, output the MSE curve and the history of
%                               consensus iterations
% do_consensus_acc_circ         For the path graph, estimate the eigenvalue of the 
%                               foundational weight matrix (if instructed) calculate the optimal value 
%                               of the mixing parameter, perform accelerated consensus 
%                               iterations, output the MSE curve and the history of
%                               consensus iterations
% do_consensus_poly             Perform accelerated consensus iterations for the 
%                               polynomial filter, output the MSE curve and the 
%                               history of consensus iterations
% do_consensus_poly_circ        For the path graph, perform accelerated consensus 
%                               iterations for the polynomial filter, output the 
%                               MSE curve and the history of consensus iterations
% estimate_l2_pwr_3             Estimate the second (by the modulus) largest eigenvalue 
%                               of a given stochastic matrix based on the baseline 
%                               Distributed Orthogonal Iteration algortithm.
% estimate_l2_pwr_4             Estimate the second (by the modulus) largest eigenvalue 
%                               of a given stochastic matrix based on the modified and 
%                               streamlined Distributed Orthogonal Iteration algortithm 
%                               (Algorithm 1 in the paper)
% generate_opt_matrices_path    Script for generating the path graph adjacency matrix, 
%                               and various weights and weight matrices. Weight matrices 
%                               generated: Metropolis-Hastings, FDLA optimal. Weights 
%                               generated: Polynomial weights for 3, 5, 7 order polynomial filter
% generate_opt_matrices_RG      Script for generating the random geometric graph adjacency matrix, 
%                               and various weights and weight matrices. Weight matrices 
%                               generated: Metropolis-Hastings, FDLA optimal. Weights 
%                               generated: Polynomial weights for 3, 5, 7 order polynomial filter
% main_consensus_circ_lin_saved_mtx
%                               main script for the path graph. script for obtaining the MSE and 
%                               averaging time characterization of the proposed algorithm. Run this 
%                               script to obtain the simulation results used to create figures in 
%                               the paper. The simulation results are saved in a separate file for 
%                               a graph with the number of nodes N. The data saved by this script 
%                               can be furtheer accessed by the plotting routines.
% main_consensus_RG_saved_mtx   main script for the random geometric graph. script for obtaining the MSE and 
%                               averaging time characterization of the proposed algorithm. Run this 
%                               script to obtain the simulation results used to create figures in 
%                               the paper. The simulation results are saved in a separate file for 
%                               a graph with the number of nodes N. The data saved by this script 
%                               can be furtheer accessed by the plotting routines.
% plot_results_PATH_Slope       Script for plotting the figures. Path
%                               graph, slope initialization. Results for
%                               the proposed algorithm.
% plot_results_PATH_Slope_all   Script for plotting the figures. Path
%                               graph, slope initialization. Results all algorithms.
% plot_results_RGG_Slope        Script for plotting the figures. Random
%                               Geometric Graph graph, slope initialization. Results for
%                               the proposed algorithm.
% plot_results_RGG_Slope_all    Script for plotting the figures. Random
%                               Geometric Graph graph, slope initialization. Results for
%                               all algorithms.
% plot_results_RGG_Spike        Script for plotting the figures. Random
%                               Geometric Graph graph, spike initialization. Results for
%                               the proposed algorithm.