#
using Plots
using CurveFit, Statistics, CPUTime, LaTeXStrings, Distributions, HypothesisTests
include("Functions.jl")
ENV["GKSwstype"]="nul"
#
#
############################### Test parameters ########################################
#
#
verbose_option = [true ; false][2]
#
x_length = 1.
y_length = 1.
z_length = 1.
#
filling_ratio = 0.4
n_array = ceil.(Int, 10.0.^[ 1.7;1.8; 1.9 ; 2; 2.1 ; 2.2; 2.3;2.4;2.5 ;2.6;2.7;2.8;2.9])
#
n_repetitions = 1#100#15000
max_time = 60
#
const sphere_filling_option   = [true ; false][1]
const code_test_option        = [true ; false][1]
const p_value_option          = [true ; false][1] 
# 
const basic_code_option       = [true ; false][1]
const joint_code_option       = [true ; false][2]
const grid_code_option        = [true ; false][1]
const grid_approx_code_option = [true ; false][1]
#
############################### Initializing elements ########################################
#
#
#Initializing values
n_length = length(n_array)
cube_values = (x_length,y_length,z_length)
#
if sphere_filling_option
    max_filling = pi/(3*sqrt(2))
    distance_func = (nn-> 2*((filling_ratio*max_filling*x_length*y_length*z_length / (4*nn*pi/3))^(1/3)) )
else
    distance_func = (xx-> ((x_length*y_length*z_length/xx)^(1/3))*filling_ratio)
end
#
distance_array = distance_func.(n_array)
length(ARGS)>=2 ? args_checked=ARGS[:] : args_checked=["_DEFAULT" ; ""]
#
#
#Initializing functions
function test_results(distance_func,sampled_array, distance_threshold, name_code ; error_option = true)
    n_value = length(sampled_array[:, 1 ])
    errors = 0
    for i in 1:n_value, j in (i + 1):n_value
        if distance_func(sampled_array[i, : ]-sampled_array[j,:]) < distance_threshold
            errors+=1
        end
    end
    if errors>0
        exception_string = "-- Problem with n_val = "*string(n_value)*"("*name_code*"). Found "*string(errors)*" exceptions."
        if error_option
            @error exception_string
            flush(stdout)
            return false
        else
            @warn exception_string
        end
    end
    return true
end
function nice_out_print(name_print, time_array, time_array_std, time_array_std_log)
    to_add_string = " "^length(name_print)
    println(name_print*" code times: ", dig_round.(time_array))
    println(to_add_string*"         +/- ", dig_round.(time_array_std))
    println(to_add_string*"   (log) +/- ", dig_round.(time_array_std_log))
end
function fit_log10(time_array, n_array, start_i_fit)
    pure_index = (x->time_array[x]!=0.0&&x>=start_i_fit).(1:length(n_array))
    fit_val = linear_fit(log10.(n_array[pure_index]), log10.(time_array[pure_index]))
    return fit_val
end
#
function p_extract(full_array, distr_func=Uniform())
    dimensions = length(full_array[1,:])
    p_value = 0
    for i in 1:dimensions
        p_value += pvalue(ExactOneSampleKSTest(full_array[:,i] ,distr_func))
    end
    return p_value/dimensions
end
#
#
#Initializing arrays
if basic_code_option
    exit_array_basic = [true for i in 1:n_length]
    time_array_basic = zeros(n_length)
    time_array_basic_std = zeros(n_length)
    time_array_basic_std_log = zeros(n_length)
    time_array_basic_list = zeros(n_length,n_repetitions)
    p_value_basic_list = zeros(n_length,n_repetitions)
    p_value_basic_mean = zeros(n_length)
    @time sample_sphere_basic(n_array[1], distance_array[1], euclidean_distance , uniform_box, cube_values, 5, false, true)
end
#
if joint_code_option
    exit_array_joint = [true for i in 1:n_length]
    time_array_joint = zeros(n_length)
    time_array_joint_std = zeros(n_length)
    time_array_joint_std_log = zeros(n_length)
    time_array_joint_list = zeros(n_length,n_repetitions)
    p_value_joint_list = zeros(n_length,n_repetitions)
    p_value_joint_mean = zeros(n_length)
    @time sample_sphere_joint(n_array[1], distance_array[1], euclidean_distance , uniform_box, cube_values, 5, false, true)
end
#
if grid_code_option
    exit_array_grid = [true for i in 1:n_length]
    time_array_grid = zeros(n_length)
    time_array_grid_std = zeros(n_length)
    time_array_grid_std_log = zeros(n_length)
    time_array_grid_list = zeros(n_length,n_repetitions)
    p_value_grid_list = zeros(n_length,n_repetitions)
    p_value_grid_mean = zeros(n_length)
    @time sample_sphere_grid(n_array[1], distance_array[1], euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values,5, false, true ; approx=true, approx_eff = true)
end
#
if grid_code_option
    exit_array_grid_approx = [true for i in 1:n_length]
    time_array_grid_approx = zeros(n_length)
    time_array_grid_approx_std = zeros(n_length)
    time_array_grid_approx_std_log = zeros(n_length)
    time_array_grid_approx_list = zeros(n_length,n_repetitions)
    p_value_grid_approx_list = zeros(n_length,n_repetitions)
    p_value_grid_approx_mean = zeros(n_length)
    @time sample_sphere_grid(n_array[1], distance_array[1], euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values,5, false, true ; approx=false, approx_eff = false)
end
#
#
#
############################### Testing the algorithms ########################################
#
#
println("The number of repetitions is ", n_repetitions,".\nThe filling fraction is ", filling_ratio,".\nThe number of points is ", n_array )
println("The filling factor is defined from the sphere volume? ", sphere_filling_option)
flush(stdout)
#
time_start = time()
#
for i_main in 1:length(n_array)
    #
    time_start_for = time()
    #
    n_main = n_array[i_main]
    distance_threshold = distance_array[i_main]
    if verbose_option
        println("Starting n_main = ", n_main, ". Distance threshold: ", distance_threshold)
    end
    for i_rep in 1:n_repetitions
        #
        #
        if verbose_option
            println("-- Starting repetition ", i_rep)
            flush(stdout)
        end
        #
        #Computing the codes & timing them
        if basic_code_option
            CPUtic()
            sampled_array_basic_val, n_value_basic_val = sample_sphere_basic(n_main, distance_threshold, euclidean_distance , uniform_box, (x_length,y_length,z_length), max_time)
            time_array_basic_list[i_main,i_rep] = CPUtoq()
            exit_array_basic[i_main] = test_results(euclidean_distance,sampled_array_basic_val, distance_threshold, "Basic")
            if p_value_option
                p_value_basic_list[i_main,i_rep] = p_extract(sampled_array_basic_val)
            end
        end
        if joint_code_option
            CPUtic()
            sampled_array_joint_val, n_value_joint_val = sample_sphere_joint(n_main, distance_threshold, euclidean_distance , uniform_box, (x_length,y_length,z_length), max_time)
            time_array_joint_list[i_main,i_rep] = CPUtoq()
            exit_array_joint[i_main] = test_results(euclidean_distance,sampled_array_joint_val, distance_threshold, "Joint")
            if p_value_option
                p_value_joint_list[i_main,i_rep] = p_extract(sampled_array_joint_val)
            end
        end
        if grid_code_option
            CPUtic()
            sampled_array_grid_val, n_value_grid_val = sample_sphere_grid(n_main, distance_threshold, euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values, max_time; approx=false)
            time_array_grid_list[i_main,i_rep] = CPUtoq()
            exit_array_grid[i_main] = test_results(euclidean_distance,sampled_array_grid_val, distance_threshold, "Grid")       
            if p_value_option
                p_value_grid_list[i_main,i_rep] = p_extract(sampled_array_grid_val)
            end
        end
        if grid_approx_code_option
            CPUtic()
            sampled_array_grid_approx_val, n_value_grid_approx_val = sample_sphere_grid(n_main, distance_threshold, euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values, max_time; approx=true)
            time_array_grid_approx_list[i_main,i_rep] = CPUtoq()
            exit_array_grid_approx[i_main] = test_results(euclidean_distance,sampled_array_grid_approx_val, distance_threshold, "Grid, approx")    
            if p_value_option
                p_value_grid_approx_list[i_main,i_rep] = p_extract(sampled_array_grid_approx_val)
            end
        end
        #
        #
    end
    #
    #
    #Saving the data
    if basic_code_option
        time_array_basic[i_main]= mean(time_array_basic_list[i_main,:])
        time_array_basic_std[i_main]= std(time_array_basic_list[i_main,:])
        time_array_basic_std_log[i_main]= std(log.(time_array_basic_list[i_main,:]))
        if p_value_option
            p_value_basic_mean[i_main] = mean(p_value_basic_list[i_main,:])
        end
    end
    if joint_code_option
        time_array_joint[i_main] = mean(time_array_joint_list[i_main,:])
        time_array_joint_std[i_main] = std(time_array_joint_list[i_main,:])
        time_array_joint_std_log[i_main] = std(log.(time_array_joint_list[i_main,:]))
        if p_value_option
            p_value_joint_mean[i_main] = mean(p_value_joint_list[i_main,:])
        end
    end
    if grid_code_option   
        time_array_grid[i_main]= mean(time_array_grid_list[i_main,:])    
        time_array_grid_std[i_main]= std(time_array_grid_list[i_main,:])
        time_array_grid_std_log[i_main]= std(log.(time_array_grid_list[i_main,:]))
        if p_value_option
            p_value_grid_mean[i_main] = mean(p_value_grid_list[i_main,:])
        end
    end
    if grid_approx_code_option  
        time_array_grid_approx[i_main]= mean(time_array_grid_approx_list[i_main,:])    
        time_array_grid_approx_std[i_main]= std(time_array_grid_approx_list[i_main,:])
        time_array_grid_approx_std_log[i_main]= std(log.(time_array_grid_approx_list[i_main,:]))
        if p_value_option
            p_value_grid_approx_mean[i_main] = mean(p_value_grid_approx_list[i_main,:])
        end
    end
    #
    #
    if verbose_option 
        println("Finished n_main = ", n_main, " in ", time()-time_start_for)
    end
    #
    #
end
#
#
#
#Printing & elaborating values
#
dig_round = x-> round(x; digits=4)
start_i_fit = 2
start_i_plot = 2
#
println("\nEvaluation concluded in ", time()-time_start," seconds.\n")
println("N values: ", n_array, ", repetitions: ",n_repetitions,", filling: ", filling_ratio,"\n")
#
if sphere_filling_option
    filling_string ="f"*L"_{\mathrm{max}}"
else
    filling_string=""
end
#
plot()
ylabel!(L"log_{10}(\mathrm{CPU}\;\mathrm{time})")
xlabel!(L"log_{10}(\mathrm{N}\;\mathrm{spheres})")
title!("\nFilling: f = "*string(filling_ratio)*filling_string*", repetitions: "*string(n_repetitions))
(x_start, x_end) = Tuple((log10.(n_array[start_i_plot:end]))[[1;end]])
#
if basic_code_option
    nice_out_print("Basic", time_array_basic, time_array_basic_std, time_array_basic_std_log)
    fit_basic = fit_log10(time_array_basic, n_array, start_i_fit) 
    println("--- fit_basic = ", fit_basic,"\n")
    plot!( log10.(n_array[start_i_plot:end]),  log10.(time_array_basic[start_i_plot:end]), label="Basic", seriestype=:scatter, framestyle = :box, color=:red)
    #
    fit_f(x) = fit_basic[1]+x*fit_basic[2]
    plot!(fit_f  ,x_start,  x_end, framestyle = :box, linestyle=:dash, color=:red, label="")
    #
end
if joint_code_option
    nice_out_print("Joint", time_array_joint, time_array_joint_std, time_array_joint_std_log)
    fit_joint = fit_log10(time_array_joint, n_array, start_i_fit) 
    println("--- fit_joint = ", fit_joint,"\n")
    plot!( log10.(n_array[start_i_plot:end]),  log10.(time_array_joint[start_i_plot:end]), label="Joint", seriestype=:scatter, framestyle = :box, color=:orange)
    #
    fit_f(x) = fit_joint[1]+x*fit_joint[2]
    plot!(fit_f  ,x_start,  x_end, framestyle = :box, linestyle=:dash, color=:orange, label="")
    #
end
if grid_code_option
    nice_out_print("Grid", time_array_grid, time_array_grid_std, time_array_grid_std_log)
    fit_grid = fit_log10(time_array_grid, n_array, start_i_fit) 
    println("--- fit_grid = ", fit_grid,"\n")
    plot!( log10.(n_array[start_i_plot:end]),  log10.(time_array_grid[start_i_plot:end]), label="Grid", seriestype=:scatter, framestyle = :box, color=:blue)
    #
    fit_f(x) = fit_grid[1]+x*fit_grid[2]
    plot!(fit_f  ,x_start,  x_end, framestyle = :box, linestyle=:dash, color=:blue, label="")
    #
end
if grid_approx_code_option
    nice_out_print("Grid, approx", time_array_grid_approx, time_array_grid_approx_std, time_array_grid_approx_std_log)
    fit_grid_approx = fit_log10(time_array_grid_approx, n_array, start_i_fit) 
    println("--- fit_grid_approx = ", fit_grid_approx,"\n")
    plot!( log10.(n_array[start_i_plot:end]),  log10.(time_array_grid_approx[start_i_plot:end]), label="Grid, approx", seriestype=:scatter, framestyle = :box, color=:green)
    #
    fit_f(x) = fit_grid_approx[1]+x*fit_grid_approx[2]
    plot!(fit_f  ,x_start,  x_end, framestyle = :box, linestyle=:dash, color=:green, label="")
    #
end
#
#
mkpath("Data")
png("Data/TIME_fill"*string(filling_ratio)*"_reps"*string(n_repetitions)*"_"*args_checked[1]*"_"*args_checked[2]*".png")
plot!()
#
if p_value_option
    plot()
    plot_type_choice = 2#1, 2
    func_plot = [log10; x->x][plot_type_choice]
    counter_func_plot = [x->x ; x->10^x][plot_type_choice]
    x_start_here,  x_end_here = counter_func_plot.((x_start,  x_end))
    ylabel!([L"log_{10}(\langle\mathrm{p}_{\mathrm{uniform}}\rangle)";L"\langle\mathrm{p}_{\mathrm{uniform}}\rangle"][plot_type_choice])
    xlabel!([L"log_{10}(\mathrm{N}\;\mathrm{spheres})";L"\mathrm{N}\;\mathrm{spheres}"][plot_type_choice])
    title!("\nFilling: f = "*string(filling_ratio)*filling_string*", repetitions: "*string(n_repetitions))
    plot!( xx->func_plot(0.05),x_start_here,  x_end_here, framestyle = :box, linestyle=:dash, color=:black, label="<p> = 0.05", ylims=(0,1))
    plot!( xx->func_plot(0.5),x_start_here,  x_end_here, framestyle = :box,   color=:black, label="<p> = 0.5", ylims=(0,1))
    #
    if basic_code_option
        plot!( func_plot.(n_array[start_i_plot:end]),  func_plot.(p_value_basic_mean[start_i_plot:end]), label="Basic", seriestype=:scatter, framestyle = :box, color=:red)
    end
    if joint_code_option
        plot!( func_plot.(n_array[start_i_plot:end]),  func_plot.(p_value_joint_mean[start_i_plot:end]), label="Joint", seriestype=:scatter, framestyle = :box, color=:orange)
    end
    if grid_code_option
        plot!( func_plot.(n_array[start_i_plot:end]),  func_plot.(p_value_grid_mean[start_i_plot:end]), label="Grid", seriestype=:scatter, framestyle = :box, color=:blue)
    end
    if grid_approx_code_option
        plot!( func_plot.(n_array[start_i_plot:end]),  func_plot.(p_value_grid_approx_mean[start_i_plot:end]), label="Grid, approx", seriestype=:scatter, framestyle = :box, color=:green)
    end
end
#
#
png("Data/P_VALUE_fill"*string(filling_ratio)*"_reps"*string(n_repetitions)*"_"*args_checked[1]*"_"*args_checked[2]*".png")
plot!()