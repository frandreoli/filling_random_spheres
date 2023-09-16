#
using Plots
using CurveFit, Statistics, CPUTime, LaTeXStrings
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
filling_ratio = 0.6#0.95
n_array = ceil.(Int, 10.0.^[ 1.7;1.8; 1.9 ; 2; 2.1 ; 2.2; 2.3;2.4 ])#[ 40 , 60 ,80, 100, 200]#; 1000]
#
n_repetitions = 100#20000
max_time = 30#180
#
approx_option = true
#
############################### Initializing elements ########################################
#
#
n_length = length(n_array)
cube_values = (x_length,y_length,z_length)
#distance_func = (xx-> ((x_length*y_length*z_length/xx)^(1/3))*filling_ratio)
#
#max_filling = pi/(3*sqrt(2))
#distance_func = (nn-> 2*((filling_ratio*max_filling*x_length*y_length*z_length / (4*nn*pi/3))^(1/3)) )
distance_func = (xx-> ((x_length*y_length*z_length/xx)^(1/3))*filling_ratio)
#
distance_array = distance_func.(n_array)
length(ARGS)>=2 ? args_checked=ARGS[:] : args_checked=["_DEFAULT" ; ""]
#
#
exit_array_joint = [true for i in 1:n_length]
exit_array_single = [true for i in 1:n_length]
exit_array_smart = [true for i in 1:n_length]
#
time_array_joint = zeros(n_length)
time_array_single = zeros(n_length)
time_array_smart = zeros(n_length)
#
time_array_joint_std = zeros(n_length)
time_array_single_std = zeros(n_length)
time_array_smart_std = zeros(n_length)
#
time_array_joint_std_log = zeros(n_length)
time_array_single_std_log = zeros(n_length)
time_array_smart_std_log = zeros(n_length)
#
time_array_joint_list = zeros(n_length,n_repetitions)
time_array_single_list = zeros(n_length,n_repetitions)
time_array_smart_list = zeros(n_length,n_repetitions)
#
@time sample_sphere_joint(n_array[1], distance_array[1], euclidean_distance , uniform_box, cube_values, 5, false, true)
@time sample_sphere_single(n_array[1], distance_array[1], euclidean_distance , uniform_box, cube_values, 5, false, true)
@time sample_sphere_smart(n_array[1], distance_array[1], euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values,5, false, true ; approx=approx_option, approx_eff = approx_option)
#
#
############################### Testing the algorithms ########################################
#
#
println("The number of repetitions is ", n_repetitions,".\nThe filling fraction is ", filling_ratio,".\nThe number of points is ", n_array )
println("The approximation option is ",  approx_option)
flush(stdout)
#
time_start = time()
#
for i_main in 1:length(n_array)
    #
    n_main = n_array[i_main]
    distance_threshold = distance_array[i_main]
    if verbose_option
        println("Starting n_main = ", n_main, ". Distance threshold: ", distance_threshold)
    end
    for i_rep in 1:n_repetitions
        if verbose_option
            println("-- Starting repetition ", i_rep)
            flush(stdout)
        end
        #
        #
        CPUtic()
        sampled_array_joint_val, n_value_joint_val = sample_sphere_joint(n_main, distance_threshold, euclidean_distance , uniform_box, (x_length,y_length,z_length), max_time)
        time_array_joint_list[i_main,i_rep] = CPUtoq()
        #
        CPUtic()
        sampled_array_single_val, n_value_single_val = sample_sphere_single(n_main, distance_threshold, euclidean_distance , uniform_box, (x_length,y_length,z_length), max_time)
        time_array_single_list[i_main,i_rep] = CPUtoq()
        #
        CPUtic()
        sampled_array_smart_val, n_value_smart_val = sample_sphere_smart(n_main, distance_threshold, euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values, max_time; approx=approx_option)
        time_array_smart_list[i_main,i_rep] = CPUtoq()
        #
        #Testing the joint code
        errors_joint = 0
        for i in 1:n_value_joint_val, j in (i + 1):n_value_joint_val
            if euclidean_distance(sampled_array_joint_val[i, : ]-sampled_array_joint_val[j,:]) < distance_threshold
                errors_joint+=1
            end
        end
        if errors_joint>0
            error("-- Problem with n_val = ", n_main," (joint). Found ", errors_joint," exceptions.")
            flush(stdout)
            exit_array_joint[i_main]= false
        end
        #
        #Testing the single code
        errors_single = 0
        for i in 1:n_value_single_val, j in (i + 1):n_value_single_val 
            if euclidean_distance(sampled_array_single_val[i, : ]-sampled_array_single_val[j,:]) < distance_threshold
                errors_single+=1
            end
        end
        if errors_single>0
            error("-- Problem with n_val = ", n_main," (single). Found ", errors_single," exceptions.")
            flush(stdout)
            exit_array_single[i_main]= false
        end
        #
        #Testing the smart code
        errors_smart = 0
        for i in 1:n_value_smart_val, j in (i + 1):n_value_smart_val 
            if euclidean_distance(sampled_array_smart_val[i, : ]-sampled_array_smart_val[j,:]) < distance_threshold
                errors_smart+=1
            end
        end
        if errors_smart>0
            error("-- Problem with n_val = ", n_main," (smart). Found ", errors_smart," exceptions.")
            flush(stdout)
            exit_array_smart[i_main]= false
        end        
        #
    end
    #
    time_array_joint[i_main] = mean(time_array_joint_list[i_main,:])
    time_array_single[i_main]= mean(time_array_single_list[i_main,:])
    time_array_smart[i_main]= mean(time_array_smart_list[i_main,:])
    time_array_joint_std[i_main] = std(time_array_joint_list[i_main,:])
    time_array_single_std[i_main]= std(time_array_single_list[i_main,:])
    time_array_smart_std[i_main]= std(time_array_smart_list[i_main,:])
    time_array_joint_std_log[i_main] = std(log.(time_array_joint_list[i_main,:]))
    time_array_single_std_log[i_main]= std(log.(time_array_single_list[i_main,:]))
    time_array_smart_std_log[i_main]= std(log.(time_array_smart_list[i_main,:]))
    #
    if verbose_option 
        println("Finished.\n")
    end
    #
    #
end
#
#
#
#Printing values
#
dig_round = x-> round(x; digits=4)
#
println("\nEvaluation concluded in ", time()-time_start," seconds.\n")
println("N values: ", n_array, ", repetitions: ",n_repetitions,", filling: ", filling_ratio)
println("joint times: ", dig_round.(time_array_joint))
println("         +/- ", dig_round.(time_array_joint_std))
println("   (log) +/- ", dig_round.(time_array_joint_std_log))
println("")
println("single times: ", dig_round.(time_array_single))
println("         +/- ", dig_round.(time_array_single_std))
println("   (log) +/- ", dig_round.(time_array_single_std_log))
println("")
println("smart times: ", dig_round.(time_array_smart))
println("         +/- ", dig_round.(time_array_smart_std))
println("   (log) +/- ", dig_round.(time_array_smart_std_log))
#
#
#
#
#Fitting the scaling
start_i_fit = 2
#
pure_index_joint = (x->time_array_joint[x]!=0.0&&x>=start_i_fit).(1:length(n_array))
fit_joint = linear_fit(log10.(n_array[pure_index_joint]), log10.(time_array_joint[pure_index_joint]))
#
pure_index_single = (x->time_array_single[x]!=0.0&&x>=start_i_fit).(1:length(n_array))
fit_single = linear_fit(log10.(n_array[pure_index_single]), log10.(time_array_single[pure_index_single]))
#
pure_index_smart = (x->time_array_joint[x]!=0.0&&x>=start_i_fit).(1:length(n_array))
fit_smart = linear_fit(log10.(n_array[pure_index_smart]), log10.(time_array_smart[pure_index_smart]))
#
println("")
println("fit_joint = ", fit_joint)
println("fit_single = ", fit_single)
println("fit_smart = ", fit_smart)
#
#
if approx_option
    approx_string = ", approx."
else
    approx_string = ""
end
#
#Final plot
start_i_plot = 2
plot( log10.(n_array[start_i_plot:end]),  log10.(time_array_joint[start_i_plot:end]), label="joint", seriestype=:scatter, framestyle = :box)#, yerror=time_array_joint_std_log[start_i_plot:end] )
plot!( log10.(n_array[start_i_plot:end]),  log10.(time_array_single[start_i_plot:end]), label="single", seriestype=:scatter, framestyle = :box)#, yerror=time_array_single_std_log[start_i_plot:end] )
plot!( log10.(n_array[start_i_plot:end]),  log10.(time_array_smart[start_i_plot:end]), label="smart", seriestype=:scatter, framestyle = :box)#, yerror=time_array_smart_std_log[start_i_plot:end] )
ylabel!(L"log_{10}(\mathrm{CPU}\;\mathrm{time})")
xlabel!(L"log_{10}(\mathrm{N}\;\mathrm{spheres})")
title!("Filling ratio: "*string(filling_ratio)*", repetitions: "*string(n_repetitions)*approx_string)
mkpath("Data")
png("Data/fill"*string(filling_ratio)*"_rep"*string(n_repetitions)*"_"*args_checked[1]*"_"*args_checked[2]*".png")
plot!()