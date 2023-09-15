#
using Plots
using CurveFit, Statistics, CPUTime
include("Functions.jl")
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
filling = 0.8#0.95
n_array = [10, 20,  40 , 60 ,80, 100, 200]#; 1000]
#
n_repetitions = 200#100000
max_time = 3000#180
#
#
############################### Initializing elements ########################################
#
#
n_length = length(n_array)
cube_values = (x_length,y_length,z_length)
distance_func = (xx-> ((x_length*y_length*z_length/xx)^(1/3))*filling)
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
sample_sphere_joint(n_array[1], distance_array[1], euclidean_distance , uniform_box, cube_values, max_time, false, true)
sample_sphere_single(n_array[1], distance_array[1], euclidean_distance , uniform_box, cube_values, max_time, false, true)
sample_sphere_smart(n_array[1], distance_array[1], euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values, max_time, false, true)
#
#
############################### Testing the algorithms ########################################
#
#
println("The number of repetitions is ", n_repetitions,".\nThe filling fraction is ", filling,".\nThe number of points is ", n_array )
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
        sampled_array_smart_val, n_value_smart_val = sample_sphere_smart(n_main, distance_threshold, euclidean_distance, x->(region_shape_cube(x,cube_values)), cube_values, max_time)
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
println("\nEvaluation concluded in ", time()-time_start," seconds.\n")
println("N values: ", n_array, ", repetitions: ",n_repetitions,", filling: ", filling)
println("joint times: ", time_array_joint)
println("         +/- ", time_array_joint_std)
println("   (log) +/- ", time_array_joint_std_log)
println("")
println("single times: ", time_array_single)
println("         +/- ", time_array_single_std)
println("   (log) +/- ", time_array_single_std_log)
println("")
println("smart times: ", time_array_smart)
println("         +/- ", time_array_smart_std)
println("   (log) +/- ", time_array_smart_std_log)
#
#
#Fitting the scaling
start_i_fit = 2
fit_joint = linear_fit(log.(n_array[start_i_fit:end]), log.(time_array_joint[start_i_fit:end]))
fit_single = linear_fit(log.(n_array[start_i_fit:end]), log.(time_array_single[start_i_fit:end]))
fit_smart = linear_fit(log.(n_array[start_i_fit:end]), log.(time_array_smart[start_i_fit:end]))
println("fit_joint = ", fit_joint)
println("fit_single = ", fit_single)
println("fit_smart = ", fit_smart)
#
#
#Final plot
start_i_plot = 2
plot(log.(n_array[start_i_plot:end]), log.(time_array_joint[start_i_plot:end]), label="joint", seriestype=:scatter, yerror=time_array_joint_std_log[start_i_plot:end])
plot!(log.(n_array[start_i_plot:end]), log.(time_array_single[start_i_plot:end]), label="single", seriestype=:scatter, yerror=time_array_single_std_log[start_i_plot:end])
plot!(log.(n_array[start_i_plot:end]), log.(time_array_smart[start_i_plot:end]), label="smart", seriestype=:scatter, yerror=time_array_smart_std_log[start_i_plot:end])
ylabel!("log(T)")
xlabel!("log(N)")
title!("Filling at "*string(filling)*" and "*string(n_repetitions)*" repetitions")
mkpath("Data")
png("Data/fill"*string(filling)*"_rep"*string(n_repetitions)*"_"*args_checked[1]*"_"*args_checked[2]*".png")