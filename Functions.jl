#
############################### Generic functions definitions ########################################
#
#
function euclidean_distance(values)
    return sqrt(sum(values.^2))
end
#
function uniform_box(n_value, x_length,y_length,z_length)
    r_atoms = Array{Float64}(undef,n_value,3)
    r_atoms[:,1] = (rand(n_value).-0.5).*x_length
    r_atoms[:,2] = (rand(n_value).-0.5).*y_length
    r_atoms[:,3] = (rand(n_value).-0.5).*z_length
    return r_atoms
end
#
function array_replace!(array,a,b)
    value_temp=array[a,:]
    array[a,:]=array[b,:]
    array[b,:]=value_temp
end
#
#
#
############################### Joint algorithm definitions ########################################
#
#
function sample_sphere_joint(n_start, distance_threshold, distance_function , sampling_function, sampling_options, time_threshold=60, print_option=false, check_option=false)
    #
    function println_if(string...)
        print_option ? println(string...) : nothing
    end
    #
    println_if("Starting evaluation core.")
    #
    n_value = n_start
    sampled_array = sampling_function(n_start, sampling_options...)
    #
    time_start=time()
    #
    #
    remain  = n_value
    i_old   = 0
    i_subs = n_value
    time_condition = true
    #
    while remain>0 && time_condition
        i_subs = n_value 
        #
        for i in 1:i_subs
            #
            start_point = max(i,i_old)
            #
            for j in (start_point + 1):i_subs
                while distance_function(sampled_array[i,:]-sampled_array[j,:])<distance_threshold && j<=i_subs            
                    array_replace!(sampled_array,j,i_subs)
                    i_subs-=1
                end
            end
        end
        #
        remaining_i = (i_subs+1):n_value
        sampled_array[remaining_i,:] = sampling_function(length(remaining_i), sampling_options...)
        remain = length(remaining_i)
        #
        println_if("Values remaining to sample: ",remain,".\nPoints already sampled: ",i_subs,".\nComputation_valuel time: ",time()-time_start)
        #
        i_old = i_subs
        #
        time_condition = (time()-time_start)<=time_threshold
        #
        if !time_condition
            @warn "Exceeded computation_value time (joint). Number of points correctly sampled: "*string(i_subs)
            n_value=i_subs
            sampled_array = sampled_array[1:n_value,:]
        end
        #
    end
    #
    core_computation_time=(time()-time_start)
    #
    if check_option
        remain=0
        println_if("Starting consistency check.")
        for i = 1:n_value
            for j = i + 1:n_value #i<j<=n_value
                if euclidean_distance(sampled_array[i,:]-sampled_array[j,:])<distance_threshold
                    remain=remain+1
                end
            end
        end
        if remain!=0
            @warn "Code failure. Number of exceptions: "*string(remain)
        end
    end
    #
    println_if("Number of exceptions: ",remain)
    println_if("Computation correctly completed.\nCore evaluation computation_valuel time: ", core_computation_time,"\nOverall computation_valuel time: ", time()-time_start)
    #
    return sampled_array, n_value 
end
#
#
#
############################### Single algorithm definitions ########################################
#
#
function sample_sphere_single(n_start, distance_threshold, distance_function , sampling_function, sampling_options, time_threshold=60, print_option=false, check_option=false)
    #
    function println_if(string...)
        print_option ? println(string...) : nothing
    end
    #    
    sampled_array = sampling_function(1, sampling_options...)
    time_condition = true
    time_start = time()
    n_value = 1
    #
    while n_value<n_start && time_condition
        keep_sampling = true
        while keep_sampling
            keep_sampling = false
            sampled_point = sampling_function(1, sampling_options...)
            for i in 1:n_value
                if distance_function(sampled_array[i,:]-sampled_point[1,:]) < distance_threshold
                    keep_sampling = true
                    break
                end
            end
            if !keep_sampling
                sampled_array = vcat(sampled_array,sampled_point)
            end
        end
        n_value = length(sampled_array[:,1])
        #
        time_condition = (time()-time_start)<=time_threshold
        if !time_condition
            @warn "Exceeded computation_value time (single). Number of atoms correctly sampled: "*string(n_start-n_value)
            sampled_array = sampled_array[1:n_value,:]
        end
    end
    #
    if check_option
        remain=0
        println_if("Starting consistency check.")
        for i = 1:n_value
            for j = i + 1:n_value 
                if euclidean_distance(sampled_array[i,:]-sampled_array[j,:])<distance_threshold
                    remain=remain+1
                end
            end
        end
        if remain!=0
            @warn "Code failure. Number of exceptions: "*string(remain)
        end
    end
    #
    return sampled_array, n_value 
end
#
#
#
############################### Smart algorithm definitions ########################################
#
#
function region_shape_cube(vec, sizes)
    dims = length(sizes)
    valid = true
    for i in 1:dims
        if vec[i]<0||vec[i]>sizes[i]
            valid = false
            break
        end
    end
    return valid
end
#
function valid_neighbour(i_neigh,dim_points )
    valid = true
    for i in 1:length(i_neigh)
        if i_neigh[i]<1 || i_neigh[i]>dim_points[i]
            valid = false
            break
        end
    end
    return valid
end
#
function sample_sphere_smart(n_start, distance_threshold, distance_function, region_shape_func, enclosing_cube::Tuple, time_threshold=60, print_option=false, check_option=false)
    #
    #
    dimensions = length(enclosing_cube)
    #
    neigh_order = [0]
    for i in 1:ceil(Int,sqrt(dimensions))
        push!(neigh_order,i)
        push!(neigh_order,-i)
    end
    #
    neighbour_list = collect(Iterators.product([neigh_order for i in 1:dimensions]... ))[2:end]
    #
    small_cube_size = distance_threshold/sqrt(dimensions)
    #
    function println_if(string...)
        print_option ? println(string...) : nothing
    end
    #
    #Construct small cubes where only one point is allowed to be
    dim_points = (x->ceil(Int, x/small_cube_size)).(enclosing_cube)
    #
    grid_points  = Array{Float64}(undef, dim_points..., 3)
    grid_points .= 0.0
    #
    grid_indices_values = collect(Iterators.product((x->collect(1:x)).(dim_points)... ))
    grid_indices_availability  = Array{Bool}(undef, dim_points...)
    grid_indices_availability .= true
    #
    n_value = 0
    #
    #
    time_condition = true
    time_start = time()
    while n_value<n_start && time_condition
        point_rejection = true
        while point_rejection && time_condition
            #
            i_list = rand((grid_indices_values[grid_indices_availability])[:])
            sampled_point = (i_list .- rand(dimensions)).*small_cube_size
            point_rejection = !region_shape_func(sampled_point)
            #
            if !point_rejection  
                #
                i_near_list = (x->i_list.+x).(neighbour_list)
                i_near_list = i_near_list[(x->valid_neighbour(x,dim_points)).(i_near_list)]
                i_near_list = i_near_list[(x->grid_points[x...,:]!=[0.0;0.0;0.0]).(i_near_list)]
                #
                for i_near in i_near_list
                    if distance_function(grid_points[i_near...,:]-sampled_point)<distance_threshold
                        point_rejection = true
                        break
                    end
                end
            end
            #
            if !point_rejection
                grid_points[i_list...,:] = sampled_point[:]
                grid_indices_availability[i_list...] = false
                n_value+=1
            end
            time_condition = (time()-time_start)<=time_threshold
        end
    end
    #
    #
    correct_indices = (grid_indices_values[.!grid_indices_availability])[:]
    #
    if !time_condition
        @warn "Exceeded computation_value time (smart). Number of points correctly sampled: "*string(n_value)
    else
        length(correct_indices[:]) != n_start || n_value!=n_start ? error("Wrong amount of sampled spheres.") : nothing
    end
    #    
    sampled_array = Array{Float64}(undef, n_value, 3)
    for i in 1:n_value
        sampled_array[i,:] = grid_points[correct_indices[i]...,:]
    end
    #
    #
    if check_option
        remain=0
        println_if("Starting consistency check.")
        for i = 1:n_value
            for j = (i + 1):n_value 
                if euclidean_distance(sampled_array[i,:]-sampled_array[j,:])<distance_threshold
                    remain+=1
                    println_if("Exception at: ", sampled_array[i,:], " and ", sampled_array[j,:])
                end
            end
        end
        if remain!=0
            @warn "Code failure. Number of exceptions: "*string(remain)
        end
        println_if("End consistency check.")
    end
    #
    return sampled_array, n_value 
end
#
#