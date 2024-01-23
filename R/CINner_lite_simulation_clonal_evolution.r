simulation_clonal_evolution <- function(t_end_time = 0,
                                        n_selective_clones = 1,
                                        vec_s_mut_time = 0,
                                        vec_s_mut_hierarchy = 0,
                                        vec_propensity = 0,
                                        vec_propensity_sum = 0,
                                        range_clonal_perc = 0,
                                        range_population = 0,
                                        t_tau_step = 0) {
    #---------------------Create one simulation for the clonal evolution
    #   Prepare the clonal populations: initial neutral clone has 1 cell
    #   and every other clone starts with 0 cell
    vec_populations <- rep(0, length = (n_selective_clones + 1))
    vec_populations[1] <- 1
    #   Prepare the record of time points
    t_previous <- 0
    t_now <- 0
    i_next_selective_mut <- 1
    t_next_selective_mut <- vec_s_mut_time[i_next_selective_mut]
    #   Prepare the record of clonal evolution
    record_t_previous <- c()
    record_t_now <- c()
    record_vec_populations <- matrix(
        nrow = 0,
        ncol = (n_selective_clones + 1)
    )
    record_vec_count_division <- matrix(
        nrow = 0,
        ncol = (n_selective_clones + 1)
    )
    record_vec_count_death <- matrix(
        nrow = 0,
        ncol = (n_selective_clones + 1)
    )
    #   Simulation for the clonal evolution
    while (t_now < t_end_time) {
        #-------------------------------------------Find next time point
        t_previous <- t_now
        t_now <- t_now + t_tau_step
        if (t_now > t_end_time) {
            #   Wrap up the simulation if final time is exceeded
            t_now <- t_end_time
            next
        } else if ((t_previous <= t_next_selective_mut) & (t_now > t_next_selective_mut)) {
            #   Next selection mutation arrives...
            #   Update clonal populations
            clone_next <- i_next_selective_mut
            vec_populations[clone_next + 1] <- 1
            # if (i_next_selective_mut > 1) {
            #     clone_next <- i_next_selective_mut - 1
            #     vec_populations[clone_next + 1] <- 1
            # }
            #   Update time point of selective sweeps
            i_next_selective_mut <- i_next_selective_mut + 1
            if (i_next_selective_mut <= length(vec_s_mut_time)) {
                t_next_selective_mut <- vec_s_mut_time[i_next_selective_mut]
            } else {
                t_next_selective_mut <- Inf
            }
        }
        #---Find the Poisson propensities of event counts for all clones
        vec_propensity_division <- t_tau_step * vec_populations * vec_propensity[1, ]
        vec_propensity_death <- t_tau_step * vec_populations * vec_propensity[2, ]
        #---Find event counts for each clone
        vec_count_division <- rep(Inf, length = (n_selective_clones + 1))
        vec_count_death <- rep(Inf, length = (n_selective_clones + 1))
        while (max(vec_count_division + vec_count_death -
            vec_populations) > 0) {
            vec_count_division <- rpois(
                n = (n_selective_clones + 1),
                lambda = vec_propensity_division
            )
            vec_count_death <- rpois(
                n = (n_selective_clones + 1),
                lambda = vec_propensity_death
            )
        }
        #---Update the clonal populations
        vec_populations <- vec_populations +
            vec_count_division - vec_count_death
        #---Update record of clonal evolution
        record_t_previous <- c(record_t_previous, t_previous)
        record_t_now <- c(record_t_now, t_now)
        record_vec_populations <- rbind(
            record_vec_populations,
            vec_populations
        )
        record_vec_count_division <- rbind(
            record_vec_count_division,
            vec_count_division
        )
        record_vec_count_death <- rbind(
            record_vec_count_death,
            vec_count_death
        )
    }
    #------------------------Decide the quality flag of clonal evolution
    if (max(vec_populations) == 0) {
        condition_population <- 0
        condition_cancer <- 0
        simulation <- list()
        output <- list()
        output$condition_population <- condition_population
        output$condition_cancer <- condition_cancer
        output$simulation <- simulation
        return(output)
    }
    vec_population_perc <- 100 * vec_populations / sum(vec_populations)
    #   Check that selective clones' populations are in wanted range
    if ((min(vec_population_perc) < range_clonal_perc[1]) |
        (max(vec_population_perc) > range_clonal_perc[2])) {
        condition_cancer <- 0
    } else {
        condition_cancer <- 1
    }
    #   Check that total population size is in wanted range
    if ((sum(vec_populations) < range_population[1]) | (sum(vec_populations) > range_population[2])) {
        condition_population <- 0
    } else {
        condition_population <- 1
    }
    #---------------------------------------------Prepare output package
    simulation <- list()
    simulation$record_t_previous <- record_t_previous
    simulation$record_t_now <- record_t_now
    simulation$record_vec_populations <- record_vec_populations
    simulation$record_vec_count_division <- record_vec_count_division
    simulation$record_vec_count_death <- record_vec_count_death

    output <- list()
    output$condition_population <- condition_population
    output$condition_cancer <- condition_cancer
    output$simulation <- simulation
    return(output)
}
