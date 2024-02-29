compute_SFS_one_iteration <- function(vec_A_and_K, vec_p, vec_SFS_positions, library_SFS_component) {
    # 	Add the neutral component
    A <- vec_A_and_K[1]
    vec_SFS_model <- A * unlist(library_SFS_component[[1, 1]])

    # print("...")
    # print(sum(unlist(library_SFS_component[[1, 1]])))

    # 	Add the binomial humps
    for (i_hump in seq_along(vec_p)) {
        p <- vec_p[i_hump]
        loc <- which(vec_SFS_positions == p)
        K <- vec_A_and_K[i_hump + 1]
        vec_SFS_model <- vec_SFS_model + K * unlist(library_SFS_component[[2, loc]])
    }
    # 	Return the full SFS
    return(vec_SFS_model)
}
