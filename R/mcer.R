
efficient_mc_er <- function(mc_seq, mc_order=1){
    unique_states <- as.character(unique(mc_seq))
    n_states <- length(unique_states)
    state_space <- data.frame(matrix(NA, nrow=n_states, ncol=mc_order))
    for(i in 1:mc_order){
        state_space[,i] <- unique_states
    }
    state_space <- c(expand.grid(state_space), sep="")
    state_space <- do.call(paste, state_space)

    char_seq <- paste(mc_seq, collapse = '')
    # print(unique_states)
    res <- mcer(char_seq, mc_order, unique_states, state_space)
    colnames(res[[2]]) <- unique_states
    rownames(res[[2]]) <- state_space
    colnames(res[[3]]) <- 'Prob'
    rownames(res[[3]]) <- state_space

    res
}

fixedwindow_lz77 <- function(mc_seq, window_size = 10){
    char_seq <- paste(mc_seq, collapse = '')
    lz77entropy(char_seq, window_size)
}

fixedquick_lz77 <- function(mc_seq, window_size = 10){
    char_seq <- paste(mc_seq, collapse = '')
    lz77entropy_quick(char_seq, window_size)
}

block_bootstrap <- function(mc_seq, block_size){
    n_blocks <- floor(length(mc_seq) / block_size)
    new_seq <- rep(NA, length = n_blocks * block_size)
    blocks <- sample(1:n_blocks, replace = T)
    for(i in 1:n_blocks){
        beg_ind <- (i-1)*block_size + 1
        end_ind <- beg_ind + block_size - 1

        block_start <- (blocks[i] - 1) * block_size + 1
        block_end <- block_start + block_size - 1
        new_seq[beg_ind:end_ind] <- mc_seq[block_start:block_end]
    }
    return(new_seq)
}

overlap_bootstrap <- function(mc_seq, block_size){
    seq_len <- length(mc_seq)
    n_blocks <- floor(seq_len / block_size)
    new_seq <- rep(NA, length = n_blocks * block_size)
    block_loc <- sample(1:(seq_len - block_size + 1), replace = T)
    for(i in 1:n_blocks){
        beg_ind <- (i-1)*block_size + 1
        end_ind <- beg_ind + block_size - 1
        block_start <- block_loc[i]
        block_end <- block_start + block_size - 1
        new_seq[beg_ind:end_ind] <- mc_seq[block_start:block_end]
    }
    return(new_seq)
}

stationary_bootstrap <- function(mc_seq, p){
    seq_len <- length(mc_seq)
    block_sizes <- rgeom(seq_len, p) + 1
    block_locs <- cumsum(block_sizes)
    n_blocks <- sum(block_locs <= seq_len) + 1
    block_starts <- block_locs - block_sizes + 1
    block_ind <- sample(1:seq_len, size = n_blocks, replace = T)
    new_seq <- rep(NA, length = cumsum(block_sizes)[n_blocks])
    for(i in 1:n_blocks){
        block_size <- block_sizes[i]
        # beg_ind <- (i-1)*block_size + 1
        beg_ind <- block_starts[i]
        end_ind <- beg_ind + block_size - 1
        block_start <- block_ind[i]
        block_end <- block_start + block_size - 1
        seq_pos <- ((block_start:block_end - 1) %% seq_len) + 1
        new_seq[beg_ind:end_ind] <- mc_seq[seq_pos]
    }
    new_seq[1:seq_len]
}

sim_mc <- function(trans_mat,
                   n_sims=100){
    n_states <- nrow(trans_mat)
    states <- seq(1, n_states)

    simulations <- matrix(0, nrow = 1, ncol = n_sims)

    stat_mat <- calc_eigen_stat(trans_mat = trans_mat)
    init_state <- sample(x = states,
                         size = 1,
                         replace = TRUE,
                         prob = stat_mat)
    simulations[1,1] <- init_state

    for (i in 2:n_sims){
        prev_step <- simulations[1,(i-1)]
        next_step <- sample(states,
                            size = 1,
                            replace = TRUE,
                            prob = trans_mat[prev_step,])
        simulations[1, i] <- next_step
    }
    return(as.vector(simulations))
}


calc_eigen_stat <- function(trans_mat){
    tm_eig <- eigen(t(trans_mat))
    if(any(round(Mod(tm_eig$values),10)==1)){
        lamb1 <- which(abs(tm_eig$values-1) == min(abs(tm_eig$values-1)))
        stat_vec <- tm_eig$vectors[,lamb1] / sum(tm_eig$vectors[,lamb1])
        return(Re(stat_vec))
    } else{
        stat_vec <- rep(0, nrow(trans_mat))
        return(stat_vec)
    }
}
