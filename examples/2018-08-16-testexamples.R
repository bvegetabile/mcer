test_tm <- list('fo_low' = matrix(c(0.01,0.99,0.98,0.02),2,2,T),
                'fo_med' = matrix(c(0.2,0.8,0.9,0.1),2,2,T),
                'fo_hi' = matrix(c(0.45,0.55,0.55,0.45),2,2,T))

chain_length <- 5000
tm <- test_tm[[1]]
sm <- calc_eigen_stat(tm)
CalcMarkovEntropyRate(
    tm,
    sm
)
chain <- sim_mc(tm, chain_length)
chain <- ifelse(chain==1, 'A', 'B')
# chain <- ifelse(rbinom(chain_length, 1, 0.25) == 1, 'C', chain)
test_er <- mcer::efficient_mc_er(chain, 2)
CalcMarkovEntropyRate(
    tm,
    sm
)
test_er
fixedwindow_lz77(chain, chain_length/2)
fixedquick_lz77(chain, chain_length/2)
# mcer::efficient_mc_er(stationary_bootstrap(chain, p = .01),2)

mcer::efficient_mc_er(chain, 1)$log_likelihood
mcer::efficient_mc_er(chain, 2)$log_likelihood
mcer::efficient_mc_er(chain, 3)$log_likelihood
mcer::efficient_mc_er(chain, 4)$log_likelihood
mcer::efficient_mc_er(chain, 5)$log_likelihood
mcer::efficient_mc_er(chain, 6)$log_likelihood


morder_transmat <- function(mc_order = 1, unique_states = 1:2){
    # mc_order <- 3
    # unique_states <- 1:8
    n_states <- length(unique_states)
    state_space <- data.frame(matrix(NA, nrow=n_states, ncol=mc_order))
    for(i in 1:mc_order){
        state_space[,i] <- unique_states
    }
    state_space <- c(expand.grid(state_space), sep="")
    state_space <- do.call(paste, state_space)
    trans_table <- matrix(NA, nrow=length(state_space), ncol=n_states)
    for(i in 1:length(state_space)){
        alphas <- rgamma(n_states, 500, 1000)
        p <- dirmult::rdirichlet(n=1, alphas)
        trans_table[i, ] <- as.vector(p)
    }

    trans_mat <- matrix(0, length(state_space), length(state_space))
    rownames(trans_mat) <- colnames(trans_mat) <- state_space
    for(i in 1:length(state_space)){
        for(j in 1:n_states){
            from_state <- state_space[i]
            to_state <- paste(state_space[i], unique_states[j], sep='')
            to_state <- substr(to_state, 2, 6)
            fs <- which(state_space == from_state)
            ts <- which(state_space == to_state)
            trans_mat[fs, ts] <- trans_table[i, j]
        }
    }
    list('trans_table' = trans_table,
         'trans_mat' = trans_mat,
         'true_er' = CalcMarkovEntropyRate(trans_mat, CalcEigenStationary(trans_mat)))
}


morder_seq <- function(mc_order, unique_states, trans_table, seq_len = 1000){
    n_states <- length(unique_states)
    state_space <- data.frame(matrix(NA, nrow=n_states, ncol=mc_order))
    for(i in 1:mc_order){
        state_space[,i] <- unique_states
    }
    state_space <- c(expand.grid(state_space), sep="")
    state_space <- do.call(paste, state_space)

    samp_seq <- rep(NA, seq_len)
    x0 <- sample(state_space, size = 1, replace = T)
    samp_seq[1:mc_order] <- as.integer(strsplit(x0, split = '')[[1]])
    start_spot <- mc_order + 1
    for(s in start_spot:seq_len){
        left_word <- paste(samp_seq[(s-mc_order):(s-1)], collapse='')
        tp <- trans_table[which(state_space == left_word), ]
        new_symbol <- sample(unique_states, size=1, replace = T, prob = tp)
        samp_seq[s] <- new_symbol
    }
    samp_seq
}

ord <- 4
testmat <- morder_transmat(ord, 1:ord)
seq4 <- morder_seq(ord, 1:ord, testmat$trans_table, 1000)

hmph <- c(
mcer::efficient_mc_er(seq4, 1)$log_likelihood,
mcer::efficient_mc_er(seq4, 2)$log_likelihood,
mcer::efficient_mc_er(seq4, 3)$log_likelihood,
mcer::efficient_mc_er(seq4, 4)$log_likelihood,
mcer::efficient_mc_er(seq4, 5)$log_likelihood,
mcer::efficient_mc_er(seq4, 6)$log_likelihood,
mcer::efficient_mc_er(seq4, 7)$log_likelihood,
mcer::efficient_mc_er(seq4, 8)$log_likelihood,
mcer::efficient_mc_er(seq4, 9)$log_likelihood
)
plot(hmph)
diff(hmph)

S <- seq(-1,1000,.1)
plot(S, pchisq(S, 108))
