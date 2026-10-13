simulate_SIR <- function(edges, beta = 0.3, gamma = 0.1, nsteps = 50,
                         initial_infected = 1, verbose = TRUE) {
  # Ensure data.table
  edges <- as.data.table(edges)
  setnames(edges, c("Var1", "Var2"))

  # Build undirected adjacency list
  adj <- rbind(
    edges[, .(node = Var1, neighbor = Var2)],
    edges[, .(node = Var2, neighbor = Var1)]
  )[, .(neighbors = list(unique(neighbor))), by = node]

  nodes <- adj$node

  # Initialize states
  state <- data.table(node = nodes, status = "S")
  state[sample(.N, initial_infected), status := "I"]

  # Record time series
  ts <- data.table(step = 0,
                   S = sum(state$status == "S"),
                   I = sum(state$status == "I"),
                   R = sum(state$status == "R"))

  # Simulation loop
  for (t in 1:nsteps) {
    infected <- state[status == "I", node]
    susceptible <- state[status == "S", node]

    new_infected <- character(0)
    for (i in infected) {
      neigh <- adj[node == i, neighbors][[1]]
      sus_neigh <- neigh[neigh %in% susceptible]
      newly <- sus_neigh[runif(length(sus_neigh)) < beta]
      new_infected <- c(new_infected, newly)
    }

    recovered <- infected[runif(length(infected)) < gamma]

    state[node %in% new_infected, status := "I"]
    state[node %in% recovered, status := "R"]

    ts <- rbind(ts, data.table(
      step = t,
      S = sum(state$status == "S"),
      I = sum(state$status == "I"),
      R = sum(state$status == "R")
    ))

    if (verbose)
      cat(sprintf("Step %d: S=%d I=%d R=%d\n",
                  t, ts$S[t+1], ts$I[t+1], ts$R[t+1]))

    if (sum(state$status == "I") == 0) break
  }

  # --- Summary measures ---
  peak_I <- max(ts$I)
  time_to_peak <- ts$I %>% which.max() - 1
  final_size <- max(ts$R)
  duration <- max(ts$step)
  total_nodes <- nrow(state)

  summary <- data.table(
    peak_I = peak_I,
    peak_I_frac = peak_I / total_nodes,
    time_to_peak = time_to_peak,
    final_size = final_size,
    final_size_frac = final_size / total_nodes,
    duration = duration
  )

  return(list(timeseries = ts, summary = summary))
}

#>    peak_I peak_I_frac time_to_peak final_size final_size_frac duration
#> 1:      3        0.50            4          6              1        8

head(res$timeseries)
