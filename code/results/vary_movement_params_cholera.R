# required packages
library(ggplot2)
library(deSolve)
library(reshape2)
library(tidyverse)
library(cowplot)
library(pracma)
library(RColorBrewer)
library(scales)

# to paralellize
library(doParallel)
library(foreach)
registerDoParallel(detectCores() - 2) # update this if you want to use more cores

# load optimal control files
source("code/optimal_control_functions.R")

chol_mod_details <- setup_model("cholera")

#### FUNCTIONS -----------------------------------------------------------------
sens_analysis_setup <- function(change_params, 
                                multiplier, 
                                id_col = NA){
  #browser()
  # create data.frame with parameters to change (start with baseline)
  test_params <- data.frame(change_params)
  # expand to have one row per parameter to change
  test_params <- test_params[rep(1, times = length(change_params) + 1),]
  # change each parameter by x-fold once
  for(i in 1:length(change_params)){
    test_params[i,i] = multiplier*test_params[i,i]
  }
  # create id column if desired
  if(!any(is.na(id_col))){
    test_params$id = id_col
  }
  # add rows for unique and uniform
  test_params <- bind_rows(
    test_params %>% mutate(control_type = "unique"),
    test_params %>% mutate(control_type = "uniform")
  )
}

#### RUN COST CHANGE OPTIMAL CONTROL CASES ----------------------------------------

# setup test parameters
test_params <- sens_analysis_setup(change_params = chol_mod_details$params[c("m1", "m2")],
                                   multiplier = 10,
                                   id_col = c("m1", "m2", "base"))
# add third case where we increase both by ten fold
test_params <- test_params %>% 
  bind_rows(data.frame(m1 = rep(unlist(chol_mod_details$params["m1"])*10,2), 
                       m2 = rep(unlist(chol_mod_details$params["m2"])*10,2), 
                       id = rep("m1_m2",2), 
                       control_type = c("unique", "uniform"))) %>%
  arrange(control_type, id)


# iterate over each parameter set in test_params
test2 <- foreach(
  i = 1:nrow(test_params),
  .packages = c("deSolve", "tidyverse", "pracma")
) %dopar% {
  oc_optim(
    model = "cholera",
    change_params = test_params[i,]
  )
}
test_params$test_case <- 1:nrow(test_params)

# print number of iterations for each run
sapply(test2, function(i){i$n_iterations})

#### MANIPULATE OUTPUT FOR PLOTTING --------------------------------------------
states <- lapply(
  1:nrow(test_params),
  function(i) { # browser();
    return(data.frame(
      test_case = i,
      reshape2::melt(test2[[i]]$trajectories, "time")
    ))
  }
)
states <- as.data.frame(do.call(rbind, states))
states <- left_join(test_params, states)

# reformat for easy plotting
states <- states %>%
  mutate(
    patch = substr(variable, 2, 2),
    variable = substr(variable, 1, 1)
  ) %>%
  mutate(
    plot_var = ifelse(control_type == "unique", paste("patch", patch), "uniform")
  )



# relative costs
j_vals <- lapply(
  1:nrow(test_params),
  function(i) {
    j <- test2[[i]][["j"]]
    j <- c(j, j_tot = sum(j[substr(names(j), 1, 1) == "j"]))
    return(data.frame(test_case = i, melt(j)))
  }
)
j_vals <- as.data.frame(do.call(rbind, j_vals))
j_vals <- left_join(test_params, j_vals, by = "test_case") %>%
  rename(variable = L1)



