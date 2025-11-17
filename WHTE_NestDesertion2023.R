# ---
#   title: "Nest Desertion Script (Following Woodcock codes example-Stephen Oppel) "
# author: "Mariem Belkadhi"
# date: "2025-10-08"
# output: pdf_document
# ---
  
  
  

library(readxl)

obs <- readxl::read_xlsx("Observations2023.xlsx") |> as.data.frame()
head(obs)



# Whiskered Tern brood desertion / departure model (NIMBLE) 
# 3 latent states: 1=two parents present, 2=one parent present, 3=zero parents
# Observations O are 1,2,3 (observed: 2,1,0 parents) with small misclassification


## Libraries
if(!requireNamespace("nimble", quietly=TRUE)) install.packages("nimble")
pkgs <- c("lubridate","tidyverse","data.table","MCMCvis","coda","ggplot2")
for(p in pkgs) if(!requireNamespace(p, quietly=TRUE)) install.packages(p)
library(lubridate)
library(tidyverse)
library(data.table)
library(MCMCvis)
library(coda)
library(nimble)

set.seed(2025)


## Input checks



if("NEST ID" %in% names(obs)) obs$NestID <- as.character(obs[["NEST ID"]])
if("NEST_ID" %in% names(obs)) obs$NestID <- as.character(obs[["NEST_ID"]])
if(!("NestID" %in% names(obs))) stop("Column 'NEST ID' (or 'NEST_ID') not found in 'obs'.")
if(!("Parents" %in% names(obs))) stop("Column 'Parents' not found in 'obs'.")
if(!("DATE" %in% names(obs))) stop("Column 'DATE' not found in 'obs'.")



## ----- Build encounter matrix y (nind x nweeks) --------------------------

obs <- obs %>% mutate(DATE = as.Date(DATE))

survey_dates <- sort(unique(obs$DATE))
nweeks <- length(survey_dates)
nest_ids <- sort(unique(obs$NestID))
nind <- length(nest_ids)

wide <- obs %>% 
  select(NestID, Parents, DATE) %>%
  distinct(NestID, DATE, .keep_all = TRUE) %>%
  mutate(DATE = as.Date(DATE)) %>%
  arrange(NestID, DATE) %>%
  tidyr::pivot_wider(names_from = DATE, values_from = Parents)


# Map Parents -> model categories:
# Parents==2 -> 1  ; Parents==1 -> 2 ; Parents==0 -> 3 ; NA -> NA (missing survey)


y_mat <- matrix(NA_integer_, nrow = nind, ncol = nweeks,
                dimnames = list(nest_ids, as.character(survey_dates)))
for(i in seq_len(nind)){
  row <- wide[i, -1]
  for(j in seq_len(nweeks)){
    v <- suppressWarnings(as.numeric(row[[j]]))
    if(is.na(v)) {
      y_mat[i,j] <- NA_integer_
    } else if(v >= 2) {
      y_mat[i,j] <- 1L
    } else if(v == 1) {
      y_mat[i,j] <- 2L
    } else { # v==0 or other <=0
      y_mat[i,j] <- 3L
    }
  }
}


# First survey index with any observation per nest 

f <- integer(nind)
for(i in seq_len(nind)){
  idx <- which(!is.na(y_mat[i,]))
  f[i] <- if(length(idx)==0) 1L else min(idx)
}

isObs <- 1L * (!is.na(y_mat))
effort <- matrix(1, nrow = nind, ncol = nweeks)
# year   <- rep(1L, nind); nyears <- 1L
# tag    <- rep(0L, nind)


##  Build *monotone* initial z (no return allowed in TRUE states)
# keeping raw y as data (with NAs), but z inits must obey 1→2→3 monotonicity.

z_init <- matrix(NA_integer_, nrow = nind, ncol = nweeks,
                 dimnames = dimnames(y_mat))
mods <- 0L
for(i in seq_len(nind)){
  
  zi <- as.integer(y_mat[i, ])
  
  obs_pos <- which(!is.na(zi))
  if(length(obs_pos) == 0){
    zi[] <- 1L
  } else {
    first <- obs_pos[1]
    
    if(first > 1) zi[1:(first-1)] <- zi[first]  ## you can leave NA before the first observation because the model ignores anything before the first observation
    
    for(t in (first+1):nweeks){
      if(is.na(zi[t])) {
        y_mat[i,t]<-3 ## intermediate NAs should be set to state 3 (no birds observed) in the observation matrix
        zi[t] <- zi[t-1]}
      if(t<nweeks){
        if(min(zi[(t+1):nweeks], na.rm=T)<zi[t]) {zi[t]<-min(zi[(t+1):nweeks], na.rm=T)}  ## the minimum state (max number of ind) is maintained until the last observation of that state
      }
    }
    
    # zi_fix <- cummax(zi)  ## THIS IS WRONG
    # mods <- mods + sum(zi_fix != zi, na.rm = TRUE)
    # zi <- zi_fix
  }
  
  zi[zi < 1] <- 1L; zi[zi > 3] <- 3L
  z_init[i, ] <- zi
}
cat("Initial z modifications to enforce monotonicity:", mods, "\n")


## NIMBLE model: NO RETURN in true states


whi.mig.model <- nimbleCode({
  
  ## Prior parameters (changed to use normal 0 prior on logit scale directly)
  
  ## nest departure is state dependent (for 2 or 1 individual separately)
  for (s in 1:3) {
    lm.mean[s] ~ dnorm(0, 1)
    b.mig.week[s] ~ dnorm(0, 1)  
  }
             
  lp.mean ~ dnorm(0, 0.1)
  b.obs.eff ~ dnorm(1, 1)
  # eps <- 1e-6            ## this is a fixed near-zero observation probability that does not make much sense?                    
  
  
  for (i in 1:nind) {
    logit.mig[i,f[i]] <- lm.mean[1] + b.mig.week[1] * week[f[i]]  ## on first occasion there are 2 parents present (=state 1)
    mig[i,f[i]] <- ilogit(logit.mig[i,f[i]])
    logit.p.obs[i,f[i]] <- lp.mean + b.obs.eff * isObs[i,f[i]]
    p.obs[i,f[i]] <- ilogit(logit.p.obs[i,f[i]])
    for (t in (f[i]+1):nweeks) {
      logit.mig[i,t] <- lm.mean[z[i,t-1]] + b.mig.week[z[i,t-1]] * week[t]  ## state dependent mig probability
      mig[i,t] <- ilogit(logit.mig[i,t])
      logit.p.obs[i,t] <- lp.mean + b.obs.eff * isObs[i,t]
      p.obs[i,t] <- ilogit(logit.p.obs[i,t])
    }
  }
  
  # Transitions & observations (no return in TRUE states)
  for (i in 1:nind) {
    for (t in f[i]:(nweeks-1)) {
      # ---- TRUE state transitions (rows=from, cols=to) ----
      # From 1 (two parents): each parent leaves independently with prob mig
      ps[1,i,t,1] <- (1 - mig[i,t]) * (1 - mig[i,t])            #  2 present
      ps[1,i,t,2] <- mig[i,t] * (1 - mig[i,t])              #  1 present ## SO deleted a 2 *
      ps[1,i,t,3] <- mig[i,t] * mig[i,t]                        #  0 present
      
      # From 2 (one parent present): that parent may leave; absent parent cannot return
      ps[2,i,t,1] <- 0
      ps[2,i,t,2] <- 1 - mig[i,t]
      ps[2,i,t,3] <- mig[i,t]
      
      # From 3 (zero present): absorbing (no return)
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # ---- OBSERVATION model (detection error allowed) ----
      # True=1 (two present): Binomial(2, p.obs)
      po[1,i,t,1] <- p.obs[i,t] * p.obs[i,t]                    # observe 2
      po[1,i,t,2] <- p.obs[i,t] * (1 - p.obs[i,t])          # observe 1  ## SO deleted a 2*
      po[1,i,t,3] <- (1 - p.obs[i,t]) * (1 - p.obs[i,t])        # observe 0
      
      # True=2 (one present): usually observe 1, allow tiny chance to see 2
      po[2,i,t,1] <- 0  ## this should by definition be 0 - you can only know that it is not 0 if birds are individually recognizable? That info is not in this model though!?
      po[2,i,t,2] <- p.obs[i,t]
      po[2,i,t,3] <- 1 - p.obs[i,t]
      
      # True=3 (zero present): mostly 0; allow tiny misclassification to 1 or 2
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1  ## if nobody is at the nest then no bird can be observed
    }
  }
  
  # Masked likelihood for missing surveys:
  #   if isObs=1 -> use po; if isObs=0 -> uniform over {1,2,3}  ## this should be in the observation process
  for (i in 1:nind) {
    # Define latent state at first capture
    z[i,f[i]] <- 1 ## on first record the nest should have had 2 parents - otherwise it is pointless to include it for estimation of nest desertion
    for (t in (f[i]+1):nweeks) {
      # State process - transition probability from one state to the next (migration)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:3])                # state process
      # for (k in 1:3) { ## I am not sure what this will actually achieve?
      #   w[i,t,k] <- isObs[i,t] * po[z[i,t], i, t-1, k] + (1 - isObs[i,t]) * (1/3)
      #   w[i,t,k] <- po[z[i,t], i, t-1, k]
      # }
      # Observation process - probability to observe a number of birds given the true state z[i,t]
      y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:3])                              # observation process
    }
  }
  
})



## Data & constants for NIMBLE 


whi.constants <- list(f = f,
                      week = seq_len(nweeks),
                      nind = nind,
                      #nyears = 1L,  ## not needed in this model
                      #year = rep(1L, nind),   ## not needed in this model
                      #effort = effort, ## not yet included - should affect p.obs, but only if effort can vary - so far it does not, replaced by isObs
                      #tag = tag, ## not needed in this model
                      isObs = isObs, ## this is effectively a binary matrix of whether a nest was observed or not?
                      nweeks = nweeks)

whi.data <- list(y = y_mat)

whi.dims <- list(ps = c(3, nind, nweeks, 3),
                 po = c(3, nind, nweeks, 3))


## Inits & parameters 


parameters <- c("lm.mean", "b.mig.week", "lp.mean","b.obs.eff")

smartInit <- list(
  z = z_init,
  lm.mean = c(0.1,0.1,0.1),
  b.mig.week = c(0.1,0.1,0.1),
  b.obs.eff = 1.5,
  lp.mean = 0
)


## MCMC settings


n.iter   <- 20000
n.burnin <- 10000
n.chains <- 3
thin     <- 2


# PRELIMINARY TEST OF NIMBLE MODEL TO IDENTIFY PROBLEMS --------------------
test <- nimbleModel(code = whi.mig.model,
                    constants = whi.constants,
                    data = whi.data,
                    inits = smartInit,
                    calculate=TRUE)

### make sure that none of the logProbs result in NA or -Inf as the model will not converge
test$calculate()  # will sum all log probs - if there is -Inf or NA then something is not properly initialised
test$initializeInfo()

# use test output as starting values or check where the NA comes from
test$logProb_b.mig.week
test$logProb_b.obs.eff
test$logProb_y ### there should not be any -Inf in this matrix
y_mat[2,6:12]
z_init[2,6:12]
isObs[2,6:12]



tic <- function() assign(".timer", proc.time(), envir=.GlobalEnv)
toc <- function() print((proc.time()-get(".timer", envir=.GlobalEnv))["elapsed"])
tic()
whi_mcmc <- nimbleMCMC(code = whi.mig.model,
                       constants = whi.constants,
                       data = whi.data,
                       inits = smartInit,
                       dimensions = whi.dims,
                       monitors = parameters,
                       thin = thin,
                       niter = n.iter,
                       nburnin = n.burnin,
                       nchains = n.chains,
                       summary = TRUE)
toc()

if(!dir.exists("output")) dir.create("output", showWarnings = FALSE)
print(whi_mcmc$summary)


## Posterior: weekly desertion probability 

MCMCout  <- do.call(rbind, whi_mcmc$samples)
parmcols <- colnames(MCMCout)
nSamps   <- nrow(MCMCout)
weeks    <- seq_len(nweeks)

lm_vec <- as.numeric(MCMCout[, "lm.mean[1]"])
bmig_vec    <- as.numeric(MCMCout[, "b.mig.week[1]"])

Wmat          <- matrix(weeks, nrow = nSamps, ncol = nweeks, byrow = TRUE)
logit_mig_mat <- matrix(lm_vec, nSamps, nweeks) + matrix(bmig_vec, nSamps, nweeks) * Wmat
mig_mat       <- plogis(logit_mig_mat)

mig_summary <- data.frame(
  week = weeks,
  mig_median = apply(mig_mat, 2, median),
  mig_lcl    = apply(mig_mat, 2, quantile, 0.025),
  mig_ucl    = apply(mig_mat, 2, quantile, 0.975),
  Date       = survey_dates
)

ggplot(mig_summary, aes(x = Date)) +
  geom_ribbon(aes(ymin = mig_lcl, ymax = mig_ucl), alpha = 0.20, fill = "firebrick") +
  geom_line(aes(y = mig_median), linewidth = 1, color = "firebrick") +
  scale_y_continuous("Weekly prob. first parent leaves", limits = c(0,1)) +
  scale_x_date("Survey date") +
  theme_minimal()



## Simulate cumulative proportion deserted (pop = 1000) 


pop0 <- 1000
left_mat    <- matrix(0L, nrow = nSamps, ncol = nweeks)
pop_remain  <- matrix(pop0, nrow = nSamps, ncol = nweeks)
for(w in 1:nweeks){
  probs  <- mig_mat[, w]
  left_w <- rbinom(nSamps, pop_remain[, w], probs)
  left_mat[, w] <- left_w
  if(w < nweeks) pop_remain[, w+1] <- pop_remain[, w] - left_w
}
prop_mig_mat <- t(apply(left_mat, 1, cumsum)) / pop0
prop_summary <- data.frame(
  week    = weeks,
  prop_med= apply(prop_mig_mat, 2, median),
  prop_lcl= apply(prop_mig_mat, 2, quantile, 0.025),
  prop_ucl= apply(prop_mig_mat, 2, quantile, 0.975),
  Date    = survey_dates
)

ggplot(prop_summary, aes(x = Date)) +
  geom_ribbon(aes(ymin = prop_lcl, ymax = prop_ucl), alpha = 0.20, fill = "firebrick") +
  geom_line(aes(y = prop_med), linewidth = 1, color = "firebrick") +
  scale_y_continuous("Cumulative proportion nests deserted (sim.)", limits = c(0,1)) +
  scale_x_date("Survey date") +
  theme_minimal()





# R-hat and ESS



library(MCMCvis)
MCMCsummary(whi_mcmc$samples, params = c("lm.mean", "b.mig.week", "lp.mean","b.obs.eff"))
MCMCplot(whi_mcmc$samples, params = c("lm.mean", "b.mig.week", "lp.mean","b.obs.eff"))

## ============================================================
## 1) 1→2 (two→one parent) transition probability over time
## ============================================================



# If mig_mat wasn't created above, rebuild it (safe guard)
if(!exists("mig_mat")) {
  weeks  <- seq_len(nweeks)
  lm_vec <- as.numeric(MCMCout[, "lm.mean"])
  bmig_vec <- as.numeric(MCMCout[, "b.mig.week"])
  Wmat <- matrix(weeks, nrow = nrow(MCMCout), ncol = nweeks, byrow = TRUE)
  logit_mig_mat <- matrix(lm_vec, nrow(MCMCout), nweeks) + 
    matrix(bmig_vec, nrow(MCMCout), nweeks) * Wmat
  mig_mat <- plogis(logit_mig_mat) # per-parent leaving prob m_t
}

# ---- Choose formulation for P(1→2) ----
use_independent_parents <- TRUE
# TRUE  -> P12 = 2 * m * (1 - m)      (recommended under independence)
# FALSE -> P12 = 1 * m * (1 - m)      (matches your current ps line exactly)

if(use_independent_parents) {
  P12_mat <- 2 * mig_mat * (1 - mig_mat)
} else {
  P12_mat <-     mig_mat * (1 - mig_mat)
}

P12_summary <- data.frame(
  Date = survey_dates,
  P12_median = apply(P12_mat, 2, median),
  P12_lcl    = apply(P12_mat, 2, quantile, 0.025),
  P12_ucl    = apply(P12_mat, 2, quantile, 0.975)
)

library(ggplot2)
ggplot(P12_summary, aes(x = Date)) +
  geom_ribbon(aes(ymin = P12_lcl, ymax = P12_ucl), alpha = 0.20) +
  geom_line(aes(y = P12_median), linewidth = 1) +
  labs(x = "Survey date",
       y = expression(P[12]~"(P(two" %->% "one) per interval)"),
       title = "Transition probability 1→2 over time") +
  theme_minimal()

## (Optional) Overlay BOTH versions in one plot for comparison:
P12_indep <- 2 * mig_mat * (1 - mig_mat)
P12_model <-     mig_mat * (1 - mig_mat)
cmp <- data.frame(
  Date = rep(survey_dates, 2),
  median = c(apply(P12_indep, 2, median), apply(P12_model, 2, median)),
  lcl    = c(apply(P12_indep, 2, quantile, 0.025), apply(P12_model, 2, quantile, 0.025)),
  ucl    = c(apply(P12_indep, 2, quantile, 0.975), apply(P12_model, 2, quantile, 0.975)),
  version = rep(c("Independent (2*m*(1-m))", "As-coded (m*(1-m))"), each = length(survey_dates))
)

ggplot(cmp, aes(x = Date, y = median, colour = version, linetype = version)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = version), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("#1f77b4", "#d62728")) +
  scale_fill_manual(values   = c("#1f77b4", "#d62728")) +
  labs(x = "Survey date",
       y = "P(1→2)",
       title = "Partial desertion (1→2) under two formulations") +
  theme_minimal() +
  theme(legend.position = "bottom")



## Full transition matrix summaries from posterior (mig_mat)



#  Choose formulation for the 1→2 cell
independent_for_matrix <- TRUE   # TRUE: P12 = 2*m*(1-m)   FALSE: P12 = m*(1-m)

# Helper to compute P-row from m for a *single* time t
.row_from_state1 <- function(m, independent = TRUE) {
  if (independent) {
    c(P11 = (1-m)^2, P12 = 2*m*(1-m), P13 = m^2)
  } else {
    c(P11 = (1-m)^2, P12 =   m*(1-m), P13 = m^2)
  }
}
.row_from_state2 <- function(m) c(P21 = 0*m, P22 = 1-m, P23 = m)
.row_from_state3 <- function(m) c(P31 = 0*m, P32 = 0*m, P33 = 1+0*m)

q3 <- function(x) c(lcl = stats::quantile(x, 0.025), med = stats::median(x), ucl = stats::quantile(x, 0.975))

# Build arrays: P_med, P_lcl, P_ucl with dims [from, to, t]
P_med <- array(NA_real_, dim = c(3,3,nweeks), dimnames = list(from = 1:3, to = 1:3, t = as.character(survey_dates)))
P_lcl <- P_med; P_ucl <- P_med

ps_long <- vector("list", nweeks)

for (tt in seq_len(nweeks)) {
  m_t <- mig_mat[, tt]               # posterior samples of m at week t
  # distributions for each cell
  P11 <- (1 - m_t)^2
  P12 <- if (independent_for_matrix) 2*m_t*(1-m_t) else m_t*(1-m_t)
  P13 <- (m_t)^2
  P21 <- 0*m_t
  P22 <- 1 - m_t
  P23 <- m_t
  P31 <- 0*m_t
  P32 <- 0*m_t
  P33 <- 1 + 0*m_t
  
  # summaries
  s11 <- q3(P11); s12 <- q3(P12); s13 <- q3(P13)
  s21 <- q3(P21); s22 <- q3(P22); s23 <- q3(P23)
  s31 <- q3(P31); s32 <- q3(P32); s33 <- q3(P33)
  

  P_med[1,1,tt] <- s11["med"]; P_lcl[1,1,tt] <- s11["lcl"]; P_ucl[1,1,tt] <- s11["ucl"]
  P_med[1,2,tt] <- s12["med"]; P_lcl[1,2,tt] <- s12["lcl"]; P_ucl[1,2,tt] <- s12["ucl"]
  P_med[1,3,tt] <- s13["med"]; P_lcl[1,3,tt] <- s13["lcl"]; P_ucl[1,3,tt] <- s13["ucl"]
  
  P_med[2,1,tt] <- s21["med"]; P_lcl[2,1,tt] <- s21["lcl"]; P_ucl[2,1,tt] <- s21["ucl"]
  P_med[2,2,tt] <- s22["med"]; P_lcl[2,2,tt] <- s22["lcl"]; P_ucl[2,2,tt] <- s22["ucl"]
  P_med[2,3,tt] <- s23["med"]; P_lcl[2,3,tt] <- s23["lcl"]; P_ucl[2,3,tt] <- s23["ucl"]
  
  P_med[3,1,tt] <- s31["med"]; P_lcl[3,1,tt] <- s31["lcl"]; P_ucl[3,1,tt] <- s31["ucl"]
  P_med[3,2,tt] <- s32["med"]; P_lcl[3,2,tt] <- s32["lcl"]; P_ucl[3,2,tt] <- s32["ucl"]
  P_med[3,3,tt] <- s33["med"]; P_lcl[3,3,tt] <- s33["lcl"]; P_ucl[3,3,tt] <- s33["ucl"]
  
  # long table for plotting/tables
  ps_long[[tt]] <- data.frame(
    Date = survey_dates[tt],
    from = rep(1:3, each = 3),
    to   = rep(1:3, times = 3),
    lcl  = c(s11["lcl"], s12["lcl"], s13["lcl"],
             s21["lcl"], s22["lcl"], s23["lcl"],
             s31["lcl"], s32["lcl"], s33["lcl"]),
    med  = c(s11["med"], s12["med"], s13["med"],
             s21["med"], s22["med"], s23["med"],
             s31["med"], s32["med"], s33["med"]),
    ucl  = c(s11["ucl"], s12["ucl"], s13["ucl"],
             s21["ucl"], s22["ucl"], s23["ucl"],
             s31["ucl"], s32["ucl"], s33["ucl"])
  )
}
ps_long <- dplyr::bind_rows(ps_long)  # requires dplyr (already loaded)


show_P <- function(date = NULL, idx = NULL, stat = c("med","lcl","ucl")) {
  stat <- match.arg(stat)
  if (!is.null(date)) {
    idx <- match(as.Date(date), as.Date(dimnames(P_med)[[3]]))
    if (is.na(idx)) stop("Date not found in survey_dates.")
  }
  if (is.null(idx)) idx <- 1L
  Mat <- switch(stat,
                med = P_med[,,idx],
                lcl = P_lcl[,,idx],
                ucl = P_ucl[,,idx])
  dimnames(Mat) <- list(from = c("1(two)","2(one)","3(zero)"),
                        to   = c("1(two)","2(one)","3(zero)"))
  return(round(Mat, 3))
}

## Examples:
# First date (median)
print(show_P(idx = 1,  stat = "med"))
# Last date (median)
print(show_P(idx = nweeks, stat = "med"))



make_wide <- function(arr, label = "med") {
  # arr: P_med / P_lcl / P_ucl
  DD <- lapply(seq_len(dim(arr)[3]), function(k) {
    M <- arr[,,k]
    data.frame(
      Date = as.Date(dimnames(arr)[[3]][k]),
      P11 = M[1,1], P12 = M[1,2], P13 = M[1,3],
      P21 = M[2,1], P22 = M[2,2], P23 = M[2,3],
      P31 = M[3,1], P32 = M[3,2], P33 = M[3,3]
    )
  })
  out <- dplyr::bind_rows(DD)
  names(out) <- c("Date", paste0(names(out)[-1], "_", label))
  out
}
wide_med <- make_wide(P_med, "med")
wide_lcl <- make_wide(P_lcl, "lcl")
wide_ucl <- make_wide(P_ucl, "ucl")
wide_all <- wide_med %>%
  dplyr::left_join(wide_lcl, by = "Date") %>%
  dplyr::left_join(wide_ucl, by = "Date")



## Plot ALL 9 cells (medians + ribbons) as a 3x3 facet grid
ps_long$cell <- factor(paste0(ps_long$from, "→", ps_long$to),
                       levels = c("1→1","1→2","1→3",
                                  "2→1","2→2","2→3",
                                  "3→1","3→2","3→3"))
ggplot(ps_long, aes(x = Date)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.20) +
  geom_line(aes(y = med), linewidth = 0.9) +
  facet_wrap(~ cell, ncol = 3, scales = "free_y") +
  labs(x = "Date", y = "Probability", title = "Time-varying transition matrix P(t)") +
  theme_minimal()




## --Sanity check & robust rebuild of mig_mat and P(1→2) summaries ---

# Stack posterior draws and sanity-check
MCMCout <- do.call(rbind, whi_mcmc$samples)

if (nrow(MCMCout) < 100) {
  stop("Too few posterior draws. Increase n.iter / reduce burn-in or check monitors.")
}
if (!all(c("lm.mean","b.mig.week") %in% colnames(MCMCout))) {
  stop("Columns 'lm.mean' and/or 'b.mig.week' not found in samples; check monitors.")
}

#  Drop any bad draws (NA/Inf) BEFORE computing
ok <- is.finite(MCMCout[,"lm.mean"]) & is.finite(MCMCout[,"b.mig.week"])
MCMCout <- MCMCout[ok, , drop = FALSE]

# Recompute per-week per-parent leaving probability m_t
weeks <- seq_len(nweeks)
Wmat  <- matrix(weeks, nrow = nrow(MCMCout), ncol = nweeks, byrow = TRUE)
logit_mig_mat <- matrix(MCMCout[,"lm.mean"], nrow(MCMCout), nweeks) +
  matrix(MCMCout[,"b.mig.week"], nrow(MCMCout), nweeks) * Wmat
mig_mat <- plogis(logit_mig_mat)  # nDraws x nweeks

# Compute P(1→2) with your preferred formulation
use_independent_parents <- TRUE  # TRUE: 2*m*(1-m)  | FALSE: m*(1-m) to match your as-coded line
P12_mat <- if (use_independent_parents) 2 * mig_mat * (1 - mig_mat) else mig_mat * (1 - mig_mat)


q3 <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(lcl=NA_real_, med=NA_real_, ucl=NA_real_))
  c(lcl = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    med = stats::median(x, na.rm = TRUE),
    ucl = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE))
}

P12_summary <- t(apply(P12_mat, 2, q3))
P12_summary <- data.frame(
  Date = survey_dates,
  lcl  = P12_summary[, "lcl"],
  med  = P12_summary[, "med"],
  ucl  = P12_summary[, "ucl"]
)


print(head(P12_summary, 10))







