library(rethinking)
library(tidyverse)
library(brms)
data("chimpanzees")
chimp_dat <- chimpanzees %>%
  mutate(treatment = 1 + prosoc_left + 2 * condition,
         treatmeent = factor(treatment),
         actor = factor(actor))

dat_list <- list(pulled_left = chimp_dat$pulled_left,
                 actor = as.integer(chimp_dat$actor),
                 treatment = as.integer(chimp_dat$treatment))

q11.4 <- quap(alist(
  pulled_left ~ dbinom(1, p),
  logit(p) <- a[actor] + b[treatment],
  a[actor] ~ dnorm(0, 1.5),
  b[treatment] ~ dnorm(0, 0.5)
), data = dat_list)

m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom(1, p),
    logit(p) <- a[actor] + b[treatment],
    a[actor] ~ dnorm(0, 1.5),
    b[treatment] ~ dnorm(0, 0.5)
  ), data = dat_list, chains = 4, log_lik = TRUE
)
b11.4 <- brm(bf(pulled_left ~ a + b,
                a ~ 0 + actor,
                b ~ 0 + treatment,
                nl = TRUE), data = chimp_dat, family = bernoulli,
             prior = c(prior(normal(0, 1.5), class = b, nlpar = a),
                       prior(normal(0, 0.5), class = b, nlpar = b)),
             iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)
q_samp <- extract.samples(q11.4)
q_draws <- bind_cols(
  q_samp$a %>% 
    as_tibble(.name_repair = ~paste0("b_a_actor", 1:ncol(q_samp$a))) %>% 
    slice_sample(n = 8000) %>% 
    rowid_to_column(var = ".draw"),
  q_samp$b %>% 
    as_tibble(.name_repair = ~paste0("b_b_treatment", 1:ncol(q_samp$b))) %>% 
    slice_sample(n = 8000)
) %>% 
  pivot_longer(-.draw, names_to = "parameter", values_to = "QUAP")
b_draws <- as_draws_df(b11.4) %>% 
  as_tibble() %>% 
  select(-lp__) %>% 
  pivot_longer(cols = -c(.chain, .iteration, .draw),
               names_to = "parameter", values_to = "MCMC")

post_comp <- full_join(b_draws, q_draws, by = c(".draw", "parameter")) %>% 
  pivot_longer(cols = c(MCMC, QUAP), names_to = "type") %>% 
  mutate(parameter = str_replace_all(parameter, "b_[a|b]_([a-z]*)([0-9])",
                                     "\\1 \\2"),
         parameter = str_to_title(parameter))

post_comp %>% 
  ggplot(aes(x = value, color = type)) +
  facet_wrap(~parameter, nrow = 3) +
  geom_density(key_glyph = "timeseries") +
  scale_color_okabeito() +
  labs(x = "Value", y = "Density", color = NULL)
post_comp %>% 
  filter(parameter == "Actor 2") %>% 
  ggplot(aes(x = value, color = type)) +
  geom_density(key_glyph = "timeseries") +
  scale_color_okabeito() +
  labs(x = "Actor 2", y = "Density", color = NULL)

data("Kline")

kline_dat <- Kline %>%
  mutate(P = standardize(log(population)))

no_hawaii <- filter(kline_dat, culture != "Hawaii")
no_hawaii$contact_id <- ifelse(no_hawaii$contact == "high", 2, 1)


dat <- list(
  T = no_hawaii$total_tools,
  P = no_hawaii$P,
  cid = no_hawaii$contact_id
)

m11.10b <- ulam(
  alist(T ~ dpois(lambda),
        log(lambda) <- a[cid] + b[cid] * P,
        a[cid] ~ dnorm(3, 0.5),
        b[cid] ~ dnorm(0, 0.2)
  ), data = dat, chains = 4, log_lik = TRUE
)
precis(m11.10b)

data(eagles)
eagle_dat <- eagles %>%
  as_tibble() %>%
  mutate(pirateL = ifelse(P == "L", 1, 0),
         victimL = ifelse(V == "L", 1, 0),
         pirateA = ifelse(A == "A", 1, 0))
eagle_quap <- quap(alist(y ~ dbinom(n, p),
                         logit(p) <- a + bP * pirateL + bV * victimL + bA * pirateA,
                         a ~ dnorm(0, 1.5),
                         bP ~ dnorm(0, 1),
                         bV ~ dnorm(0, 1),
                         bA ~ dnorm(0, 1)),
                   data = eagle_dat)
dat <- list(pirateL = eagle_dat$pirateL, victimL = eagle_dat$victimL, pirateA = eagle_dat$pirateA, y = eagle_dat$y, n = eagle_dat$n)
eagle_ulam <- ulam(
  alist(y ~ dbinom(n, p),
        logit(p) <- a + bP * pirateL + bV * victimL + bA * pirateA,
        a ~ dnorm(0, 1.5),
        bP ~ dnorm(0, 1),
        bV ~ dnorm(0, 1),
        bA ~ dnorm(0, 1)
        ), data = dat
  )
precis(eagle_quap)
precis(eagle_ulam)

eagle_brms <- brm(y | trials(n) ~ 1 + pirateL + victimL + pirateA,
                  data = eagle_dat, family = binomial,
                  prior = c(prior(normal(0, 1.5), class = Intercept),
                            prior(normal(0, 1), class = b)),
                  iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234
                  )

eagle_dat %>% 
  add_linpred_draws(eagle_brms) %>% 
  mutate(prob = inv_logit_scaled(.linpred),
         label = paste0(P, A, V)) %>% 
  ggplot(aes(x = label, y = prob)) +
  stat_pointinterval(aes(color = "Posterior"), .width = 0.89, size = 5) +
  geom_point(data = eagle_dat, size = 2,
             aes(x = paste0(P, A, V), y = y / n, color = "Observed")) +
  scale_color_manual(values = c("Posterior" = "#009FB7",
                                "Observed" = "#272727"),
                     name = NULL) +
  labs(x = "Case", y = "Probability")

eagle_dat %>% 
  add_epred_draws(eagle_brms) %>% 
  mutate(label = paste0(P, A, V)) %>% 
  ggplot(aes(x = label, y = .epred)) +
  stat_pointinterval(aes(color = "Posterior"), .width = 0.89, size = 5) +
  geom_point(data = eagle_dat, size = 2,
             aes(x = paste0(P, A, V), y = y, color = "Observed")) +
  scale_color_manual(values = c("Posterior" = "#009FB7",
                                "Observed" = "#272727"),
                     name = NULL) +
  labs(x = "Case", y = "Successes")

data("salamanders")
salamander_dat <- salamanders %>%
  mutate(cov = standardize(PCTCOVER),
         age = standardize(FORESTAGE))

sal_quap <- quap(alist(SALAMAN ~ dpois(lambda),
                       log(lambda) <- a + bC * cov,
                       a ~ dnorm(0, 1),
                       bC ~ dnorm(0, 0.5)),
                 data = salamander_dat)
sal_brms <- brm(SALAMAN ~ 1 + cov, data = salamander_dat, family = poisson,
                prior = c(prior(normal(0, 1), class = Intercept),
                          prior(normal(0, 0.5), class = b)),
                iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)
precis(sal_quap)
fixef(sal_brms, probs = c(0.055, 0.945))
set.seed(219)

epreds <- tibble(cov = seq(-2, 1.5, by = 0.05)) %>% 
  add_epred_draws(sal_brms)

preds <- tibble(cov = seq(-2, 1.5, by = 0.01)) %>% 
  add_predicted_draws(sal_brms)

ggplot() +
  stat_lineribbon(data = preds, aes(x = cov, y = .prediction,
                                    fill = "Prediction (89%)"),
                  .width = c(0.89), size = 0) +
  stat_lineribbon(data = epreds, aes(x = cov, y = .epred),
                  .width = c(0.67, 0.89, 0.97), size = 0.5) +
  geom_point(data = salamander_dat, aes(x = cov, y = SALAMAN)) +
  scale_fill_manual(values = c("#F0F0F0", ramp_blue(seq(1, 0.2, length.out = 3))),
                    breaks = c("Prediction (89%)", 0.67, 0.89, 0.97)) +
  labs(x = "Ground cover (standardized)", y = "Observed Salamanders",
       fill = "Interval")
sal_brms2 <- brm(SALAMAN ~ 1 + cov + age, data = salamander_dat, family = poisson,
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.5), class = b)),
                 iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)

summary(sal_brms2)
sal_brms3 <- brm(SALAMAN ~ 1 + age, data = salamander_dat, family = poisson,
                 prior = c(prior(normal(0, 1), class = Intercept),
                           prior(normal(0, 0.5), class = b)),
                 iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)

summary(sal_brms3)

data("NWOGrants")

nwo_dat <- NWOGrants %>%
  mutate(gender = factor(gender, levels = c("m", "f")))

b11h4_total <- brm(awards | trials(applications) ~ 0 + gender, data = nwo_dat,
                   family = binomial,
                   prior = prior(normal(0, 1.5), class = b),
                   iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)
                   
b11h4_direct <- brm(bf(awards | trials(applications) ~ g + d + i,
                       g ~ 0 + gender,
                       d ~ 0 + discipline,
                       i ~ 0 + gender:discipline,
                       nl = TRUE), data = nwo_dat, family = binomial,
                       prior = c(prior(normal(0, 1.5), nlpar = g),
                                 prior(normal(0, 1.5), nlpar = d),
                                 prior(normal(0, 1.5), nlpar = i)),
                       iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)

as_draws_df(b11h4_total) %>%
  mutate(diff_male = inv_logit_scaled(b_genderm) - inv_logit_scaled(b_genderf)) %>%
  ggplot(aes(x = diff_male)) + 
  stat_halfeye(.width = c(0.67, 0.89, 0.97), fill = "#009FB7") +
  labs(x = "&beta;<sub>M</sub> &minus; &beta;<sub>F</sub>", y = "Density")
 
apps_per_dept <- nwo_dat %>%
  group_by(discipline) %>%
  summarize(applications = sum(applications))

male_dat <- apps_per_dept %>%
  mutate(gender = "m") %>%
  uncount(applications) %>%
  mutate(applications = 1L)

female_dat <- apps_per_dept %>%
  mutate(gender = "f") %>%
  uncount(applications) %>%
  mutate(applications = 1L)

marg_eff <- bind_rows(add_epred_draws(male_dat, b11h4_direct),
                      add_epred_draws(female_dat, b11h4_direct)) %>%
  pivot_wider(names_from = "gender", values_from = ".epred") %>%
  mutate(diff = m - f)

ggplot(marg_eff, aes(x = diff)) +
  stat_halfeye(.width = c(0.67, 0.89, 0.97), fill = "#009FB7") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Difference in Awards (Male - Female)", y = "Density")

n <- 1000
g <- rbernoulli(n, p = 0.5)
s <- rbernoulli(n, p = 0.5)
d <- rbernoulli(n, p = inv_logit_scaled(2 * g - s))
a <- rbernoulli(n, p = inv_logit_scaled(0 * g + d + s - 2))

dat <- tibble(g, d, a) %>% 
  mutate(across(everything(), as.integer),
         across(everything(), as.factor))

mod <- brm(a ~ 1 + d + g, data = dat, family = bernoulli,
           prior = c(prior(normal(0, 1), class = Intercept),
                     prior(normal(0, 1), class = b)),
           iter = 4000, warmup = 2000, chains = 4, cores = 4)

as_draws_df(mod, variable = "b_g1") %>% 
  mean_hdi(b_g1, .width = 0.89)

data("Primates301")

primate_dat <- Primates301 %>%
  as_tibble() %>%
  select(social_learning, genus, species, brain, research_effort) %>%
  drop_na(everything()) %>%
  mutate(log_brain = standardize(log(brain)),
         log_effort = log(research_effort)) %>%
  rowid_to_column()

b11h6a <- brm(social_learning ~ 1 + log_brain, data = primate_dat,
              family = poisson,
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 0.5), class = b)),
              iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)

summary(b11h6a)

preds <- primate_dat %>%
  add_predicted_draws(b11h6a)

preds %>%
  filter(rowid %in% 1:50) %>%
  ggplot(aes(x = rowid)) +
  stat_pointinterval(aes(y = .prediction, color = "Posterior"), .width=0.89) +
  geom_point(data = filter(primate_dat, rowid %in% 1:50),
             aes(y = social_learning, color = "Observed")) +
  scale_color_manual(values = c("Posterior" = "#272727", "Observed" = "#009F87")) +
  expand_limits(y = c(0, 200)) +
  labs(x = "Cases", y = "Social learning", color = NULL)

preds %>% 
  filter(rowid %in% 51:100) %>% 
  ggplot(aes(x = rowid)) +
  stat_pointinterval(aes(y = .prediction, color = "Posterior"), .width = 0.89) +
  geom_point(data = filter(primate_dat, rowid %in% 51:100),
             aes(y = social_learning, color = "Observed")) +
  scale_color_manual(values = c("Posterior" = "#272727", "Observed" = "#009FB7")) +
  expand_limits(y = c(0, 200)) +
  labs(x = "Case", y = "Social Learning", color = NULL)

preds %>% 
  filter(rowid %in% 101:150) %>% 
  ggplot(aes(x = rowid)) +
  stat_pointinterval(aes(y = .prediction, color = "Posterior"), .width = 0.89) +
  geom_point(data = filter(primate_dat, rowid %in% 101:150),
             aes(y = social_learning, color = "Observed")) +
  scale_color_manual(values = c("Posterior" = "#272727", "Observed" = "#009FB7")) +
  expand_limits(y = c(0, 200)) +
  labs(x = "Case", y = "Social Learning", color = NULL)

b11h6b <- brm(social_learning ~ 1 + log_brain + log_effort, data = primate_dat,
              family = poisson,
              prior = c(prior(normal(0, 1), class = Intercept),
                        prior(normal(0, 0.5), class = b)),
              iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234,
              )

summary(b11h6b)

b11h6a <- add_criterion(b11h6a, criterion = "loo")
b11h6b <- add_criterion(b11h6a, criterion = "loo")

set.seed(220)

library(gghighlight)
bind_cols(
  primate_dat, 
  as_tibble(b11h6a$criteria$loo$pointwise) %>%
    select(loo1 = elpd_loo),
  as_tibble(b11h6b$criteria$loo$pointwise) %>%
    select(loo2 = elpd_loo)
) %>%
  mutate(diff = loo2 - loo1) %>%
  ggplot(aes(x = diff, y = log_effort)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point() +
  gghighlight(n = 1, diff > 15, label_key = genus, max_highlight = 10) +
  labs(x = "LOO<sub>2</sub> - LOO<sub>1</sub>", y = "Research Effort(log")

nwo_dat %>% 
  group_by(discipline) %>% 
  summarize(f = sum(applications[which(gender == "f")]),
            m = sum(applications[which(gender == "m")]),
            total_apps = sum(applications),
            total_awards = sum(awards)) %>% 
  mutate(female_pct = f / sum(f),
         male_pct = m / sum(m),
         award_pct = total_awards / total_apps) %>% 
  ggplot(aes(x = female_pct, y = male_pct)) +
  geom_point(aes(size = award_pct, color = abs(female_pct - male_pct) > 0.05)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_label_repel(data = ~filter(.x, abs(female_pct - male_pct) > 0.05),
                   aes(label = discipline),
                   max.overlaps = Inf, nudge_y = 0.03) +
  scale_size(breaks = seq(0.1, 0.3, by = 0.04)) +
  scale_color_manual(values = c("black", "#009FB7")) +
  guides(color = "none") +
  expand_limits(x = c(0, 0.35), y = c(0, 0.35)) +
  coord_equal() +
  labs(x = "Female Application Rate", y = "Male Application Rate",
       size = "Award Rate") +
  theme(legend.position = "right")