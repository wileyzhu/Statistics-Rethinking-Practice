# 📘 Statistical Rethinking Practice

This repository contains my personal code practice, notes, and model implementations based on:

> **Statistical Rethinking: A Bayesian Course with Examples in R and Stan**  
> by *Richard McElreath*

It also incorporates practice insights and modeling exercises adapted from:

> 🔗 [**SR2 Solutions**](https://sr2-solutions.wjakethompson.com/) by Jake Thompson — a structured and annotated solution set to the 2nd edition of Statistical Rethinking.

---

## 🎯 Purpose

The main goal of this repository is to:

- Reinforce Bayesian modeling concepts through practical coding
- Explore real and simulated examples of important ideas (e.g. post-treatment bias, multicollinearity, spurious relationships)
- Transition from intuition → modeling → checking assumptions and predictions
- Prepare for more advanced model-building and research use cases

---

## 🗂️ Repository Overview

Instead of strictly following chapters, this repo is organized by **conceptual topics**, each drawn from the book or from SR2-inspired exercises.

| File                      | Focus Area                                          |
|---------------------------|-----------------------------------------------------|
| `Chapter11_Practice.R`    | Practice based on Chapter 11 (e.g. collider bias)   |
| `Monte_Carlos.R`          | Simulation and uncertainty via Monte Carlo          |
| `Masked_Relationship.R`   | Confounding and unobserved variables                |
| `Post_treatment_bias.R`   | Consequences of post-treatment adjustment           |
| `Multicollinearity.R`     | Identifiability and correlated predictors           |
| `Spurious_Waffle.R`       | Spurious correlation (Waffle/Divorce example)       |
| `Conditioning.R`          | Impact of conditioning on different variables       |
| `Categorical_Variables.R` | Encoding and modeling of categorical predictors     |
| `Uleyss_Compass.R`        | Toy example exploring DAG logic and identifiability |

---

## 🧰 Tools Used

- 📦 [`rethinking`](https://github.com/rmcelreath/rethinking): Bayesian modeling (`quap`, `ulam`)
- 📦 `rstan`: MCMC backend for full Bayesian inference
- 📦 `ggplot2`, `dagitty`, and `tidyverse` for visualization and data prep
- 🧪 R script files with inline comments and simulations

---

## 📈 Learning Highlights

- Priors, likelihoods, and Bayesian updating
- Posterior simulations and predictive checks
- Multilevel models (hierarchical structure and partial pooling)
- DAG-based thinking for causality
- Quadratic approximation vs. full MCMC

---

## 🚧 Work in Progress

- [ ] Expand modeling practice to more chapters
- [ ] Add markdown summaries per concept (`.md` or `.Rmd`)
- [ ] Add posterior predictive plots
- [ ] Begin full side project (e.g. **clutch performance in esports**)

---

## 🤝 Acknowledgements

- 📘 Richard McElreath for the Statistical Rethinking book and lectures  
- 🧠 Jake Thompson for the [SR2 Solutions](https://sr2-solutions.wjakethompson.com/) reference site  
- 🙏 The R and Bayesian modeling community

---

> _“Bayesian modeling: making uncertainty visible, one posterior at a time.”_
