# Load required libraries
library(tidyverse)
library(pwr)
library(ggplot2)
library(cowplot)

# ----------------------------
# 1. Power Analysis: Microbiome Shifts (Paired t-test)
# ----------------------------
alpha <- 0.05
power_target <- 0.80
effect_size <- 0.7  # Cohen's d
n_microbiome <- pwr.t.test(
  d = effect_size,
  sig.level = alpha,
  power = power_target,
  type = "paired"
)$n

# Adjust for repeated measures (ICC = 0.3)
n_adjusted <- ceiling(n_microbiome / (1 - 0.3))  # Final n=20

# ----------------------------
# 2. Power Analysis: FMT Tumor Reduction (Mixed Model)
# ----------------------------
set.seed(123)
n_mice_per_group <- 20
effect_tumor <- 0.4  # 40% reduction
sim_power <- mean(replicate(1000, {
  fmt_effect <- rnorm(n_mice_per_group, mean = 1 - effect_tumor, sd = 0.25)
  control <- rnorm(n_mice_per_group, mean = 1, sd = 0.25)
  t.test(control, fmt_effect)$p.value < alpha
}))

# ----------------------------
# 3. Power vs. Effect Size Curves (Separate Plot)
# ----------------------------
effect_sizes <- seq(0.3, 1.0, by = 0.1)
power_curve <- map_dfr(effect_sizes, ~{
  data.frame(
    EffectSize = .x,
    Microbiome = pwr.t.test(d = .x, n = n_adjusted, sig.level = alpha, type = "paired")$power,
    Tumor = pwr.t.test(d = .x, n = n_mice_per_group, sig.level = alpha, type = "two.sample")$power
  )
}) %>% pivot_longer(-EffectSize, names_to = "Endpoint", values_to = "Power")

# Plot 1: Power Curves
plot_power <- ggplot(power_curve, aes(x = EffectSize, y = Power, color = Endpoint)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = power_target, linetype = "dashed", color = "red") +
  annotate("text", x = 0.7, y = 0.85, label = "Target Power (80%)", color = "red", size = 5) +
  labs(
    title = "Statistical Power by Effect Size",
    x = "Effect Size (Cohen's d)",
    y = "Power",
    color = "Endpoint"
  ) +
  scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ----------------------------
# 4. Sample Size Justification (Separate Plot)
# ----------------------------
df_sample <- data.frame(
  Endpoint = c("Microbiome", "Tumor FMT"),
  N = c(n_adjusted, n_mice_per_group * 5)  # 5 groups
)

plot_sample <- ggplot(df_sample, aes(x = Endpoint, y = N, fill = Endpoint)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_text(aes(label = N), vjust = -0.5, size = 6, fontface = "bold") +
  labs(
    title = "Sample Size Justification",
    x = "",
    y = "Total Sample Size"
  ) +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e")) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  ylim(0, max(df_sample$N) * 1.1)

# ----------------------------
# Save Plots Separately
# ----------------------------
ggsave("power_curve.png", plot_power, width = 8, height = 6, dpi = 300)
ggsave("sample_size.png", plot_sample, width = 6, height = 6, dpi = 300)

# Print results
cat("=== Key Results ===\n")
cat("1. Microbiome: n =", n_adjusted, "participants (80% power for d = 0.7)\n")
cat("2. Tumor FMT: 100 mice (power =", round(sim_power * 100, 1), "% for 40% reduction)\n")