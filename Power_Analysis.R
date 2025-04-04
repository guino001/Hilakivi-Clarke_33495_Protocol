# Load required packages
install.packages("pwr")

library(pwr)
library(ggplot2)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# --- Power Analysis for Changes in the Gut Microbiome (Human Participants) ---

# Parameters
alpha_microbiome <- 0.01  # Bonferroni-corrected alpha (0.05/5)
power_target <- 0.80     # Target power
n_microbiome <- 20       # Number of participants in crossover design
effect_sizes_microbiome <- seq(0.2, 1.0, by = 0.1)  # Range of effect sizes to explore

# Power calculation for paired t-test across effect sizes
power_microbiome <- sapply(effect_sizes_microbiome, function(d) {
  pwr.t.test(n = n_microbiome, 
             d = d, 
             sig.level = alpha_microbiome, 
             type = "paired", 
             alternative = "two.sided")$power
})

# Create data frame for plotting
df_microbiome <- data.frame(
  EffectSize = effect_sizes_microbiome,
  Power = power_microbiome
)

# Plot Power vs Effect Size for Microbiome
plot_microbiome <- ggplot(df_microbiome, aes(x = EffectSize, y = Power)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "green") +
  labs(title = "Power vs Effect Size for Gut Microbiome Changes",
       subtitle = "n = 20 participants, α = 0.01 (Bonferroni-corrected), Paired t-test",
       x = "Effect Size (Cohen's d)",
       y = "Power") +
  annotate("text", x = 0.5, y = 0.85, label = "d = 0.5\nPower ≈ 0.85", color = "green") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(0.2, 1.0, by = 0.2))

# Save the plot
ggsave("Power_Microbiome.png", plot_microbiome, width = 8, height = 6, dpi = 300)

# Verify sample size for d = 0.5
sample_size_microbiome <- pwr.t.test(d = 0.5, 
                                     sig.level = alpha_microbiome, 
                                     power = power_target, 
                                     type = "paired", 
                                     alternative = "two.sided")$n
# cat("Required sample size for microbiome (d = 0.5, power = 0.80, α = 0.01):", ceiling(sample_size_microbiome), "\n")

# --- Power Analysis for Breast Cancer Risk (Mouse Model) ---

# Parameters
alpha_mouse <- 0.0125  # Bonferroni-corrected alpha (0.05/4)
power_target <- 0.80   # Target power
n_mouse <- 100         # Number of mice per group
effect_sizes_mouse <- seq(0.2, 1.2, by = 0.1)  # Range of effect sizes to explore

# Power calculation for two-sample t-test across effect sizes
power_mouse <- sapply(effect_sizes_mouse, function(d) {
  pwr.t.test(n = n_mouse, 
             d = d, 
             sig.level = alpha_mouse, 
             type = "two.sample", 
             alternative = "two.sided")$power
})

# Create data frame for plotting
df_mouse <- data.frame(
  EffectSize = effect_sizes_mouse,
  Power = power_mouse
)

# Plot Power vs Effect Size for Mouse Model
plot_mouse <- ggplot(df_mouse, aes(x = EffectSize, y = Power)) +
  geom_line(color = "purple", size = 1) +
  geom_point(color = "purple", size = 2) +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.65, linetype = "dashed", color = "green") +
  labs(title = "Power vs Effect Size for Breast Cancer Risk in Mice",
       subtitle = "n = 100 mice per group, α = 0.0125 (Bonferroni-corrected), Two-sample t-test",
       x = "Effect Size (Cohen's d)",
       y = "Power") +
  annotate("text", x = 0.65, y = 0.85, label = "d = 0.65\nPower ≈ 0.95", color = "green") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(breaks = seq(0.2, 1.2, by = 0.2))

# Save the plot
ggsave("Power_Mouse.png", plot_mouse, width = 8, height = 6, dpi = 300)

# Verify sample size for d = 0.65
sample_size_mouse <- pwr.t.test(d = 0.65, 
                                sig.level = alpha_mouse, 
                                power = power_target, 
                                type = "two.sample", 
                                alternative = "two.sided")$n
# cat("Required sample size for mouse model (d = 0.65, power = 0.80, α = 0.0125):", ceiling(sample_size_mouse), "\n")

# --- Combined Operating Characteristics Table ---
operating_characteristics <- data.frame(
  Endpoint = rep(c("Microbiome", "Breast Cancer Risk"), each = 5),
  EffectSize = c(0.3, 0.4, 0.5, 0.65, 0.8, 0.4, 0.5, 0.65, 0.8, 1.0),
  Power = c(
    pwr.t.test(n = n_microbiome, d = 0.3, sig.level = alpha_microbiome, type = "paired")$power,
    pwr.t.test(n = n_microbiome, d = 0.4, sig.level = alpha_microbiome, type = "paired")$power,
    pwr.t.test(n = n_microbiome, d = 0.5, sig.level = alpha_microbiome, type = "paired")$power,
    pwr.t.test(n = n_microbiome, d = 0.65, sig.level = alpha_microbiome, type = "paired")$power,
    pwr.t.test(n = n_microbiome, d = 0.8, sig.level = alpha_microbiome, type = "paired")$power,
    pwr.t.test(n = n_mouse, d = 0.4, sig.level = alpha_mouse, type = "two.sample")$power,
    pwr.t.test(n = n_mouse, d = 0.5, sig.level = alpha_mouse, type = "two.sample")$power,
    pwr.t.test(n = n_mouse, d = 0.65, sig.level = alpha_mouse, type = "two.sample")$power,
    pwr.t.test(n = n_mouse, d = 0.8, sig.level = alpha_mouse, type = "two.sample")$power,
    pwr.t.test(n = n_mouse, d = 1.0, sig.level = alpha_mouse, type = "two.sample")$power
  )
)

# Round power values for readability
operating_characteristics$Power <- round(operating_characteristics$Power, 2)

# Print table
print("Operating Characteristics Table:")
print(operating_characteristics)

# Save table as CSV for attachment
# write.csv(operating_characteristics, "Operating_Characteristics.csv", row.names = FALSE)

# Display plots
print(plot_microbiome)
print(plot_mouse)