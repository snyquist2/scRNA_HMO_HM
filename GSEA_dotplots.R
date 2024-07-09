
library(ggplot2)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

TypeI_LC2 <- read.csv("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/TypeIvsTypeII_final/Secretory_Lactocytes/GOBP_LC2_TypeI.csv")
TypeII_LC2 <- read.csv("/Users/f006f9w/Dropbox (Dartmouth College)/Mac/Desktop/HMO_data/TypeIvsTypeII_final/Secretory_Lactocytes/GOBP_LC2_TypeII.csv")

TypeI_LC2$enrichment_direction <- ifelse(TypeI_LC2$NES > 0, "Positively Enriched", "Negatively Enriched")
TypeII_LC2$enrichment_direction <- ifelse(TypeII_LC2$NES > 0, "Positively Enriched", "Negatively Enriched")

TypeI_LC2_subset <- TypeI_LC2[, c("Description", "NES", "p.adjust")]
TypeI_LC2_subset$Type <- "Type-I"

TypeII_LC2_subset <- TypeII_LC2[, c("Description", "NES", "p.adjust")]
TypeII_LC2_subset$Type <- "Type-II"

merged_df <- rbind(TypeI_LC2_subset, TypeII_LC2_subset)

## Positive Bubble Plot 
positive_df <- merged_df %>% 
  filter(NES > 0) %>%
  arrange(desc(NES))

df_top_40 <- positive_df[1:40, ]
df_remaining <- positive_df[41:nrow(positive_df), ]

top_40_descriptions <- df_top_40$Description
additional_df <- df_remaining[df_remaining$Description %in% top_40_descriptions, ]
positive_df <- rbind(df_top_40, additional_df)
positive_df <- df_top_40

p1 <- ggplot(positive_df, aes(x = Type, y = Description, size = -log(p.adjust), color = NES)) +
  geom_point() +
  scale_size_continuous(range = c(2, 4)) +
  scale_color_continuous() +
  scale_color_continuous("Enrichment Score", limits = c(min(positive_df$NES), max(positive_df$NES)), 
                         type="viridis") +
  theme_bw() +
  labs(x = "Type", y = "Description", size = "-log10(P-value)", color = "Enrichment Score")+
  theme(axis.text.y = element_text(color = "black", face = "bold.italic", size = 8), 
        axis.title = element_text(size = 0))

p1


### Negative bubble plot 
negative_df <- merged_df %>%
  filter(NES < 0) %>%
  arrange(NES)

df_top_40 <- negative_df[1:40, ]
df_remaining <- negative_df[41:nrow(positive_df), ]

top_40_descriptions <- df_top_40$Description
additional_df <- df_remaining[df_remaining$Description %in% top_40_descriptions, ]
negative_df <- rbind(df_top_40, additional_df)

p2 <- ggplot(negative_df, aes(x = Type, y = Description, size = -log(p.adjust), color = NES)) +
  geom_point() +
  scale_size_continuous(range = c(2, 4)) +
  scale_color_continuous("Enrichment Score", limits = c(min(negative_df$NES), max(negative_df$NES)), 
                         type = "viridis") +
  theme_bw() +
  labs(x = "Type", y = "Description", size = "-log10(P-value)", color = "Enrichment Score")+
  theme(axis.text.y = element_text(color = "black", face = "bold.italic", size = 8), 
        axis.title = element_text(size = 0))

p2

ggsave("plots_combined_LC1.pdf", arrangeGrob(p1, p2), device = "pdf", height = 7, width = 6)




