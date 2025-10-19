#Hypotheses presentation
# Load required packages
library(tidyverse)
library(gridExtra)
library(ggExtra)
library(Rmisc)     #data manipulation

###Functional trait space variation ####

# Set seed for reproducibility
set.seed(123)

#Species from island and mainland - two traits shown
n_I <- 100
n_M <- 1000

Trait1_I <- rnorm(n = n_I, mean = 0, 1)
Trait2_I <- rnorm(n = n_I, mean = 0, 1)

Trait1_M <- rnorm(n = n_M, mean = 0, 1)
Trait2_M <- rnorm(n = n_M, mean = 0, 1)


data.Trait<- data.frame(Trait.Value = c(Trait1_I, Trait2_I,
                                        Trait1_M, Trait2_M),
                        Trait.Label = c(rep("Trait1",n_I),rep("Trait2",n_I),
                                        rep("Trait1",n_M),rep("Trait2",n_M)),
                        Group = c(rep("Island", (2*n_I)),
                                  rep("Mainland", (2*n_M))))

# Mainland alone - Not used ####

M <- data.Trait %>%
  filter(Group == "Mainland") %>%
  mutate(id = rep(seq(1,n_M,1),2)) %>% 
  pivot_wider(names_from = Trait.Label, values_from = Trait.Value)

Fig_1 <- ggplot(M, aes(x = Trait1, y = Trait2, color = Group))+
  geom_point(alpha = 0.25, size = 1)+
  geom_density2d()+
  theme_bw()+
  scale_color_manual(values = c("gray"))+
  scale_x_continuous(limits = c(-5,9.5))+
  scale_y_continuous(limits = c(-5,11.5))+
  labs(x = "Dimension 1", y = "Dimension 2", subtitle = "\nMainland functional trait space")+
  geom_point(x = mean(M$Trait1), y = mean(M$Trait2), size = 3, color = "black")+
  theme(legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.text = element_blank())

Fig_A <- ggMarginal(Fig_1, type = "density",
                    groupColour = TRUE, 
                    groupFill = TRUE)

#Hypothesis 1 ####
I_T1 <- data.Trait %>%
  filter(Group == "Island") %>%
  mutate(Trait.Value = Trait.Value*2,
         id = seq(1,n_I,1))

H1 <- data.Trait %>%
  filter(Group != "Island") %>%
  mutate(id = rep(seq(1,n_M,1),2)) %>% 
  full_join(I_T1)%>%
  pivot_wider(names_from = Trait.Label, values_from = Trait.Value)

H1_ICen <- H1 %>% 
  filter(Group == "Island") %>% 
  summarise(Trait1_mu = mean(Trait1), 
            Trait2_mu = mean(Trait2))

H1_MCen <- H1 %>% 
  filter(Group == "Mainland") %>% 
  summarise(Trait1_mu = mean(Trait1), 
            Trait2_mu = mean(Trait2))

Fig2 <- ggplot(H1, aes(x = Trait1, y = Trait2, color = Group))+
  geom_point(alpha = 0.25)+
  geom_density2d()+
  scale_size_manual(values = c(2,1))+
  scale_color_manual(values = c("#2980B9", "gray"))+
  scale_x_continuous(limits = c(-5,9.5))+
  scale_y_continuous(limits = c(-5,11.5))+
  geom_point(x = H1_ICen$Trait1_mu, y = H1_ICen$Trait2_mu, size = 5, color = "#20638f")+
  geom_point(x = H1_MCen$Trait1_mu, y = H1_MCen$Trait2_mu, size = 3, color = "black")+
  labs(x = "Dimension 1", 
       y = "Dimension 2", 
       title = "",
       subtitle = expression(~bold("A")~"    H1 - Functional expansion"),
       color = "")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.75, 0.75), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))

Fig_B <- ggMarginal(Fig2, type = "density",
                    groupColour = TRUE, 
                    groupFill = TRUE)


#Hypothesis 2 ####
I_T2 <- data.Trait %>%
  filter(Group == "Island") %>%
  filter(Trait.Label == "Trait1") %>%
  mutate(Trait.Value = Trait.Value+5,
         id = seq(1,n_I,1))

H2 <- data.Trait %>%
  filter(Group != "Island" | Trait.Label != "Trait1") %>%
  mutate(id = c(seq(1,n_I,1),rep(seq(1,n_M,1),2))) %>% 
  full_join(I_T2)%>%
  pivot_wider(names_from = Trait.Label, values_from = Trait.Value)

H2_ICen <- H2 %>% 
  filter(Group == "Island") %>% 
  summarise(Trait1_mu = mean(Trait1), 
            Trait2_mu = mean(Trait2))

H2_MCen <- H2 %>% 
  filter(Group == "Mainland") %>% 
  summarise(Trait1_mu = mean(Trait1), 
            Trait2_mu = mean(Trait2))

Fig3 <- ggplot(H2, aes(x = Trait1, y = Trait2, color = Group))+
  geom_point(alpha = 0.25)+
  geom_density2d()+
  scale_size_manual(values = c(2,1))+
  scale_color_manual(values = c("#2980B9", "gray"))+
  scale_x_continuous(limits = c(-5,9.5))+
  scale_y_continuous(limits = c(-5,11.5))+
  geom_point(x = H2_ICen$Trait1_mu, y = H2_ICen$Trait2_mu, size = 5, color = "#20638f")+
  geom_point(x = H2_MCen$Trait1_mu, y = H2_MCen$Trait2_mu, size = 3, color = "black")+
  labs(x = "Dimension 1", y = "Dimension 2", 
       title = "",
       subtitle = expression(bold("B")~"    H2 - Functional shift"))+
  theme(legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))

Fig_C <- ggMarginal(Fig3, type = "density",
                    groupColour = TRUE, 
                    groupFill = TRUE)

#Hypothesis 3 ####
I_T3 <- data.Trait %>%
  filter(Group == "Island") %>% 
  mutate(Trait.Value = (Trait.Value*2)+5,
         id = rep(seq(1,n_I,1),2))

H3 <- data.Trait %>%
  filter(Group != "Island") %>%
  mutate(id = rep(seq(1,n_M,1),2)) %>% 
  full_join(I_T3)%>%
  pivot_wider(names_from = Trait.Label, values_from = Trait.Value)

H3_ICen <- H3 %>% 
  filter(Group == "Island") %>% 
  summarise(Trait1_mu = mean(Trait1), 
            Trait2_mu = mean(Trait2))

H3_MCen <- H3 %>% 
  filter(Group == "Mainland") %>% 
  summarise(Trait1_mu = mean(Trait1), 
            Trait2_mu = mean(Trait2))

Fig4 <- ggplot(H3, aes(x = Trait1, y = Trait2, color = Group))+
  geom_point(alpha = 0.25)+
  geom_density2d()+
  scale_size_manual(values = c(2,1))+
  scale_color_manual(values = c("#2980B9", "gray"))+
  scale_x_continuous(limits = c(-5,9.5))+
  scale_y_continuous(limits = c(-5,11.5))+
  geom_point(x = H3_ICen$Trait1_mu, y = H3_ICen$Trait2_mu, size = 5, color = "#20638f")+
  geom_point(x = H3_MCen$Trait1_mu, y = H3_MCen$Trait2_mu, size = 3, color = "black")+
  labs(x = "Dimension 1", y = "Dimension 2", 
       title = "",
       subtitle = expression(bold("C")~"    H3 - Functional shift and expansion"))+
  theme(legend.position = "none", 
        axis.ticks = element_blank(), 
        axis.text = element_blank(), 
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(fill=NA,colour="grey50"))

Fig_D <- ggMarginal(Fig4, type = "density",
                    groupColour = TRUE, 
                    groupFill = TRUE)

#Combine figure ####

Fig_conceptual <- grid.arrange(Fig_B, Fig_C, Fig_D, ncol = 3)
ggsave("Figure 1 - Conceptual.jpg", plot = Fig_conceptual,
       dpi = 450, width = 300, height = 100, units = "mm")
