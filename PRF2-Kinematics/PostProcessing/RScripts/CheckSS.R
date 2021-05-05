library(tidyverse)

#Data Import -------------
posData1_0 <- read_csv("data/SS_RSL1.0_tol_1.csv")
posData1_1 <- read_csv("data/SS_RSL1.1_tol_1.csv")
posData1_2 <- read_csv("data/SS_RSL1.2_tol_1.csv")
posData1_3 <- read_csv("data/SS_RSL1.3_tol_1.csv")
posData1_4 <- read_csv("data/SS_RSL1.4_tol_1.csv")
posData1_5 <- read_csv("data/SS_RSL1.5_tol_1.csv")
posData1_6 <- read_csv("data/SS_RSL1.6_tol_1.csv")
posData1_8 <- read_csv("data/SS_RSL1.8_tol_1.csv")
posData2_0 <- read_csv("data/SS_RSL2.0_tol_1.csv")

posData <-mutate(posData2_0, visc = 2.0*2.0*pi*10.0*amp*0.15*0.8/Re)
posData <- arrange(posData,desc(amp),desc(Re))
posData <- select(posData, Re, amp, err, t_ss, v_ss)
A0.4 <- filter(posData, amp == 0.12) %>%
  select(Re,err,t_ss,v_ss) %>%
  arrange(Re)
A0.5 <- filter(posData, amp == 0.15) %>%
  select(Re,err,t_ss,v_ss) %>%
  arrange(Re)
A0.6 <- filter(posData, amp == 0.18) %>%
  select(Re,err,t_ss,v_ss) %>%
  arrange(Re)
A0.7 <- filter(posData, amp == 0.21) %>%
  select(Re,err,t_ss,v_ss) %>%
  arrange(Re)
A0.8 <- filter(posData, amp == 0.24) %>%
  select(Re,err,t_ss,v_ss) %>%
  arrange(Re)
A0.8
