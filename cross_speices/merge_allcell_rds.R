rm(list = ls())
setwd('/data/liangxian/project/MHD_hypoMap/down_sample/')

library(Seurat)
library(tidyverse)
library(patchwork)

mice.neu <- readRDS('data/mus_neu_smaple.rds')
monkey.neu <- readRDS('data/monkey_neu_smaple.rds')
mice.nn <- readRDS('data/mus_nonneu_smaple.rds')
monkey.nn <- readRDS('data/monkey_nonneu_smaple.rds')

mice.nn$integrated_level <- 'mice_non_neuron'
monkey.nn$integrated_level <- 'monkey_non_neuron'
mice.neu$integrated_level <- 'mice_neuron'
monkey.neu$integrated_level <- 'monkey_neuron'

merge.neu <- merge(mice.neu, y = c(monkey.neu, monkey.nn, mice.nn))

saveRDS(merge.neu, 'data/allcell_seu.rds')
