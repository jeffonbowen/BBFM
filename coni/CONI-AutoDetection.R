#############################
## ARU Automated Detection ##
#############################

## Calculate precision, recall and F-score for varying levels of score to simulate score threshold.

library(tidyverse)

# Import benchmark validation results
bench.dat <- read_csv("bench-validation.csv")
bench.dat <- bench.dat %>% rename(species = "TOP1MATCH*")

#bench.dat <- dplyr::filter(bench.dat, species == "CONI")  # make sure only CONI records shown
#bench.dat <- dplyr::filter(bench.dat, CHANNEL == 1)


table(bench.dat$species)      # quick check of the data
table(bench.dat$validation)

#####
## Distribution of TP and FP
#####

rm(list = ls())   # freshstart
library(tidyverse)

# Import benchmark validation results
bench.dat <- read_csv("bench-validation.csv")
bench.dat <- bench.dat %>% rename(species = "TOP1MATCH*")
# Add standardized score
bench.dat <- dplyr::filter(bench.dat, species == "CONI")  # make sure only CONI records shown

bench.dat$score = 1-((bench.dat$TOP1DIST-min(bench.dat$TOP1DIST))/
                       (max(bench.dat$TOP1DIST)-min(bench.dat$TOP1DIST)))
summary(bench.dat$score)


bench.dat <- dplyr::filter(bench.dat, CHANNEL == 0)

# Freq histogram of true positive and false positive CONI hits.
g.dist <- ggplot(data=bench.dat, aes(x=TOP1DIST, colour=validation))+
  geom_histogram(fill="white", bins=50, position="dodge")+
  theme_bw()+
  xlab(label = "Distance to Cluster Centre")+
  ylab(label = "CONI Recognnizer Hits")+
  theme(legend.position = c(0.45, 0.7))+
  scale_color_manual(name="Validation", labels=c("False Positive", "True Positive"), values=c("blue", "red"))+
  theme(axis.title=element_text(size=9))
g.dist

ggsave("CONI dist distribution.png", g.dist, width = 5, height = 3.5)


#####
## Simple Calc With All Distances
#####

rm(list = ls())   # freshstart
library(tidyverse)

# Import benchmark validation results
bench.dat <- read_csv("bench-validation.csv")
bench.dat <- bench.dat %>% rename(species = "TOP1MATCH*")
# Add standardized score
bench.dat$score = 1-((bench.dat$TOP1DIST-min(bench.dat$TOP1DIST))/
                       (max(bench.dat$TOP1DIST)-min(bench.dat$TOP1DIST)))
summary(bench.dat$score)

#bench.dat <- dplyr::filter(bench.dat, species == "CONI")  # make sure only CONI records shown

bench.dat <- dplyr::filter(bench.dat, CHANNEL == 1)

# Make a blank table to append results to.
metrics <- tibble(thr.dist = double(),
                  tp = integer(), 
                  fp = integer(), 
                  fn = integer(),
                  precision = double(), 
                  recall = double(),
                  fscore = double())

beta <- 1
i <- 1
tp <- 0
fp <- 0
fn <- 0

  valid.cnt <- table(bench.dat$validation)
  tp <- valid.cnt[4]
  fp <- valid.cnt[2]
  fn <- valid.cnt[1]
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  fscore <- ((beta^2+1) * precision * recall) / (beta^2*precision+recall)
  thr.dist <- i
  metrics.temp <- tibble(thr.dist, tp, fp, fn, precision, recall, fscore)
  metrics <- bind_rows(metrics,metrics.temp)



###########
##Using Distance
###########

rm(list = ls())   # freshstart
library(tidyverse)

# Import benchmark validation results
bench.dat <- read_csv("bench-validation.csv")
bench.dat <- bench.dat %>% rename(species = "TOP1MATCH*")
bench.dat$validation <- factor(bench.dat$validation)
bench.dat <- dplyr::filter(bench.dat, CHANNEL == 0)
bench.dat <- dplyr::filter(bench.dat, species == "CONI")

summary(bench.dat$validation)

summary(bench.dat$TOP1DIST)

bench.dat.thresh <- bench.dat %>% filter(TOP1DIST <= 0.3)
valid.cnt <- table(bench.dat.thresh$validation)
table(bench.dat.thresh$validation)
valid.cnt[names(valid.cnt) == "tp"]
valid.cnt[names(valid.cnt) == "fp"]

table(bench.dat$validation)

# Add standardized score
bench.dat$invdist = 1/bench.dat$TOP1DIST
bench.dat$score = (bench.dat$invdist-min(bench.dat$invdist))/(max(bench.dat$invdist)-min(bench.dat$invdist))

# Make a blank table to append results to.
metrics <- tibble(thr.dist = double(),
                  tp = integer(), 
                  fp = integer(), 
                  fn = integer(),
                  precision = double(), 
                  HL = integer(),
                  recall = double(),
                  fscore = double())

# Freq histogram of true positive and false positive CONI hits.
g.dist <- ggplot(data=bench.dat, aes(x=TOP1DIST, colour=validation))+
  geom_histogram(fill="white", bins=50, position="dodge")+
  theme_bw()+
  xlab(label = "Distance")+
  ylab(label = "CONI Recognnizer Hits")+
  theme(legend.position = c(0.2, 0.8))+
  scale_color_manual(name="Validation", labels=c("False Positive", "True Positive"), values=c("blue", "red"))
g.dist

HL <- 3814    # change this later
beta <- 1

dist.range <- seq(1.1, 0.2, -0.01)
i <- 0
tp <- 0
fp <- 0
fn <- 978

for (i in dist.range)
{
  bench.dat.thresh <- bench.dat %>% filter(TOP1DIST <= i)
  valid.cnt <- table(bench.dat.thresh$validation)
  tp <- valid.cnt[names(valid.cnt) == "tp"]
  fp <- valid.cnt[names(valid.cnt) == "fp"]
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  fscore <- ((beta^2+1) * precision * recall) / (beta^2*precision+recall)
  thr.dist <- i
  metrics.temp <- tibble(thr.dist, tp, fp, fn, precision, HL, recall, fscore)
  metrics <- bind_rows(metrics,metrics.temp)
}

# metrics <- metrics %>% replace_na(list(tp=0, fp=0, precision=0, recall=0, fscore=0))

g.precision <- ggplot(metrics, aes(x=thr.dist, y=precision))+
  geom_line(size=0.75)+
  xlab(label = "Distance") +
  ylab(label = "Precision")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.precision

g.recall <- ggplot(metrics, aes(x=thr.dist, y=recall))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "Distance") +
  ylab(label = "Recall")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.recall

g.fscore <- ggplot(metrics, aes(x=thr.dist, y=fscore))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "Distance") +
  ylab(label = "F-Score")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.fscore

g.pr <- ggplot(metrics, aes(x=recall, y=precision))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "Recall (relative to Human Listening") +
  ylab(label = "Precision")+
  xlim(0,1)+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.pr

library(gridExtra)
g.arrange.distance <- grid.arrange(g.precision, g.recall, g.fscore, nrow = 1)

ggsave("Metrics by distance.png", g.arrange.distance, width = 7, height = 4)


###########
##Using Score
###########


# Make a blank table to append results to.
metrics <- tibble(thr.dist = double(),
                  tp = integer(), 
                  fp = integer(), 
                  fn = integer(),
                  precision = double(), 
                  HL = integer(),
                  recall = double(),
                  fscore = double())

HL <- 3814    # change this later
beta <- 1

dist.range <- seq(0, 1, 0.01)
i <- 0
tp <- 0
fp <- 0
fn <- 978

for (i in dist.range)
{
  bench.dat.thresh <- bench.dat %>% filter(score >= i)
  valid.cnt <- table(bench.dat.thresh$validation)
  tp <- valid.cnt[names(valid.cnt) == "tp"]
  fp <- valid.cnt[names(valid.cnt) == "fp"]
  precision <- tp/(tp+fp)
  recall <- tp/(tp+fn)
  fscore <- ((beta^2+1) * precision * recall) / (beta^2*precision+recall)
  thr.dist <- i
  metrics.temp <- tibble(thr.dist, tp, fp, fn, precision, HL, recall, fscore)
  metrics <- bind_rows(metrics,metrics.temp)
}

# metrics <- metrics %>% replace_na(list(tp=0, fp=0, precision=0, recall=0, fscore=0))

g.precision2 <- ggplot(metrics, aes(x=thr.dist, y=precision))+
  geom_line(size=0.75)+
  xlab(label = "Score") +
  ylab(label = "Precision")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.precision2

g.recall2 <- ggplot(metrics, aes(x=thr.dist, y=recall))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "Score") +
  ylab(label = "Recall")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.recall2

g.fscore2 <- ggplot(metrics, aes(x=thr.dist, y=fscore))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "Score") +
  ylab(label = "F-Score")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.fscore2

g.pr2 <- ggplot(metrics, aes(x=recall, y=precision))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "Recall (relative to Human Listening") +
  ylab(label = "Precision")+
  xlim(0,1)+
  ylim(0,1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title=element_text(size=9))
g.pr2

library(gridExtra)
g.arrange.score <- grid.arrange(g.precision2, g.recall2, g.fscore2, nrow = 1)

ggsave("Metrics by distance.png", g.arrange.score, width = 7, height = 4)

g.arrange.all <- grid.arrange(g.precision, g.recall, g.fscore, g.precision2, g.recall2, g.fscore2, nrow = 2)

ggsave("Metrics by both.png", g.arrange.all, width = 7, height = 6)


#############
##PRROC
#############

library(PRROC)
library(ROCR)


## TEST
x <- rnorm( 1000 );
y <- rnorm( 1000, -1 );
# compute area under PR curve for the hard-labeled case
pr <- pr.curve( x, y, curve=TRUE );
print( pr );
plot(pr)


## Malke ROC Using PRROC

rm(list = ls())   # freshstart
library(tidyverse)

# Import benchmark validation results
bench.dat <- read_csv("bench-validation.csv")
bench.dat <- bench.dat %>% rename(species = "TOP1MATCH*")
bench.dat <- dplyr::filter(bench.dat, CHANNEL == 0)
# Add standardized score
bench.dat$score = 1-((bench.dat$TOP1DIST-min(bench.dat$TOP1DIST))/
                       (max(bench.dat$TOP1DIST)-min(bench.dat$TOP1DIST)))
summary(bench.dat$score)

fg.1 <- dplyr::filter(bench.dat, species == "CONI")
fg <- fg.1$TOP1DIST

bg.1 <- dplyr::filter(bench.dat, species == "Other")
bg <- bg.1$TOP1DIST

pr <- pr.curve(fg, bg, curve=TRUE )
print(pr)
plot(pr)

roc <- roc.curve(fg, bg, curve=TRUE )
print(roc)
plot(roc)

roc.dat <- as.data.frame(roc$curve)

g.roc2 <- ggplot(roc.dat, aes(x=V1, y=V2))+
  geom_line()+
  geom_line(size=0.75)+
  xlab(label = "False Positive Rate") +
  ylab(label = "True Positive Rate")+
  ylim(0, 1.0)+
  theme_bw()+
  theme(legend.position = "none")+
  annotate(geom="text", x=.75, y=0.30, label="Area Under Curve = 0.6498",
         color="black")
g.roc2

ggsave("ROC.png", g.roc2, width = 5, height = 3.5)


## Using ROCR

rm(list = ls())   # freshstart
library(ROCR)
library(tidyverse)

# Import benchmark validation results
bench.dat <- read_csv("bench-validation.csv")
bench.dat <- bench.dat %>% rename(species = "TOP1MATCH*")
bench.dat <- dplyr::filter(bench.dat, CHANNEL == 1)
bench.dat <- dplyr::filter(bench.dat, species == "CONI")  # make sure only CONI records shown

# Add standardized score
bench.dat$score = 1-((bench.dat$TOP1DIST-min(bench.dat$TOP1DIST))/
                       (max(bench.dat$TOP1DIST)-min(bench.dat$TOP1DIST)))

prediction <- bench.dat$score

labels <- bench.dat$validation



labels <- recode(labels, tp = "CONI", fp = "other")

table(prediction)

table(labels)

p <- ROCR::prediction(prediction, labels)

q <- performance(p, "prec")

plot(q)



p <- prediction(ROCR.simple$predictions, ROCR.simple$labels)
q <- performance(p, "acc")
plot(q)


