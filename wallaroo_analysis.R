library(dplyr)
library(ggplot2)
library(kableExtra)
library(lubridate)
library(Rmisc)
library(bookdown)
library(knitr)
library(tinytex)
library(citr)
library(reshape2)
data = read.csv("wallaroo.csv")
data[,1] = NULL
data[is.na(data)] = 0


data = data %>%
  mutate(minutes = 0,
         date = as.POSIXct(data$date, format =  "%m/%d/%Y"))

data = data %>%
  mutate(time = paste(data$hour, data$minutes, sep=":"))


data$newdate =  as.POSIXct(paste(data$date, data$time),
                           format = "%Y-%m-%d %H:%M",
                           tz = "GMT")

space.data = data %>%
  mutate(day.no = yday(newdate),
         species = as.factor(Species))


####2. SPI ####
# Kangs only

spi.t1 = space.data %>%
  filter(day.no >= 341 & day.no <= 350) %>%
  group_by(day.no) %>%
  dplyr::summarise(Z1 = sum(Z1),
            Z2 = sum(Z2),
            Z3 = sum(Z3),
            Z4 = sum(Z4),
            Z5 = sum(Z5),
            Z6 = sum(Z6),
            total = sum(Z1, Z2, Z3, Z4, Z5, Z6),
            species = first(species),
            day.no = first(day.no),
            e1 = total/6,
            e2 = total/6,
            e3 = total/6,
            e4 = total/6,
            e5 = total/6,
            e6 = total/6,
            f1 = Z1-e1,
            f2 = Z2-e2,
            f3 = Z3-e3,
            f4 = Z4-e4,
            f5 = Z5-e5,
            f6 = Z6-e6,
            sum = sum(abs(f1),abs(f2),abs(f3),abs(f4),abs(f5),abs(f6)),
            base = 2*(total-e2),
            spi = sum / base,
            z1e = Z1/total,
            z2e = Z2/total,
            z3e = Z3/total,
            z4e = Z4/total,
            z5e = Z5/total,
            z6e = Z6/total) %>%
  select(day.no, species, z1e,z2e,z3e,z4e,z5e,z6e, spi) %>%
  mutate(treatment = as.numeric(1))

# Kangs and BW
spi.t2 = space.data %>%
  filter(day.no >= 351 & day.no <= 360) %>%
  group_by(species, day.no) %>%
  dplyr::summarise(Z1 = sum(Z1),
            Z2 = sum(Z2),
            Z3 = sum(Z3),
            Z4 = sum(Z4),
            Z5 = sum(Z5),
            Z6 = sum(Z6),
            total = sum(Z1, Z2, Z3, Z4, Z5, Z6),
            species = first(species),
            day.no = first(day.no),
            e1 = total/6,
            e2 = total/6,
            e3 = total/6,
            e4 = total/6,
            e5 = total/6,
            e6 = total/6,
            f1 = Z1-e1,
            f2 = Z2-e2,
            f3 = Z3-e3,
            f4 = Z4-e4,
            f5 = Z5-e5,
            f6 = Z6-e6,
            sum = sum(abs(f1),abs(f2),abs(f3),abs(f4),abs(f5),abs(f6)),
            base = 2*(total-e2),
            spi = sum / base,
            z1e = Z1/total,
            z2e = Z2/total,
            z3e = Z3/total,
            z4e = Z4/total,
            z5e = Z5/total,
            z6e = Z6/total) %>%
  select(day.no, species, z1e,z2e,z3e,z4e,z5e,z6e, spi) %>%
  mutate(treatment = as.numeric(2))


# Kangs, BW and SW


spi.t3 = space.data %>%
  filter(day.no >= 58 & day.no <= 64) %>%
  group_by(species,day.no) %>%
  dplyr::summarise(Z1 = sum(Z1),
            Z2 = sum(Z2),
            Z3 = sum(Z3),
            Z4 = sum(Z4),
            Z5 = sum(Z5),
            Z6 = sum(Z6),
            total = sum(Z1, Z2, Z3, Z4, Z5, Z6),
            species = first(species),
            day.no = first(day.no),
            e1 = total/6,
            e2 = total/6,
            e3 = total/6,
            e4 = total/6,
            e5 = total/6,
            e6 = total/6,
            f1 = Z1-e1,
            f2 = Z2-e2,
            f3 = Z3-e3,
            f4 = Z4-e4,
            f5 = Z5-e5,
            f6 = Z6-e6,
            sum = sum(abs(f1),abs(f2),abs(f3),abs(f4),abs(f5),abs(f6)),
            base = 2*(total-e2),
            spi = sum / base,
            z1e = Z1/total,
            z2e = Z2/total,
            z3e = Z3/total,
            z4e = Z4/total,
            z5e = Z5/total,
            z6e = Z6/total) %>%
  select(day.no, species, z1e,z2e,z3e,z4e,z5e,z6e, spi) %>%
  mutate(treatment = as.numeric(3))

spi.df = rbind(spi.t1, spi.t2, spi.t3)

spi.df %>%
  group_by(species, treatment) %>%
  dplyr::summarise(CI = CI(spi)[2]-CI(spi)[3],
                   mean = mean(spi),
                   upper = mean + CI,
                   lower = mean - CI)

#### 3. Electivity index ####

## Kangs only
ei.t1 = space.data %>%
  filter(day.no >= 341 & day.no <= 350)%>%
  group_by(day.no) %>%
  dplyr::summarise(z1 = sum(Z1),
            z2 = sum(Z2),
            z3 = sum(Z3),
            z4 = sum(Z4),
            z5 = sum(Z5),
            z6 = sum(Z6),
            total = sum(Z1, Z2, Z3, Z4, Z5, Z6),
            species = first(species),
            day.no = first(day.no),
            e1 = total/6,
            e2 = total/6,
            e3 = total/6,
            e4 = total/6,
            e5 = total/6,
            e6 = total/6,
            w1 = z1 / e1,
            w2 = z2 / e2,
            w3 = z3 / e3,
            w4 = z4 / e4,
            w5 = z5 / e5,
            w6 = z6/e6,
            w1a = w1/sum(w1,w2,w3,w4,w5,w6),
            w2a = w2/sum(w1,w2,w3,w4,w5,w6),
            w3a = w3/sum(w1,w2,w3,w4,w5,w6),
            w4a = w4/sum(w1,w2,w3,w4,w5,w6),
            w5a = w5/sum(w1,w2,w3,w4,w5,w6),
            w6a = w6/sum(w1,w2,w3,w4,w5,w6),
            ew1 = (w1a-(1/6))/(w1a+(1/6)),
            ew2 = (w2a-(1/6))/(w2a+(1/6)),
            ew3 = (w3a-(1/6))/(w3a+(1/6)),
            ew4 = (w4a-(1/6))/(w4a+(1/6)),
            ew5 = (w5a-(1/6))/(w5a+(1/6)),
            ew6 = (w6a-(1/6))/(w6a+(1/6))) %>%
  select(species,day.no, w1a,w2a,w3a,w4a,w5a,w6a,ew1,ew2,ew3,ew4,ew5,ew6) %>%
  mutate(treatment = as.numeric(1))

## Kangs + BW
ei.t2 = space.data %>%
  filter(day.no >= 351 & day.no <= 360)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(z1 = sum(Z1),
                   z2 = sum(Z2),
                   z3 = sum(Z3),
                   z4 = sum(Z4),
                   z5 = sum(Z5),
                   z6 = sum(Z6),
                   total = sum(Z1, Z2, Z3, Z4, Z5, Z6),
                   species = first(species),
                   day.no = first(day.no),
                   e1 = total/6,
                   e2 = total/6,
                   e3 = total/6,
                   e4 = total/6,
                   e5 = total/6,
                   e6 = total/6,
                   w1 = z1 / e1,
                   w2 = z2 / e2,
                   w3 = z3 / e3,
                   w4 = z4 / e4,
                   w5 = z5 / e5,
                   w6 = z6/e6,
                   w1a = w1/sum(w1,w2,w3,w4,w5,w6),
                   w2a = w2/sum(w1,w2,w3,w4,w5,w6),
                   w3a = w3/sum(w1,w2,w3,w4,w5,w6),
                   w4a = w4/sum(w1,w2,w3,w4,w5,w6),
                   w5a = w5/sum(w1,w2,w3,w4,w5,w6),
                   w6a = w6/sum(w1,w2,w3,w4,w5,w6),
                   ew1 = (w1a-(1/6))/(w1a+(1/6)),
                   ew2 = (w2a-(1/6))/(w2a+(1/6)),
                   ew3 = (w3a-(1/6))/(w3a+(1/6)),
                   ew4 = (w4a-(1/6))/(w4a+(1/6)),
                   ew5 = (w5a-(1/6))/(w5a+(1/6)),
                   ew6 = (w6a-(1/6))/(w6a+(1/6))) %>%
  select(species,day.no, w1a,w2a,w3a,w4a,w5a,w6a,ew1,ew2,ew3,ew4,ew5,ew6) %>%
  mutate(treatment = as.numeric(2))

## Kangs + BW + SW
ei.t3 = space.data %>%
  filter(day.no >= 58 & day.no <= 64)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(z1 = sum(Z1),
                   z2 = sum(Z2),
                   z3 = sum(Z3),
                   z4 = sum(Z4),
                   z5 = sum(Z5),
                   z6 = sum(Z6),
                   total = sum(Z1, Z2, Z3, Z4, Z5, Z6),
                   species = first(species),
                   day.no = first(day.no),
                   e1 = total/6,
                   e2 = total/6,
                   e3 = total/6,
                   e4 = total/6,
                   e5 = total/6,
                   e6 = total/6,
                   w1 = z1 / e1,
                   w2 = z2 / e2,
                   w3 = z3 / e3,
                   w4 = z4 / e4,
                   w5 = z5 / e5,
                   w6 = z6/e6,
                   w1a = w1/sum(w1,w2,w3,w4,w5,w6),
                   w2a = w2/sum(w1,w2,w3,w4,w5,w6),
                   w3a = w3/sum(w1,w2,w3,w4,w5,w6),
                   w4a = w4/sum(w1,w2,w3,w4,w5,w6),
                   w5a = w5/sum(w1,w2,w3,w4,w5,w6),
                   w6a = w6/sum(w1,w2,w3,w4,w5,w6),
                   ew1 = (w1a-(1/6))/(w1a+(1/6)),
                   ew2 = (w2a-(1/6))/(w2a+(1/6)),
                   ew3 = (w3a-(1/6))/(w3a+(1/6)),
                   ew4 = (w4a-(1/6))/(w4a+(1/6)),
                   ew5 = (w5a-(1/6))/(w5a+(1/6)),
                   ew6 = (w6a-(1/6))/(w6a+(1/6))) %>%
  select(species,day.no, w1a,w2a,w3a,w4a,w5a,w6a,ew1,ew2,ew3,ew4,ew5,ew6) %>%
  mutate(treatment = as.numeric(3))


ei.df = rbind(ei.t1, ei.t2, ei.t3)

#### 4. Activity budget ####

behave = data %>%
  mutate(total = Total.1,
         species = as.factor(Species),
         day.no = yday(newdate))

behave = behave[,c(5:16,26:29)]


# K only
ab.t1 = behave %>%
  filter(day.no >= 341 & day.no <= 350)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(Re = sum(Re),
            Vi = sum(Vi),
            Lo = sum(Lo),
            Dr = sum(Dr),
            Fe = sum(Fe),
            Gr = sum(Gr),
            N = sum(N),
            P = sum(P),
            Nso = sum(Nso),
            Pso = sum(PSo),
            O = sum(O),
            total = sum(Re, Vi, Lo, Dr, Fe, Gr, N, P, Nso, Pso, O),
            Re.prop = sum(Re) / total,
            Vi.prop = sum(Vi)/total,
            Lo.prop = sum(Lo)/total,
            Dr.prop = sum(Dr)/total,
            Fe.prop = sum(Fe)/total,
            Gr.prop = sum(Gr)/total,
            N.prop = sum(N)/total,
            P.prop = sum(P)/total,
            Nso.prop = sum(Nso)/total,
            Pso.prop = sum(PSo)/total,
            O.prop = sum(O)/total,
            species = first(species)) %>%
  mutate_if(is.numeric, ~ . * 100) %>%
  select(Re.prop, Vi.prop, Lo.prop, Dr.prop,
         Fe.prop,Gr.prop, N.prop, P.prop, Nso.prop,
         Pso.prop, O.prop, species) %>%
  mutate(treatment = as.numeric(1))


# K + BW
ab.t2 = behave %>%
  filter(day.no >= 351 & day.no <= 360)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(Re = sum(Re),
            Vi = sum(Vi),
            Lo = sum(Lo),
            Dr = sum(Dr),
            Fe = sum(Fe),
            Gr = sum(Gr),
            N = sum(N),
            P = sum(P),
            Nso = sum(Nso),
            Pso = sum(PSo),
            O = sum(O),
            total = sum(Re, Vi, Lo, Dr, Fe, Gr, N, P, Nso, Pso, O),
            Re.prop = sum(Re) / total,
            Vi.prop = sum(Vi)/total,
            Lo.prop = sum(Lo)/total,
            Dr.prop = sum(Dr)/total,
            Fe.prop = sum(Fe)/total,
            Gr.prop = sum(Gr)/total,
            N.prop = sum(N)/total,
            P.prop = sum(P)/total,
            Nso.prop = sum(Nso)/total,
            Pso.prop = sum(PSo)/total,
            O.prop = sum(O)/total,
            species = first(species)) %>%
  mutate_if(is.numeric, ~ . * 100) %>%
  select(Re.prop, Vi.prop, Lo.prop, Dr.prop,
         Fe.prop,Gr.prop, N.prop, P.prop, Nso.prop,
         Pso.prop, O.prop, species) %>%
  mutate(treatment = as.numeric(2))
ab.t2[is.na(ab.t2)] = 0


# K + BW + SW
ab.t3 = behave %>%
  filter(day.no >= 58 & day.no <= 64)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(Re = sum(Re),
            Vi = sum(Vi),
            Lo = sum(Lo),
            Dr = sum(Dr),
            Fe = sum(Fe),
            Gr = sum(Gr),
            N = sum(N),
            P = sum(P),
            Nso = sum(Nso),
            Pso = sum(PSo),
            O = sum(O),
            total = sum(Re, Vi, Lo, Dr, Fe, Gr, N, P, Nso, Pso, O),
            Re.prop = sum(Re) / total,
            Vi.prop = sum(Vi)/total,
            Lo.prop = sum(Lo)/total,
            Dr.prop = sum(Dr)/total,
            Fe.prop = sum(Fe)/total,
            Gr.prop = sum(Gr)/total,
            N.prop = sum(N)/total,
            P.prop = sum(P)/total,
            Nso.prop = sum(Nso)/total,
            Pso.prop = sum(PSo)/total,
            O.prop = sum(O)/total,
            species = first(species)) %>%
  mutate_if(is.numeric, ~ . * 100) %>%
  select(Re.prop, Vi.prop, Lo.prop, Dr.prop,
         Fe.prop,Gr.prop, N.prop, P.prop, Nso.prop,
         Pso.prop, O.prop, species) %>%
  mutate(treatment = as.numeric(3))
ab.t3[is.na(ab.t3)] = 0

ab.df = rbind(ab.t1, ab.t2, ab.t3)

#### 5. Shannon's Behavioural index ####

# Rest is not included a recommended by Miller et al. 2020
# OOS is removed; removal of these two results in a smaller "total"
# Because NSO PSO and O are almost non-existent, they are also removed.
# Only retain Vi, Lo, Fe, Gr

# Kangs only
shannon.t1 = behave %>%
  filter(day.no >= 341 & day.no <= 350)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(Vi = sum(Vi),
            Lo = sum(Lo),
            Fe = sum(Fe),
            Gr = sum(Gr),
            total = sum(Vi, Lo, Fe, Gr),
            Vi.prop = sum(Vi)/total,
            Lo.prop = sum(Lo)/total,
            Fe.prop = sum(Fe)/total,
            Gr.prop = sum(Gr)/total,
            Vi.H = Vi.prop*log(Vi.prop),
            Lo.H = Lo.prop*log(Lo.prop),
            Fe.H = Fe.prop*log(Fe.prop),
            Gr.H = Gr.prop*log(Gr.prop))%>%
  mutate(treatment = as.numeric(1))

shannon.t1[is.na(shannon.t1)] = 0


# Kangs + BW only
shannon.t2 = behave %>%
  filter(day.no >= 351 & day.no <= 360)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(Vi = sum(Vi),
                   Lo = sum(Lo),
                   Fe = sum(Fe),
                   Gr = sum(Gr),
                   total = sum(Vi, Lo, Fe, Gr),
                   Vi.prop = sum(Vi)/total,
                   Lo.prop = sum(Lo)/total,
                   Fe.prop = sum(Fe)/total,
                   Gr.prop = sum(Gr)/total,
                   Vi.H = Vi.prop*log(Vi.prop),
                   Lo.H = Lo.prop*log(Lo.prop),
                   Fe.H = Fe.prop*log(Fe.prop),
                   Gr.H = Gr.prop*log(Gr.prop))%>%
  mutate(treatment = as.numeric(2))

shannon.t2[is.na(shannon.t2)] = 0

# Kangs + BW + SW
shannon.t3 = behave %>%
  filter(day.no >= 58 & day.no <= 64)%>%
  group_by(species, day.no) %>%
  dplyr::summarise(Vi = sum(Vi),
                   Lo = sum(Lo),
                   Fe = sum(Fe),
                   Gr = sum(Gr),
                   total = sum(Vi, Lo, Fe, Gr),
                   Vi.prop = sum(Vi)/total,
                   Lo.prop = sum(Lo)/total,
                   Fe.prop = sum(Fe)/total,
                   Gr.prop = sum(Gr)/total,
                   Vi.H = Vi.prop*log(Vi.prop),
                   Lo.H = Lo.prop*log(Lo.prop),
                   Fe.H = Fe.prop*log(Fe.prop),
                   Gr.H = Gr.prop*log(Gr.prop))%>%
  mutate(treatment = as.numeric(3))

shannon.t3[is.na(shannon.t3)] = 0

si.df = rbind(shannon.t1, shannon.t2, shannon.t3)


