#14 threads
library(ggplot2)
serial_iter_thousandth <- c(50,50,50,50,50,50,50,50)
parallel_iter_thousandth <- c(7,7,7,7,7,7,7,7,7)
serial_time_thousandth <- c(0.01,0.02,0.01,0.01,0.01,0.02,0.02,0.01,0.01,0.01,0.02,0.02,0.02,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02)
parallel_time_thousandth <- c(0.02,0.03,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
serial_iter_ten_th <- c(54,54)
parallel_iter_ten_th <- c(14,14)
serial_time_ten_th <- c(0.01,0.02,0.02,0.01,0.01,0.02,0.02,0.03,0.01,0.02,0.02,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.01,0.02,0.02)
parallel_time_ten_th <- c(0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.01,0.01,0.01,0.01,0.02,0.01,0.01)
s_i_hund <- c(121)
p_i_hund <- c(47)
s_t_hund <- c(0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.04,0.03,0.03,0.04,0.04,0.04,0.04,0.03,0.04,0.04,0.04,0.03,0.03,0.04,0.04,0.04,0.03,0.03,0.03,0.03,0.03,0.04)
p_t_hund <- c(0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.02,0.02,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.02,0.02,0.01,0.01,0.02,0.01)
s_i_mill <- c(437)
p_i_mill <- c(150,153,163,156,152,156,152,154,154,153,152,151,155,158,154,156,152,154,153,155,158,153,155,155,157,156,154,154,152,155,153,153)
t_4 <-(10000*29)/31
t_6 <-(10000*25)/31
t_8 <-(10000*28)/31
t_10<-10000
t_12<-10000
t_14<-10000
t_16<-(4*20000+27*10000)/31
t_18<-(13*20000+17*10000)/31
i_4 <-c(259,258,272,262,258,265,263,266,259,260,261,269,264,260,260,265,265,270,260,262,267,260,267,268,257,268,258,261,262,263,256,261,265)
i_6 <-c(223,225,220,225,225,221,222,228,225,223,224,225,223,224,226,221,222,220,221,225,219,224,223,224,222,227,219,222,219,223,225,224,223)
i_8 <-c(192,193,195,192,192,193,191,193,195,192,191,195,191,193,189,196,193,192,194,191,191,193,194,194,194,193,197,193,192,195,194,194,191)
i_10<-c(175,177,178,179,178,178,179,175,178,177,179,179,178,179,176,177,181,179,180,178,182,177,178,180,183,176,181,181,176,179,179,179,180)
i_12<-c(169,168,170,166,169,171,169,170,170,174,169,173,170,165,170,169,167,171,170,173,167,171,170,165,174,171,171,169,167,170,167,172,168)
i_16<-c(148,150,145,146,148,150,149,149,154,148,155,148,151,152,155,154,147,146,151,148,151,148,146,151,155,156,144,147,148,148,149,152,150)
i_18<-c(142,138,143,148,146,138,145,140,144,147,149,148,143,143,145,144,143,145,145,144,148,144,142,143,146,143,141,140,138,134,134,139,137)
i_20<-c(145,135,139,138,137,145,137,141,141,134,122,138,136,146,140,136,130,141,134,144,135,140,141,132,145,138,149,138,142,142,136,136,135)
t_20<-c(30000,470000,100000,500000,480000,500000,110000,490000,70000,16730000,210000,240000,350000,520000,280000,20000,500000,130000,16540000,250000,350000,350000,60000,490000,220000,260000,20000,350000,20000,210000,190000,590000,240000)
s_t_mill <- c(0.12,0.13,0.12,0.13,0.13,0.13,0.13,0.14,0.12,0.13,0.13,0.14,0.12,0.12,0.13,0.13,0.13,0.12,0.12,0.12,0.13,0.12,0.13,0.13,0.13,0.14,0.14,0.14,0.13,0.13,0.12)
p_t_mill <- c(0.01,0.02,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.02,0.01,0.01,0.01,0.01,0.01,0.01)
#times <- data.frame(ser_t_thou=mean(serial_time_thousandth),par_t_thou=mean(parallel_time_thousandth),ser_t_ten = mean(serial_time_ten_th),par_t_ten=mean(parallel_time_ten_th),ser_t_hund=mean(s_t_hund),par_t_hund=mean(p_t_hund),ser_t_mill=mean(s_t_mill),par_t_mil=mean(p_t_mill), row.names = NULL, check.rows = FALSE, check.names = TRUE, stringsAsFactors = default.stringsAsFactors())
#titles <-c("S 0.001","P 0.001","S 0.0001","P 0.0001","S 0.00001","P 0.00001","S 0.000001","P 0.000001")
titles <-c("S 10^-3","P 10^-3","S 10^-4","P 10^-4","S 10^-5","P 10^-5","S 10^-6","P 10^-6")
i_vals <-c(50,7,54,14,121,47,437,154)
times <-data.frame( names = titles,
    values = c(mean(serial_time_thousandth),mean(parallel_time_thousandth),mean(serial_time_ten_th),mean(parallel_time_ten_th),mean(s_t_hund),mean(p_t_hund),mean(s_t_mill),mean(p_t_mill))
    )
iters <-data.frame(names_i=titles,values_i=i_vals)
titles_i <- c("4","6","8","10","12","14","16","18","20")
mean_t_c <- c(mean(t_4),mean(t_6),mean(t_8),mean(t_10),mean(t_12),mean(t_14),mean(t_16),mean(t_18),mean(t_20))
mean_i_c <-c(mean(i_4),mean(i_6),mean(i_8),mean(i_10),mean(i_12),mean(p_i_mill),mean(i_16),mean(i_18),mean(i_20))
mean_i_c_t<-c(262.7576, 223.0909, 192.9697 ,178.5152, 169.5455, 154.3125, 149.6667 ,142.6970, 138.4242)
mean_t_c_t<-c( 9354.8 ,   8064.5  ,  9032.3  ,10000  , 10000 ,  10000 ,11290  , 13871,1268788)
core_i <-data.frame(names=titles_i,values=mean_i_c_t)
core_t <-data.frame(names=titles_i,values=mean_t_c_t)
#ggplot(core_t,aes(x=factor(names,level=names),y=values,fill=as.factor(names)))+
#    geom_bar(stat = "identity",show.legend = FALSE)+
#    scale_fill_manual(values = c(rgb(1,0.5,0.9,1),rgb(1,0.5,0.9,1), rgb(1,0.5,0.9,1),rgb(1,0.5,0.9,1), rgb(1,0.5,0.9,1), rgb(1,0,0,1),rgb(1,0.5,0.9,1),rgb(1,0.5,0.9,1),rgb(1,0.5,0.9,1))) +
#    geom_text(aes(x=names, y=values, label=values))+
#    labs(title="Parallel Stochastic Gradient Descent Core Counts", subtitle="Cycles Plot", y="Average Cycles", x="Increasing Core Count on 0.000001", caption="PDC")
#mean_t_c
#ggplot(core_i,aes(x=factor(names,level=names),y=values,fill=as.factor(names)))+
#    geom_bar(stat = "identity",show.legend = FALSE)+
#    scale_fill_manual(values = c(rgb(1,0.5,0.4,1),rgb(1,0.5,0.5,1), rgb(1,0.5,0.6,1),rgb(1,0.5,0.7,1), rgb(1,0.5,0.8,1), rgb(1,0.5,0.9,1),rgb(1,0.5,0.1,1),rgb(1,0.5,0.2,1),rgb(1,0.5,0.3,1))) +
#    geom_text(aes(x=names, y=values, label=values))+
#    labs(title="Parallel Stochastic Gradient Descent Core Counts", subtitle="Iteration Plot", y="Average Iterations", x="Increasing Core Count on 0.000001", caption="PDC")
ggplot(iters,aes(x=factor(names_i,level=names_i),y=values_i,fill=as.factor(names_i)))+
    geom_bar(stat = "identity",show.legend = FALSE)+
    scale_fill_manual(values = c(rgb(1,0.5,1,1),rgb(1,0.5,1,1), rgb(1,0.5,1,1),rgb(1,0.5,1,1), rgb(1,0.5,0.2,1), rgb(1,0.5,0.2,1),rgb(1,0.5,0.2,1),rgb(1,0.5,0.2,1))) +
    geom_text(aes(x=names_i, y=values_i, label=values_i))+
    labs(title="Serial vs Parallel Stochastic Gradient Descent", subtitle="Iteration Plot", y="Average Iterations", x="Serial/Parallel with decreasing Learning Step", caption="PDC")

#trunc <-c(.01354,.01161,.01516,.01161,.03677,.01290,.12838,.01129)
#ggplot(times,aes(x=factor(names,level=names), y=values,fill=as.factor(names)))+geom_bar(stat = "identity",show.legend = FALSE)+
#    scale_fill_manual(values = c(rgb(1,0.5,1,1),rgb(1,0.5,1,1), rgb(1,0.5,1,1),rgb(1,0.5,1,1), rgb(1,0.5,0.2,1), rgb(1,0.5,0.2,1),rgb(1,0.5,0.2,1),rgb(1,0.5,0.2,1))) +
#    geom_text(aes(x=names, y=values, label=trunc))+
#    labs(title="Serial vs Parallel Stochastic Gradient Descent", subtitle="Runtime Plot", y="Average Time in seconds", x="Serial/Parallel with decreasing Learning Step", caption="PDC")

mean(serial_time_thousandth)
mean(parallel_time_thousandth)
mean(serial_time_ten_th)
mean(parallel_time_ten_th)
mean(s_t_hund)
mean(p_t_hund)
mean(p_i_mill)
mean(s_t_mill)
mean(p_t_mill)


