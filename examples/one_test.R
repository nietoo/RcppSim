library(pacman)
p_load(tidyverse, spatstat, future, promises)
library(MathBioSim)

result_dir = './test/'
dir.create(result_dir, showWarnings = FALSE)
dir.create(paste0(result_dir,"pop/"), showWarnings = FALSE)
dir.create(paste0(result_dir,"pcfs/"), showWarnings = FALSE)

plan(sequential)

n_samples = 2000
initial_population = 5000
time_limit = 6000
cell_count_x = 100

params <-
  data.frame(sm=0.96,
             sw=0.09,
             b=0.01,d=0.01,dd=3/(2*log(2)),
             samples=n_samples,
             start_pop=initial_population,
             seed=1234)%>%
  mutate(area=pmax(sm, sw) * 1000)

start_time = Sys.time()
params$death_kernel_r = 10 * params$sw
params$death_kernel_nodes = 1001
x_grid_death = seq(0,params$death_kernel_r,
                   length.out = params$death_kernel_nodes)


params$birth_kernel_r = 10 * params$sm
params$birth_kernel_nodes = 1001
x_grid_birth = seq(0,params$birth_kernel_r,
                   length.out = params$birth_kernel_nodes)

params$area_length_x = params$area
params$init_density = params$start_pop / params$area_length_x

sim_params <-
  list("area_length_x" = params$area_length_x,
       "cell_count_x" = cell_count_x,
       
       "b" = params$b,
       "d" = params$d,
       "dd" = params$dd,
       
       "seed" = params$seed,
       "init_density" = params$init_density,
       
       "death_kernel_r" = params$death_kernel_r,
       "death_kernel_y" = dnorm(x_grid_death, sd = params$sw),
       
       "birth_kernel_r" = params$birth_kernel_r,
       "birth_kernel_y" = dnorm(x_grid_birth, sd = params$sm),
       
       "spline_precision" = 1e-9
  )

pcf_grid = seq(0,max(c(params$sw,params$sm))*10,length.out = 1001)

sim<-new(poisson_1d,sim_params)

pop<-numeric(n_samples)
time<-numeric(n_samples)
pcf_estimate<-list()
calculated_limit = n_samples

for(j in 1:n_samples){
  sim$run_events(sim$total_population)
  pop[j]=sim$total_population
  time[j]=sim$time
  points<-unique.ppp(ppp(sim$get_all_coordinates(),
                         rep(0,length(sim$get_all_coordinates())),
                         c(0,sim$area_length_x),
                         c(-sim$area_length_x/2,sim$area_length_x/2)))
  
  K_estimate<-Kest(points,r=pcf_grid,correction="Ripley")
  
  pcf_estimate[[j]]=data.frame(Kest=K_estimate$iso/2,x=pcf_grid)%>%
    mutate(pfc=(Kest-lag(Kest))/(pcf_grid-lag(pcf_grid))/sim$area_length_x)%>%
    pull(pfc)
  
  if (Sys.time()-start_time>time_limit){
    calculated_limit = j
    break
  }
}

pcf_est_av<-numeric(length(pcf_grid))
for(j in 1:length(pcf_estimate[[1]])){
  jrow=numeric(n_samples)
  for (k in 1:calculated_limit){
    jrow[k]=pcf_estimate[[k]][j]
  }
  pcf_est_av[j]=mean(jrow)
}

pcfs<-data.frame(id=1,r=pcf_grid,y=pcf_est_av)
pops<-data.frame(id=1,time=time,pop=pop)

ggplot(pops %>% select(2:3), aes(x=time, y=pop)) + geom_point(size=0.1)
ggsave("pop.pdf", path = paste0(result_dir, "pop/"))
ggplot(pcfs %>% select(2:3), aes(x=r, y=y)) + geom_point(size=0.1)
ggsave("pcf.pdf", path = paste0(result_dir, "pcfs/"))

write_csv(pops,paste0(result_dir,"pop/pop.csv"))
write_csv(pcfs,paste0(result_dir,"pcfs/pcf.csv"))
write_csv(params,paste0(result_dir,"params.csv"))
