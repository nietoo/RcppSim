library(pacman)                                    
p_load(tidyverse,spatstat,future,promises,listenv)  
library(MathBioSim)                                 

result_dir = './16.05.2020_grid/'
dir.create(result_dir, showWarnings = FALSE)
dir.create(paste0(result_dir,"pop/"), showWarnings = FALSE)
dir.create(paste0(result_dir,"pop_sma/"), showWarnings = FALSE)
dir.create(paste0(result_dir,"delta/"), showWarnings = FALSE)
dir.create(paste0(result_dir,"pcf/"), showWarnings = FALSE)

plan(sequential)                                 

initial_population = 5000                           
time_limit = 1000000
min_delta = Inf
epsilon = 1e-6
cell_count_x = 100
params_all<-list()
id=1

for (dd in c(0.02,0.2,0.4)){
  var<-list(
    sm = c(0.96,0.45,0.09),
    sw = c(0.09,0.45,0.96)
  )
  for (t in 1:3){
    params_all[[id]] <-
      data.frame(id=id,
                 sm=var$sm[t],                 
                 sw=var$sw[t],                 
                 b=0.2,
                 d=0.02,dd=dd,
                 start_pop=initial_population,
                 seed=1234)%>%
      mutate(area=pmax(sm, sw) * 1000)
    
    params=params_all[[id]]
    
    start_time = Sys.time()
    params$death_kernel_r = 10 * params$sw
    #params$death_kernel_r = 15
    params$death_kernel_nodes = 1001
    x_grid_death = seq(0,params$death_kernel_r,
                       length.out = params$death_kernel_nodes)
    
    
    params$birth_kernel_r = 10 * params$sm
    #params$birth_kernel_r = 15
    params$birth_kernel_nodes = 1001
    x_grid_birth = seq(0,params$birth_kernel_r,
                       length.out = params$birth_kernel_nodes)
    
    params$area_length_x = params$area
    params$init_density = params$start_pop / params$area_length_x
    
    sim_params <-
      list("area_length_x"=params$area_length_x,
           "cell_count_x"=cell_count_x,
           
           "b"=params$b,
           "d"=params$d,
           "dd"=params$dd,
           
           "seed"=params$seed,
           "init_density"=params$init_density,
           
           "death_kernel_r"=params$death_kernel_r,
           "death_kernel_y"=dnorm(x_grid_death, sd = params$sw),  
           
           "birth_kernel_r"=params$birth_kernel_r,
           "birth_kernel_y"=dnorm(x_grid_birth, sd = params$sm),                  
           
           "spline_precision" = 1e-9
      )
    
    pcf_grid = seq(0,max(c(params$sw,params$sm))*10,length.out = 1001)
    
    sim<-new(poisson_1d,sim_params)
    
    pop<-numeric()
    pcf_estimate<-list()
    pop_sma<-numeric()
    time<-numeric()
    calculated_limit = 0
    delta<-numeric()
    pop_cell<-list()
    
    j = 1
    while(1){
      sim$run_events(sim$total_population)
      time[j]=sim$time
      pop[j]=sim$total_population
      
      cur_pop_cell<-numeric(cell_count_x)
      cur_pop_cell = sim$cell_population
      pop_cell[[j]]<-numeric(cell_count_x)
      
      if(j > 1) {
        for (k in 1:cell_count_x)
          pop_cell[[j]][k] = (pop_cell[[(j-1)]][k] * (j-1) + cur_pop_cell[k])/j
      }
      else {
        pop_cell[[j]] = sim$cell_population
      }
      
      if(j > 1)
        pop_sma[j] = (pop_sma[j-1] * (j-1) + pop[j]) / j
      else
        pop_sma[j] = pop[j]
      
      delta[j] = 0.0
      if(j>1)
      {
        for(k in 1:cell_count_x){
          delta[j] = delta[j] + (pop_cell[[j]][k] - pop_cell[[(j-1)]][k])^2
        }
      }
      points<-unique.ppp(ppp(sim$get_all_coordinates(),
                             rep(0,length(sim$get_all_coordinates())),
                             c(0,sim$area_length_x),
                             c(-sim$area_length_x/2,sim$area_length_x/2)))
      
      K_estimate<-Kest(points,r=pcf_grid,correction="Ripley")
      
      pcf_estimate[[j]]=data.frame(Kest=K_estimate$iso/2,x=pcf_grid)%>%
        mutate(pfc=(Kest-lag(Kest))/(pcf_grid-lag(pcf_grid))/sim$area_length_x)%>%
        pull(pfc)
      
      if(j > 1 && delta[j] < min_delta) {
        min_delta = delta[j]
      }
      if (pop[j] == 0)
      {
        calculated_limit = j
        break
      }
      if(j > 1 && delta[j] < epsilon){
        calculated_limit = j
        break
      }
      if (time[j]>time_limit){
        calculated_limit = j
        break
      }
      j = j + 1
    }
    
    pcf_est_av<-numeric(length(pcf_grid))
    for(j in 1:length(pcf_estimate[[1]])){
      jrow=numeric()
      for (k in 1:calculated_limit){
        jrow[k]=pcf_estimate[[k]][j]
      }
      pcf_est_av[j]=mean(jrow)
    }
    
    pcfs<-data.frame(id=id,r=pcf_grid,y=pcf_est_av)
    pops<-data.frame(id=id,time=time,pop=pop)
    pops_sma<-data.frame(id=id, time=time,pop_sma=pop_sma)
    deltas<-data.frame(id=id,time=time,delta=delta)
    
    dir.create(paste0(result_dir,"pop/", id, "/"), showWarnings = FALSE)
    dir.create(paste0(result_dir,"pcf/", id, "/"), showWarnings = FALSE)
    dir.create(paste0(result_dir,"pop_sma/", id, "/"), showWarnings = FALSE)
    dir.create(paste0(result_dir,"delta/", id, "/"), showWarnings = FALSE)
    
    ggplot(pops %>% select(2:3), aes(x=time,y=pop)) + geom_point(size=0.1) + stat_smooth(size=0.1)
    ggsave(paste0(id,".pdf"), path=paste0(result_dir,"pop/", id, "/"))
    ggplot(pcfs %>% select(2:3), aes(x=r,y=y)) + geom_point(size=0.1) + stat_smooth(size=0.1)
    ggsave(paste0(id,".pdf"), path=paste0(result_dir,"pcf/", id, "/"))
    ggplot(data.frame(time=time,pop_sma=pop_sma), aes(x=time,y=pop_sma)) + geom_point(size=0.1)
    ggsave(paste0(id,".pdf"), path=paste0(result_dir,"pop_sma/", id, "/"))
    ggplot(deltas %>% select(2:3), aes(x=time,y=delta)) + geom_point(size=0.1)
    ggsave(paste0(id,".pdf"), path=paste0(result_dir,"delta/", id, "/"))
    
    write_csv(pops,paste0(result_dir,"pop/",id,"/",id,".csv"))
    write_csv(pcfs,paste0(result_dir,"pcf/",id,"/",id,".csv"))
    write_csv(pops_sma,paste0(result_dir,"pop_sma/",id,"/",id,".csv"))
    write_csv(deltas,paste0(result_dir,"delta/",id,"/",id,".csv"))
    
    id = id + 1
  }
}


params_all<-data.frame(bind_rows(params_all))
write_csv(params_all,paste0(result_dir,"params.csv"))
