library(pacman)
p_load(tidyverse,spatstat,future,promises,listenv)
library(MathBioSim)

result_dir = './30_03_2020/'
dir.create(result_dir, showWarnings = FALSE)
dir.create(paste0(result_dir,"pop/"), showWarnings = FALSE)
dir.create(paste0(result_dir,"pcfs/"), showWarnings = FALSE)

plan(sequential)

n_samples = 2000
initial_population = 10000
time_limit = 36000
min_norm = Inf
cell_count_x = 100
epsilon = 1e-6

params_all <-
  data.frame(id=c(1,2,3),
             sm=c(0.84,0.96,0.56),
             sw=c(0.09,0.09,0.13),
             b=1,d=0.1,dd=3/(2*(log(2))),
             samples=n_samples,
             start_pop=initial_population,
             seed=1234)%>%
  mutate(area=pmax(sm, sw) * 1000)

all_runs = listenv()

for (i in 1:nrow(params_all)) {
  params=params_all[i,]
  all_runs[[i]]%<-%
    {
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
      
      pop<-numeric(n_samples)
      time<-numeric(n_samples)
      pcf_estimate<-list()
      calculated_limit = n_samples
      
      for(j in 1:n_samples){
        sim$run_events(sim$total_population)
        pop[j]=sim$total_population
        time[j]=sim$time
        points<-unique.ppp(ppp(sim$get_all_coordinates(),
                               rep(0,sim$total_population),
                               c(0,sim$area_length_x),
                               c(-sim$area_length_x/2,sim$area_length_x/2)))
        
        K_estimate<-Kest(points,r=pcf_grid,correction="Ripley")
        
        pcf_estimate[[j]]=data.frame(Kest=K_estimate$iso/2,x=pcf_grid)%>%
          mutate(pfc=(Kest-lag(Kest))/(pcf_grid-lag(pcf_grid))/sim$area_length_x)%>%
          pull(pfc)
        
        l2_norm = 0.0
        
        populations = sim$get_cell_populations()
        if (j > 1){
          for(k in 1:cell_count_x){
            #cat(j,' ', k, ' ', populations[k], ' ', populations_old[k], '\n')
            l2_norm = l2_norm + (populations[k] - populations_old[k])^2
          }
        }
        populations_old = populations
        
        if (j > 1 && l2_norm < min_norm) {
          min_norm = l2_norm
          cat('parameter id #', i, ', simulation #', j, '; current norm is ', l2_norm, '\n')
        }
        
        if ((j > 1) && (l2_norm < epsilon)) {
          calculated_limit = j
          break
        }
        
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
      
      
      pcfs<-data.frame(id=i,r=pcf_grid,y=pcf_est_av)
      pops<-data.frame(id=i,time=time,pop=pop)
      
      write_csv(pops,paste0(result_dir,"pop/",i,".csv"))
      write_csv(pcfs,paste0(result_dir,"pcfs/",i,".csv"))
    }%stdout%TRUE
}

all_runs%>%as.list()

write_csv(params_all,paste0(result_dir,"params.csv"))

list.files(paste0(result_dir,"pop/"),full.names = TRUE)%>%
  map_dfr(read_csv)%>%
  write_csv(paste0(result_dir,"pop.csv"))

list.files(paste0(result_dir,"pcfs/"),full.names = TRUE)%>%
  map_dfr(read_csv)%>%
  write_csv(paste0(result_dir,"pcfs.csv"))
