#square world functions:
## Last updated: 3-11-16

make_world<-function(rows, columns, ecotype){
  expand.rows<- 1: rows  # expand the number of rows from 1 to number specified
  expand.cols<- 1: columns #expand the number of columns
  
  world<-expand.grid(y = expand.rows, x = expand.cols)  # make all combinations
  
  world$cell_id<- paste("id_", 1:nrow(world), sep="") # add a column of cell ids
  world$cell_prob<- NA # Default a column of probs with NAs
  
  if(ecotype == "sage"){
    world$cell_prob<- runif(nrow(world)) # If random = TRUE replace NAs with random uniform
  }
  if (ecotype == "prairie"){
    world$cell_prob<- runif((nrow(world)), 0.1, 0.3) # If random = TRUE replace NAs with random uniform between 0.1 and 0.3
  }
  if (ecotype == "desert"){
    world$cell_prob<- runif((nrow(world)), 0.7, 0.9) # If random = TRUE replace NAs with random uniform between 0.7 and 0.9
  }
  return(world)  
}

make_world2<-function(rows, columns){
  expand.rows<- 1: rows  # expand the number of rows from 1 to number specified
  expand.cols<- 1: columns #expand the number of columns
  
  world<-expand.grid(y = expand.rows, x = expand.cols)  # make all combinations
  
  world$cell_id<- paste("id_", 1:nrow(world), sep="") # add a column of cell ids
  world$cell_prob<- NA # Default a column of probs with NAs
  return(world)  
}


populate_world<-function(num_organisms, worldToPopulate, num_cell, replace = FALSE){
  require(dplyr)
  potential_orgs<-worldToPopulate[rep(1:nrow(worldToPopulate), each=num_cell),]
  org_locations<-potential_orgs[sample(1:nrow(potential_orgs), num_organisms, replace=replace),]
  org_locations$org<-1
  loc_totals<-org_locations %>%
    group_by(x,y) %>%
    summarise(N=sum(org))
  square_world2<-full_join(worldToPopulate, loc_totals, by=c("x", "y"))
  square_world2$N[is.na(square_world2$N)]<-0
  square_world2
}

run.walk<- function(Nsteps, x1, y1, homerange.size, mu, rho){
  require(circular)
  
  steplength<- rweibull(Nsteps, 2, 350)
  thetaz<- suppressWarnings(rwrappedcauchy(Nsteps, mu= mu, rho= rho))
  uniformz <- runif(Nsteps, 0,1)
  
  walk.valz<- data.frame(step = 1:Nsteps,steplength = steplength, thetaz = thetaz, detect.prob = uniformz, x = NA, y =NA)
  walk.valz <- rbind(data.frame(step = 0,steplength = 0, thetaz = 0, detect.prob = runif(1, 0,1), x = x1, y =y1), walk.valz)
  
  walk.valz <- find.step(walk.valz, homerange.size)
  
  return(walk.valz)
}

find.distance<-function(x1,y1,x2,y2){
  distance<- sqrt((x1 - x2)^2 + (y1 - y2)^2)
  return(distance)
}

find.step<- function(walk.valz,homerange.size){
  walk.valz$new.dist<- NA
  walk.valz$prop.dist <- NA
  
  for(i in 2:nrow(walk.valz)){
    
    dX <- walk.valz$steplength[i]*cos(walk.valz$thetaz[i])
    dY <- walk.valz$steplength[i]*sin(walk.valz$thetaz[i])
    
    potential.step.x <- as.numeric(walk.valz$x[i - 1] + dX)
    potential.step.y <- as.numeric(walk.valz$y[i - 1] + dY)
    
    walk.valz$new.dist[i] <- find.distance(x1 =walk.valz$x[1] ,y1=walk.valz$y[1], potential.step.x, potential.step.y)
    walk.valz$prop.dist[i] <- walk.valz$new.dist[i]/homerange.size
    
    probz <- (0.9999 / (1 + exp(4.741 + -9.407*walk.valz$prop.dist[i])))
    
    if(probz > runif(1,0,1)){
      #  If probability exceeds a random uniform, organism turns back toward origin point
      walk.valz$thetaz[i]<- atan2((walk.valz$y[1]-potential.step.y), (walk.valz$x[1]-potential.step.x))
    }
    
    
    walk.valz$x[i]<-as.numeric(walk.valz$x[i - 1] + walk.valz$steplength[i] * cos(walk.valz$thetaz[i]))
    walk.valz$y[i]<-as.numeric(walk.valz$y[i - 1] + walk.valz$steplength[i] * sin(walk.valz$thetaz[i]))
  }
  
  return(walk.valz[,c("step","detect.prob","x","y")])
}

move_critters<-function(pop_world, myworld, world.type="closed",Nsteps, homerange.type="fixed", homerange.size, mu, rho){
  require(dplyr)
  require(circular)
  
  pop_world.red<-pop_world[pop_world$N>0,] #reduced population to cells with organisms
  
  pop_world.red<-pop_world.red[rep(1:nrow(pop_world.red),pop_world.red$N),c("x","y","cell_id")] # expands the cells for each organism
  
  rownames(pop_world.red)<-NULL
  
  num.org<-nrow(pop_world.red)  #calculate the number of unique organisms
  
  org_stor<-  as.data.frame(matrix(NA,Nsteps*num.org,5))  # create a data.frame to stor organism movement
  names(org_stor)<-c("step","detect.prob","x","y","org_id") # rename column headers
  
  for(i in 1:num.org){
    if(homerange.type=="fixed"){
      run_steps<-run.walk(Nsteps, x1=pop_world.red$x[i], y1=pop_world.red$y[i], homerange.size=homerange.size, mu = mu, rho=rho)  # homerange fixed for all organisms
      
    }
    if(homerange.type=="random"){
      run_steps<-run.walk(Nsteps, x1=pop_world.red$x[i], y1=pop_world.red$y[i], homerange.size=rpois(1,homerange.size), mu = mu, rho=rho) # mean homerange size is equal to homerange size (poisson distribution)
      
    }
    
    run_steps$org_id<-paste("org",i,sep="_")
    
    org_stor[(Nsteps*(i-1)+i):((Nsteps*i +i)),]<-run_steps
    
  }
  
  if(world.type=="closed"){
    min.y<- min(myworld$y)   #  Organism stays on edge, can not leave myworld
    max.y<- max(myworld$y)
    min.x<- min(myworld$x)
    max.x<- max(myworld$x)
    
    org_stor$y[org_stor$y<min.y]<-min.y
    org_stor$y[org_stor$y>max.y]<-max.y
    org_stor$x[org_stor$x<min.x]<-min.x
    org_stor$x[org_stor$x>max.x]<-max.x
    
  }
  
  org_stor$cell_x<-floor(org_stor$x) # round down to find cell number in x and y
  org_stor$cell_y<-floor(org_stor$y)
  
  org_stor<-left_join(org_stor,myworld, by=c("cell_x"="x", "cell_y"="y"))  # join to myworld to get cell specific detection
  
  if(world.type=="open"){  
    org_stor$cell_prob[is.na(org_stor$cell_prob)]<-0 # sets cell specific detection to 0 if organism moves off myworld
  }
  
  org_stor<-org_stor %>% 
    mutate(do.detect = ifelse(detect.prob<cell_prob,1,0))  # create a binary value for detection. 1 = detect
  
  return(org_stor)
  
}


set_ecotype <- function(hr.world, ecotype){
  
  uniq_cells<-unique(hr.world$cell_id)
  eco_detect<-data.frame(cell_id = uniq_cells, cell_prob=NA)
  
  if(ecotype == "sage"){
    eco_detect$cell_prob<- runif(nrow(eco_detect)) # If random = TRUE replace NAs with random uniform
  }
  if (ecotype == "prairie"){
    eco_detect$cell_prob<- runif((nrow(eco_detect)), 0.1, 0.3) # If random = TRUE replace NAs with random uniform between 0.1 and 0.3
  }
  if (ecotype == "desert"){
    eco_detect$cell_prob<- runif((nrow(eco_detect)), 0.7, 0.9) # If random = TRUE replace NAs with random uniform between 0.7 and 0.9
  }
  
  head(hr.world)
  
  hr.world$cell_prob<- eco_detect$cell_prob[match(hr.world$cell_id,eco_detect$cell_id)]
  
  
  hr.world<-hr.world %>% 
    mutate(do.detect = ifelse(detect.prob<cell_prob,1,0))  # create a binary value for detection. 1 = detect
  
  return(hr.world)
  
  }



sample_myworld<-function(world,n.samples,size.x, size.y, max.iter = 1000,seed=NULL){
  #creates a function to generate random sampling grids in myworld 
  if(!is.null(seed)){
    set.seed(seed) # if seed is not null, set the seed to seed
  }
  
  s = 1 # counter for while loop
  t = 1 # counter for total iterations
  
  available.cells<-world$cell_id # create object of available points. In this case is all the rows that we have initially
  
  sample_stor<-data.frame(matrix(NA,(size.x*size.y)*n.samples,5)) #creates storage for sample cells
  
  names(sample_stor)<-c(names(world),"sample_id") # rename columns
  
  while(s<=n.samples){
    start.point<-world[sample(1:nrow(world), size=1, replace = F),] # sample an initial row
    
    if(start.point$cell_id %in% available.cells){ 
      sample.dir<-data.frame(UD=sample(c('U','D'), size=1), LR=sample(c('L','R'), size=1))
      
      if(sample.dir$UD=='U'& sample.dir$LR=='R'){
        poss.x<-seq(start.point$x, start.point$x + (size.x-1), by=1) 
        poss.y<-seq(start.point$y, start.point$y + (size.y-1), by=1) 
      }
      
      if(sample.dir$UD=='U'& sample.dir$LR=='L'){
        poss.x<-seq(start.point$x, start.point$x - (size.x-1), by=-1) 
        poss.y<-seq(start.point$y, start.point$y + (size.y-1), by=1) 
      }
      
      if(sample.dir$UD=='D'& sample.dir$LR=='R'){
        poss.x<-seq(start.point$x, start.point$x + (size.x-1), by=1) 
        poss.y<-seq(start.point$y, start.point$y - (size.y-1), by=-1) 
      }
      
      if(sample.dir$UD=='D'& sample.dir$LR=='L'){
        poss.x<-seq(start.point$x, start.point$x - (size.x-1), by=-1) 
        poss.y<-seq(start.point$y, start.point$y - (size.y-1), by=-1) 
      }
      
      poss.cells<-expand.grid(x=poss.x, y=poss.y) 
      poss.cells<-left_join(poss.cells, world, by=c("x","y")) 
      
      if(all(poss.cells$cell_id %in% available.cells)){
        sample.cells<-poss.cells
        sample.cells$sample_id<-s
        
        sample_stor[(((s-1)*(size.x*size.y))+1): (s*size.x*size.y),]<- sample.cells
        
        # start ((s-1)*(size.x*size.y))+1
        # end s*size.x*size.y
        
        available.cells<-available.cells[-which(available.cells %in% sample.cells$cell_id)]
        
        s<-s+1
      } 
    }    
    t = t+1
    if(t == max.iter){
      stop("You have reach the maximun number of available iterations")
    }
    
  }
  
  return(sample_stor)
}

datez_needed<-function(Nsteps,j, by_sample_val){
  
  if(j <= 6){
    dates_needed<-seq(1,Nsteps,by=by_sample_val)
  }
  if(j >6){
    
    x<-as.matrix(seq(1, Nsteps, by = 12))
    x1<-cbind(x[,rep(1,by_sample_val)])
    vals<- t(as.matrix(0:(by_sample_val-1)))[rep(1,nrow(x1)),]
    
    dates_needed<-sort(as.numeric(unlist(x1+vals)))
  }
  return(dates_needed)
}

try_loco_1 <- function (xy, arange,percent=percent, unin , unout) {
  out <- tryCatch(LoCoH.a.area(xy, arange, percent=percent, unin , unout), error = function(e) NA)
  return(out)
}

try_locoh <-function(coords,max_dist_all, percent = 90,unin, unout){
  # 
  # try_out<- try_loco_1(xy, arange = 2*max_dist_all,percent=percent, unin , unout)
  
  results<-list()
  
  possible_a_down<- seq(1.9,1,by=-0.1) 
  possible_a_up<- seq(2.1,3,by=0.1)
  
  poss_values <- data.frame(vals = c(2,possible_a_down,possible_a_up), order=c(1,seq(from=2,by=2,length.out = length(possible_a_down)),seq(from= 3,by=2,length.out = length(possible_a_up))))
  
  poss_values <-poss_values$vals[order(poss_values$order)]
  
  stop <- 0
  s = 1
  while(stop ==0){
    # print(s)
    try_out<- try_loco_1(xy= coords, arange = poss_values[s]*max_dist_all,percent=percent, unin , unout)
    s = s+1
    
    if(!is.na(try_out)){stop = 1; aval = poss_values[s-1]}
    if(s>length(poss_values)){stop = 1; aval = NA}
  }
  
    if(!is.na(aval)){
      results[["HR_calculate"]] <- 1
      results[["hr"]] <- as.numeric(try_out)
      results[["a"]] <- aval
      
    }else{
      results[["HR_calculate"]] <- 0 
      results[["hr"]] <- NA
      results[["a"]] <- NA
    }
    
  return(results)
}

