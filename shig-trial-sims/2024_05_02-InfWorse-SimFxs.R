# Basic simulation functions
# basic functions to generate infections
# these functions assume weak infection-conferred immunity
infection_time <- function(parm_start, parm_rate) {
  
  Ti <- floor(((-log(1-runif(1)))+(parm_rate*(parm_start)))/parm_rate) 
  
  return(Ti)
}

infection_severity <- function(parm_progress, parm_severe) {
  
  p_progress <- if_else(parm_progress < 0, 0, parm_progress)
  p_severe <- if_else(parm_severe < 0, 0, parm_severe)
  
  if (rbinom(1,1,p_progress)==1) {
    
    if (rbinom(1,1,p_severe)==1) {severity <- 3
    } else {severity <- 2}
  } else {severity <- 1}
  
  return(severity)
}

modify_hr <- function(parm_VEs, parm_arm, parm_age, parm_infN, parm_infC) {
  
  stopifnot(parm_VEs >= 0, parm_arm %in% c(0,1), parm_age >= 0, 
            parm_infN >= 0, parm_infC %in% c(0,1,2))
  
  update <- 1*((1-parm_VEs)^parm_arm)*parm_age
  
  if(update > 0) {return(update)} else {return(0)}
}

modify_rr_symptoms <- function(parm_VEp, parm_arm, parm_age, parm_infN, parm_infC) {
  
  stopifnot(parm_VEp >= 0, parm_arm %in% c(0,1), parm_age >= 0, 
            parm_infN >= 0, parm_infC %in% c(0,1,2))
  
  cap <- (1-parm_VEp)*((1-0.10)^parm_infC)
  
  if(parm_arm==0) {
    imm <- max(cap, (1-0.10)^parm_infN)
  } else {
    imm <- max(cap, (1-parm_VEp)*((1-0.10)^parm_infN))
  }
  update <- 1*imm*parm_age
  if(update > 0) {return(update)} else {return(0)}
}

modify_rr_severe <- function(parm_VEh, parm_arm, parm_age, parm_infN, parm_infC) {
  
  stopifnot(parm_VEh >= 0, parm_arm %in% c(0,1), parm_age >= 0, 
            parm_infN >= 0, parm_infC %in% c(0,1,2))
  
  cap <- (1-parm_VEh)*((1-0.20)^parm_infC)
  
  if(parm_arm==0) {
    imm <- max(cap, (1-0.20)^parm_infN)
  } else {
    imm <- max(cap, (1-parm_VEh)*((1-0.20)^parm_infN))
  }
  update <- 1*imm*parm_age
  if(update > 0) {return(update)} else {return(0)}
}