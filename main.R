folder <- getwd()
functions <- list.files(folder)
functions <- grep("main.R", functions, invert = T, value = T)
functions <- paste(folder, functions, sep = "/")
for (i in functions){
  source(i)
}

###########################################################################################################
######################################### Heterogeneous Effect ############################################
###########################################################################################################
data.cont.cont            <- makeData.cont.eff.cont.paper(N = 500, n.test = 1000)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:400, ]
data.validation.cont.cont <- data.used.full.cont.cont[401:500, ]

###########################################################################################################
######################################### CV 1 Node Specific ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.estns.cv1.cont.cont <- create.sequence(data.used     = data.used.cont.cont, 
                                                   etemp.used    = etemp.ns, 
                                                   stemp.used    = stemp.ns, 
                                                   itemp.used    = itemp,
                                                   need.cond.exp = F, 
                                                   type.var      = "cont",
                                                   use.var       = "true")

# seq.created.estns.cv1.cont.cont$tree.list
final.tree.estns.cv1.cont.cont <- EstNs.CvMethod1(data.used   = data.used.cont.cont, 
                                                  tree.list   = seq.created.estns.cv1.cont.cont$tree.list, 
                                                  lambda.used = 4, 
                                                  val.sample  = data.validation.cont.cont,
                                                  type.var    = "cont")
# final.tree.estns.cv1.cont.cont
t1 <- Sys.time()
eval.final.estns.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.estns.cv1.cont.cont[[1]], 
                                                    test.data    = data.cont.cont$test.data,
                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                    noise.var    = data.cont.cont$noise.var,
                                                    corr.split   = data.cont.cont$corr.split)
eval.final.estns.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 1 Estimator 1 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est1.cv1.cont.cont <- create.sequence(data.used     = data.used.cont.cont, 
                                                  etemp.used    = etemp.ms, 
                                                  stemp.used    = stemp.ms,
                                                  itemp.used    = itemp,
                                                  need.cond.exp = F, 
                                                  type.var      = "cont",
                                                  use.var       = "true")

# seq.created.est1.cv1.cont.cont$tree.list
final.tree.est1.cv1.cont.cont <- Est1.CvMethod1(data.used   = data.used.cont.cont, 
                                                tree.list   = seq.created.est1.cv1.cont.cont$tree.list, 
                                                lambda.used = 4, 
                                                val.sample  = data.validation.cont.cont,
                                                type.var    = "cont",
                                                use.var     = "true")
#final.tree.est1.cv1.cont.cont
t1 <- Sys.time()
eval.final.est1.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.est1.cv1.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est1.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
######################################## CV 1 True Estimator 1 ############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est1.cv1.cont.cont <- create.sequence.TruePaper(data.used     = data.used.cont.cont, 
                                                                 etemp.used    = etemp.ms.TruePaper, 
                                                                 stemp.used    = stemp.ms.TruePaper, 
                                                                 itemp.used    = itemp, 
                                                                 need.cond.exp = F, 
                                                                 type.var      = "cont", 
                                                                 use.var       = "true",
                                                                 eff           = T)

# seq.created.est1.cv1.cont.cont$tree.list
final.tree.true.est1.cv1.cont.cont <- Est1.CvMethod1.TruePaper(data.used   = data.used.cont.cont, 
                                                               tree.list   = seq.created.true.est1.cv1.cont.cont$tree.list, 
                                                               lambda.used = 4, 
                                                               val.sample  = data.validation.cont.cont,
                                                               type.var    = "cont",
                                                               use.var     = "true",
                                                               eff         = T)
#final.tree.true.est1.cv1.cont.cont
t1 <- Sys.time()
eval.final.true.est1.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est1.cv1.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est1.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 1 Estimator 2 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est2.cv1.cont.cont <- create.sequence(data.used     = data.used.cont.cont, 
                                                  etemp.used    = etemp.da,
                                                  stemp.used    = stemp.da, 
                                                  itemp.used    = itemp,
                                                  need.cond.exp = T,
                                                  type.var      = "cont", 
                                                  cond.exp.used = "GAM",
                                                  use.var       = "true")

# seq.created.est2.cv1.cont.cont$tree.list
final.tree.est2.cv1.cont.cont <- Est2.CvMethod1(data.used     = data.used.cont.cont, 
                                                tree.list     = seq.created.est2.cv1.cont.cont$tree.list, 
                                                lambda.used   = 4, 
                                                val.sample    = data.validation.cont.cont,
                                                type.var      = "cont", 
                                                cond.exp.used = "GAM",
                                                use.var       = "true")
# final.tree.est2.cv1.cont.cont
t1 <- Sys.time()
eval.final.est2.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.est2.cv1.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est2.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
######################################### CV 1 True Estimator 2 ###########################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est2.cv1.cont.cont <- create.sequence.TruePaper(data.used     = data.used.cont.cont, 
                                                                 etemp.used    = etemp.da, 
                                                                 stemp.used    = stemp.da,
                                                                 itemp.used    = itemp,
                                                                 need.cond.exp = T, 
                                                                 type.var      = "cont", 
                                                                 cond.exp.used = "GAM",
                                                                 use.var       = "true",
                                                                 eff           = T)

# seq.created.true.est2.cv1.cont.cont$tree.list
final.tree.true.est2.cv1.cont.cont <- Est2.CvMethod1.TruePaper(data.used     = data.used.cont.cont, 
                                                               tree.list     = seq.created.true.est2.cv1.cont.cont$tree.list, 
                                                               lambda.used   = 4, 
                                                               val.sample    = data.validation.cont.cont,
                                                               type.var      = "cont", 
                                                               cond.exp.used = "GAM",
                                                               use.var       = "true",
                                                               eff           = T)
# final.tree.true.est2.cv1.cont.cont
t1 <- Sys.time()
eval.final.true.est2.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est2.cv1.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est2.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 2 Node Specific ############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.estns.cv2.cont.cont <- create.sequence(data.used     = data.used.full.cont.cont, 
                                                   etemp.used    = etemp.ns, 
                                                   stemp.used    =  stemp.ns, 
                                                   itemp.used    =  itemp,
                                                   need.cond.exp = F, 
                                                   type.var      = "cont",
                                                   use.var       = "true")

# seq.created.estns.cv2.cont.cont$tree.list
final.tree.estns.cv2.cont.cont <- EstNs.CvMethod2(data.used = data.used.full.cont.cont, 
                                                  tree.list = seq.created.estns.cv2.cont.cont$tree.list, 
                                                  type.var  = "cont")
# final.tree.estns.cv2.cont.cont
t1 <- Sys.time()
eval.final.estns.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.estns.cv2.cont.cont[[1]], 
                                                    test.data    = data.cont.cont$test.data,
                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                    noise.var    = data.cont.cont$noise.var,
                                                    corr.split   = data.cont.cont$corr.split)
eval.final.estns.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 2 Estimator 1 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est1.cv2.cont.cont <- create.sequence(data.used     = data.used.full.cont.cont, 
                                                  etemp.used    = etemp.ms, 
                                                  stemp.used    = stemp.ms, 
                                                  itemp.used    = itemp,
                                                  need.cond.exp = F, 
                                                  type.var      = "cont",
                                                  use.var       = "true")

# seq.created.estns.cv2.cont.cont$tree.list
final.tree.est1.cv2.cont.cont <- Est1.CvMethod2(data.used = data.used.full.cont.cont, 
                                                tree.list = seq.created.est1.cv2.cont.cont$tree.list, 
                                                type.var  = "cont")
# final.tree.est1.cv2.cont.cont
t1 <- Sys.time()
eval.final.est1.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.est1.cv2.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est1.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
######################################### CV 2 True Estimator 1 ###########################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est1.cv2.cont.cont <- create.sequence.TruePaper(data.used     =  data.used.full.cont.cont, 
                                                                 etemp.used    =  etemp.ms.TruePaper, 
                                                                 stemp.used    =  stemp.ms.TruePaper, 
                                                                 itemp.used    =  itemp,
                                                                 need.cond.exp = F, 
                                                                 type.var      = "cont", 
                                                                 use.var       = "true",
                                                                 eff           = T)

# seq.created.true.estns.cv2.cont.cont$tree.list
final.tree.true.est1.cv2.cont.cont <- Est1.CvMethod2.TruePaper(data.used = data.used.full.cont.cont, 
                                                               tree.list = seq.created.true.est1.cv2.cont.cont$tree.list, 
                                                               type.var  = "cont",
                                                               eff       = T)
# final.tree.true.est1.cv2.cont.cont
t1 <- Sys.time()
eval.final.true.est1.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est1.cv2.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est1.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 2 Estimator 2 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est2.cv2.cont.cont <- create.sequence(data.used     = data.used.full.cont.cont, 
                                                  etemp.used    = etemp.da, 
                                                  stemp.used    = stemp.da, 
                                                  itemp.used    = itemp,
                                                  need.cond.exp = T, 
                                                  type.var      = "cont", 
                                                  cond.exp.used = "GAM",
                                                  use.var       = "true")
# seq.created.est2.cv2.cont.cont$tree.list
tree.list.cont.cont   <- seq.created.est2.cv2.cont.cont$tree.list
final.tree.est2.cv2.cont.cont <- Est2.CvMethod2(data.used = data.used.full.cont.cont, 
                                                tree.list = tree.list.cont.cont, 
                                                type.var = "cont",
                                                cond.exp.used = "GAM")
# final.tree.est2.cv2.cont.cont
t1 <- Sys.time()
eval.final.est2.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.est2.cv2.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est2.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################## CV 2 True Estimator 2 ##########################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est2.cv2.cont.cont <- create.sequence.TruePaper(data.used     = data.used.full.cont.cont, 
                                                                 etemp.used    = etemp.da, 
                                                                 stemp.used    = stemp.da, 
                                                                 itemp.used    = itemp,
                                                                 need.cond.exp = T, 
                                                                 type.var      = "cont", 
                                                                 cond.exp.used = "GAM", 
                                                                 use.var       = "true",
                                                                 eff           = T)
# seq.created.true.est2.cv2.cont.cont$tree.list
tree.list.true.cont.cont   <- seq.created.true.est2.cv2.cont.cont$tree.list
final.tree.true.est2.cv2.cont.cont <- Est2.CvMethod2.TruePaper(data.used     = data.used.full.cont.cont, 
                                                               tree.list     = tree.list.true.cont.cont,
                                                               type.var      = "cont",
                                                               cond.exp.used = "GAM",
                                                               eff           = T)
# final.tree.true.est2.cv2.cont.cont
t1 <- Sys.time()
eval.final.true.est2.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est2.cv2.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est2.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
################################################## MOB ####################################################
###########################################################################################################
t0 <- Sys.time()
fit.mob = mob.fit(data.used = data.used.full.cont.cont)
t1 <- Sys.time()
eval.mob = eval.measures.eff.mob(final.tree   = fit.mob, 
                                 test.data    = data.cont.cont$test.data, 
                                 true.trt.eff = data.cont.cont$true.trt.eff, 
                                 data.used    = data.used.full.cont.cont,
                                 noise.var    = data.cont.cont$noise.var,
                                 corr.split   = data.cont.cont$corr.split)
eval.mob$t <- as.numeric(difftime(t1, t0, units = "secs"))

t0 <- Sys.time()
fit.mob.2 = mob.fit.2(data.used = data.used.full.cont.cont)
t1 <- Sys.time()
eval.mob.2 = eval.measures.eff.mob(final.tree   = fit.mob.2, 
                                   test.data    = data.cont.cont$test.data, 
                                   true.trt.eff = data.cont.cont$true.trt.eff, 
                                   data.used    = data.used.full.cont.cont,
                                   noise.var    = data.cont.cont$noise.var,
                                   corr.split   = data.cont.cont$corr.split)
eval.mob.2$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
################################################# VT ######################################################
###########################################################################################################
t0 <- Sys.time()
vt.tree = vt.sim.cont(data.used = data.used.full.cont.cont)
t1 <- Sys.time()
eval.vt = eval.measures.eff(final.tree   = vt.tree$vt.tree, 
                            test.data    = data.cont.cont$test.data,
                            true.trt.eff = data.cont.cont$true.trt.eff,
                            noise.var    = data.cont.cont$noise.var,
                            corr.split   = data.cont.cont$corr.split)
eval.vt$t <- as.numeric(difftime(t1, t0, units = "secs"))

performance.eff <- list(ns.cv1        = eval.final.estns.cv1.cont.cont,
                        est1.cv1      = eval.final.est1.cv1.cont.cont,
                        true.est1.cv1 = eval.final.true.est1.cv1.cont.cont,
                        est2.cv1      = eval.final.est2.cv1.cont.cont,
                        true.est2.cv1 = eval.final.true.est2.cv1.cont.cont,
                        ns.cv2        = eval.final.estns.cv2.cont.cont,
                        est1.cv2      = eval.final.est1.cv2.cont.cont,
                        true.est1.cv2 = eval.final.true.est1.cv2.cont.cont,
                        est2.cv2      = eval.final.est2.cv2.cont.cont,
                        true.est2.cv2 = eval.final.true.est2.cv2.cont.cont,
                        mob           = eval.mob,
                        mob.2         = eval.mob.2, 
                        vt            = eval.vt)



###########################################################################################################
########################################### Homogeneous Effect ############################################
###########################################################################################################
data.cont.cont            <- makeData.cont.noeff.cont.paper(N = 500, n.test = 1000)
data.used.full.cont.cont  <- data.cont.cont$data.used
data.used.cont.cont       <- data.used.full.cont.cont[1:400, ]
data.validation.cont.cont <- data.used.full.cont.cont[401:500, ]

###########################################################################################################
######################################### CV 1 Node Specific ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.estns.cv1.cont.cont <- create.sequence(data.used     = data.used.cont.cont, 
                                                   etemp.used    = etemp.ns, 
                                                   stemp.used    = stemp.ns, 
                                                   itemp.used    = itemp,
                                                   need.cond.exp = F, 
                                                   type.var      = "cont",
                                                   use.var       = "true")

# seq.created.estns.cv1.cont.cont$tree.list
final.tree.estns.cv1.cont.cont <- EstNs.CvMethod1(data.used   = data.used.cont.cont, 
                                                  tree.list   = seq.created.estns.cv1.cont.cont$tree.list, 
                                                  lambda.used = 4, 
                                                  val.sample  = data.validation.cont.cont,
                                                  type.var    = "cont")
# final.tree.estns.cv1.cont.cont
t1 <- Sys.time()
eval.final.estns.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.estns.cv1.cont.cont[[1]], 
                                                    test.data    = data.cont.cont$test.data,
                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                    noise.var    = data.cont.cont$noise.var,
                                                    corr.split   = data.cont.cont$corr.split)
eval.final.estns.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 1 Estimator 1 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est1.cv1.cont.cont <- create.sequence(data.used     = data.used.cont.cont, 
                                                  etemp.used    = etemp.ms, 
                                                  stemp.used    = stemp.ms,
                                                  itemp.used    = itemp,
                                                  need.cond.exp = F, 
                                                  type.var      = "cont",
                                                  use.var       = "true")

# seq.created.est1.cv1.cont.cont$tree.list
final.tree.est1.cv1.cont.cont <- Est1.CvMethod1(data.used   = data.used.cont.cont, 
                                                tree.list   = seq.created.est1.cv1.cont.cont$tree.list, 
                                                lambda.used = 4, 
                                                val.sample  = data.validation.cont.cont,
                                                type.var    = "cont",
                                                use.var     = "true")
#final.tree.est1.cv1.cont.cont
t1 <- Sys.time()
eval.final.est1.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.est1.cv1.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est1.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
######################################## CV 1 True Estimator 1 ############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est1.cv1.cont.cont <- create.sequence.TruePaper(data.used     = data.used.cont.cont, 
                                                                 etemp.used    = etemp.ms.TruePaper, 
                                                                 stemp.used    = stemp.ms.TruePaper, 
                                                                 itemp.used    = itemp, 
                                                                 need.cond.exp = F, 
                                                                 type.var      = "cont", 
                                                                 use.var       = "true",
                                                                 eff           = F)

# seq.created.est1.cv1.cont.cont$tree.list
final.tree.true.est1.cv1.cont.cont <- Est1.CvMethod1.TruePaper(data.used   = data.used.cont.cont, 
                                                               tree.list   = seq.created.true.est1.cv1.cont.cont$tree.list, 
                                                               lambda.used = 4, 
                                                               val.sample  = data.validation.cont.cont,
                                                               type.var    = "cont",
                                                               use.var     = "true",
                                                               eff         = F)
#final.tree.true.est1.cv1.cont.cont
t1 <- Sys.time()
eval.final.true.est1.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est1.cv1.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est1.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 1 Estimator 2 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est2.cv1.cont.cont <- create.sequence(data.used     = data.used.cont.cont, 
                                                  etemp.used    = etemp.da,
                                                  stemp.used    = stemp.da, 
                                                  itemp.used    = itemp,
                                                  need.cond.exp = T,
                                                  type.var      = "cont", 
                                                  cond.exp.used = "GAM",
                                                  use.var       = "true")

# seq.created.est2.cv1.cont.cont$tree.list
final.tree.est2.cv1.cont.cont <- Est2.CvMethod1(data.used     = data.used.cont.cont, 
                                                tree.list     = seq.created.est2.cv1.cont.cont$tree.list, 
                                                lambda.used   = 4, 
                                                val.sample    = data.validation.cont.cont,
                                                type.var      = "cont", 
                                                cond.exp.used = "GAM",
                                                use.var       = "true")
# final.tree.est2.cv1.cont.cont
t1 <- Sys.time()
eval.final.est2.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.est2.cv1.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est2.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
######################################### CV 1 True Estimator 2 ###########################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est2.cv1.cont.cont <- create.sequence.TruePaper(data.used     = data.used.cont.cont, 
                                                                 etemp.used    = etemp.da, 
                                                                 stemp.used    = stemp.da,
                                                                 itemp.used    = itemp,
                                                                 need.cond.exp = T, 
                                                                 type.var      = "cont", 
                                                                 cond.exp.used = "GAM",
                                                                 use.var       = "true",
                                                                 eff           = F)

# seq.created.true.est2.cv1.cont.cont$tree.list
final.tree.true.est2.cv1.cont.cont <- Est2.CvMethod1.TruePaper(data.used     = data.used.cont.cont, 
                                                               tree.list     = seq.created.true.est2.cv1.cont.cont$tree.list, 
                                                               lambda.used   = 4, 
                                                               val.sample    = data.validation.cont.cont,
                                                               type.var      = "cont", 
                                                               cond.exp.used = "GAM",
                                                               use.var       = "true",
                                                               eff           = F)
# final.tree.true.est2.cv1.cont.cont
t1 <- Sys.time()
eval.final.true.est2.cv1.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est2.cv1.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est2.cv1.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 2 Node Specific ############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.estns.cv2.cont.cont <- create.sequence(data.used     = data.used.full.cont.cont, 
                                                   etemp.used    = etemp.ns, 
                                                   stemp.used    =  stemp.ns, 
                                                   itemp.used    =  itemp,
                                                   need.cond.exp = F, 
                                                   type.var      = "cont",
                                                   use.var       = "true")

# seq.created.estns.cv2.cont.cont$tree.list
final.tree.estns.cv2.cont.cont <- EstNs.CvMethod2(data.used = data.used.full.cont.cont, 
                                                  tree.list = seq.created.estns.cv2.cont.cont$tree.list, 
                                                  type.var  = "cont")
# final.tree.estns.cv2.cont.cont
t1 <- Sys.time()
eval.final.estns.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.estns.cv2.cont.cont[[1]], 
                                                    test.data    = data.cont.cont$test.data,
                                                    true.trt.eff = data.cont.cont$true.trt.eff,
                                                    noise.var    = data.cont.cont$noise.var,
                                                    corr.split   = data.cont.cont$corr.split)
eval.final.estns.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 2 Estimator 1 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est1.cv2.cont.cont <- create.sequence(data.used     = data.used.full.cont.cont, 
                                                  etemp.used    = etemp.ms, 
                                                  stemp.used    = stemp.ms, 
                                                  itemp.used    = itemp,
                                                  need.cond.exp = F, 
                                                  type.var      = "cont",
                                                  use.var       = "true")

# seq.created.estns.cv2.cont.cont$tree.list
final.tree.est1.cv2.cont.cont <- Est1.CvMethod2(data.used = data.used.full.cont.cont, 
                                                tree.list = seq.created.est1.cv2.cont.cont$tree.list, 
                                                type.var  = "cont")
# final.tree.est1.cv2.cont.cont
t1 <- Sys.time()
eval.final.est1.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.est1.cv2.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est1.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
######################################### CV 2 True Estimator 1 ###########################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est1.cv2.cont.cont <- create.sequence.TruePaper(data.used     =  data.used.full.cont.cont, 
                                                                 etemp.used    =  etemp.ms.TruePaper, 
                                                                 stemp.used    =  stemp.ms.TruePaper, 
                                                                 itemp.used    =  itemp,
                                                                 need.cond.exp = F, 
                                                                 type.var      = "cont", 
                                                                 use.var       = "true",
                                                                 eff           = F)

# seq.created.true.estns.cv2.cont.cont$tree.list
final.tree.true.est1.cv2.cont.cont <- Est1.CvMethod2.TruePaper(data.used = data.used.full.cont.cont, 
                                                               tree.list = seq.created.true.est1.cv2.cont.cont$tree.list, 
                                                               type.var  = "cont",
                                                               eff       = F)
# final.tree.true.est1.cv2.cont.cont
t1 <- Sys.time()
eval.final.true.est1.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est1.cv2.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est1.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################### CV 2 Estimator 2 ##############################################
###########################################################################################################
t0 <- Sys.time()
seq.created.est2.cv2.cont.cont <- create.sequence(data.used     = data.used.full.cont.cont, 
                                                  etemp.used    = etemp.da, 
                                                  stemp.used    = stemp.da, 
                                                  itemp.used    = itemp,
                                                  need.cond.exp = T, 
                                                  type.var      = "cont", 
                                                  cond.exp.used = "GAM",
                                                  use.var       = "true")
# seq.created.est2.cv2.cont.cont$tree.list
tree.list.cont.cont   <- seq.created.est2.cv2.cont.cont$tree.list
final.tree.est2.cv2.cont.cont <- Est2.CvMethod2(data.used = data.used.full.cont.cont, 
                                                tree.list = tree.list.cont.cont, 
                                                type.var = "cont",
                                                cond.exp.used = "GAM")
# final.tree.est2.cv2.cont.cont
t1 <- Sys.time()
eval.final.est2.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.est2.cv2.cont.cont[[1]], 
                                                   test.data    = data.cont.cont$test.data,
                                                   true.trt.eff = data.cont.cont$true.trt.eff,
                                                   noise.var    = data.cont.cont$noise.var,
                                                   corr.split   = data.cont.cont$corr.split)
eval.final.est2.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
########################################## CV 2 True Estimator 2 ##########################################
###########################################################################################################
t0 <- Sys.time()
seq.created.true.est2.cv2.cont.cont <- create.sequence.TruePaper(data.used     = data.used.full.cont.cont, 
                                                                 etemp.used    = etemp.da, 
                                                                 stemp.used    = stemp.da, 
                                                                 itemp.used    = itemp,
                                                                 need.cond.exp = T, 
                                                                 type.var      = "cont", 
                                                                 cond.exp.used = "GAM", 
                                                                 use.var       = "true",
                                                                 eff           = F)
# seq.created.true.est2.cv2.cont.cont$tree.list
tree.list.true.cont.cont   <- seq.created.true.est2.cv2.cont.cont$tree.list
final.tree.true.est2.cv2.cont.cont <- Est2.CvMethod2.TruePaper(data.used     = data.used.full.cont.cont, 
                                                               tree.list     = tree.list.true.cont.cont,
                                                               type.var      = "cont",
                                                               cond.exp.used = "GAM",
                                                               eff           = F)
# final.tree.true.est2.cv2.cont.cont
t1 <- Sys.time()
eval.final.true.est2.cv2.cont.cont <- eval.measures.eff(final.tree   = final.tree.true.est2.cv2.cont.cont[[1]], 
                                                        test.data    = data.cont.cont$test.data,
                                                        true.trt.eff = data.cont.cont$true.trt.eff,
                                                        noise.var    = data.cont.cont$noise.var,
                                                        corr.split   = data.cont.cont$corr.split)
eval.final.true.est2.cv2.cont.cont$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
################################################## MOB ####################################################
###########################################################################################################
t0 <- Sys.time()
fit.mob = mob.fit(data.used = data.used.full.cont.cont)
t1 <- Sys.time()
eval.mob = eval.measures.eff.mob(final.tree   = fit.mob, 
                                 test.data    = data.cont.cont$test.data, 
                                 true.trt.eff = data.cont.cont$true.trt.eff, 
                                 data.used    = data.used.full.cont.cont,
                                 noise.var    = data.cont.cont$noise.var,
                                 corr.split   = data.cont.cont$corr.split)
eval.mob$t <- as.numeric(difftime(t1, t0, units = "secs"))

t0 <- Sys.time()
fit.mob.2 = mob.fit.2(data.used = data.used.full.cont.cont)
t1 <- Sys.time()
eval.mob.2 = eval.measures.eff.mob(final.tree   = fit.mob.2, 
                                   test.data    = data.cont.cont$test.data, 
                                   true.trt.eff = data.cont.cont$true.trt.eff, 
                                   data.used    = data.used.full.cont.cont,
                                   noise.var    = data.cont.cont$noise.var,
                                   corr.split   = data.cont.cont$corr.split)
eval.mob.2$t <- as.numeric(difftime(t1, t0, units = "secs"))

###########################################################################################################
################################################# VT ######################################################
###########################################################################################################
t0 <- Sys.time()
vt.tree = vt.sim.cont(data.used = data.used.full.cont.cont)
t1 <- Sys.time()
eval.vt = eval.measures.eff(final.tree   = vt.tree$vt.tree, 
                            test.data    = data.cont.cont$test.data,
                            true.trt.eff = data.cont.cont$true.trt.eff,
                            noise.var    = data.cont.cont$noise.var,
                            corr.split   = data.cont.cont$corr.split)
eval.vt$t <- as.numeric(difftime(t1, t0, units = "secs"))

performance.noeff <- list(ns.cv1        = eval.final.estns.cv1.cont.cont,
                          est1.cv1      = eval.final.est1.cv1.cont.cont,
                          true.est1.cv1 = eval.final.true.est1.cv1.cont.cont,
                          est2.cv1      = eval.final.est2.cv1.cont.cont,
                          true.est2.cv1 = eval.final.true.est2.cv1.cont.cont,
                          ns.cv2        = eval.final.estns.cv2.cont.cont,
                          est1.cv2      = eval.final.est1.cv2.cont.cont,
                          true.est1.cv2 = eval.final.true.est1.cv2.cont.cont,
                          est2.cv2      = eval.final.est2.cv2.cont.cont,
                          true.est2.cv2 = eval.final.true.est2.cv2.cont.cont,
                          mob           = eval.mob,
                          mob.2         = eval.mob.2, 
                          vt            = eval.vt)