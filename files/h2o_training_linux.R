#!/usr/bin/env Rscript

library(h2o)
library('stringr')

#read.delim=","
#temp = list.files(path="C:\\Users\\Mike\\Desktop\\training_data\\input_not_empty_combined",pattern="*.csv")
#myfiles = lapply(temp,read.csv)
#do_called = do.call(rbind, myfiles)

# Get the files names
hhv6_files <- list.files(
  path="/slipstream/home/mmariani/hhv6_ann/training_data/input_not_empty_combined",
  pattern="hhv6",full.names=TRUE)
# First apply read.csv, then rbind
my_hhv6 <- do.call(cbind, lapply(hhv6_files, function(x) read.csv(x, stringsAsFactors = FALSE)))

hg38_files<- list.files(
  path="/slipstream/home/mmariani/hhv6_ann/training_data/input_not_empty_combined",
  pattern="hg38",full.names=TRUE)
my_hg38<- do.call(cbind, lapply(hg38_files, function(x) read.csv(x, stringsAsFactors = FALSE)))

onp_files<- list.files(
  path="/slipstream/home/mmariani/hhv6_ann/training_data/input_not_empty_combined",
  pattern="onp",full.names=TRUE)
my_onp<- do.call(cbind, lapply(onp_files, function(x) read.csv(x, stringsAsFactors = FALSE)))

my_hg38$category=rep("human",times=nrow(my_hg38))
my_hhv6$category=rep("virus",times=nrow(my_hhv6))
my_onp$category=rep("human",times=nrow(my_onp))

training_data<-rbind(my_hg38,my_hhv6)
test_data<-my_onp

#There appears to be at least 1 duplicate column, so let's remove the duplicates:
training_data <- training_data[, !duplicated(colnames(training_data))]
test_data <- test_data[, !duplicated(colnames(test_data))]



#write.csv(training_data, file="C:\\Users\\Mike\\Desktop\\train.csv", row.names=FALSE)
#write.csv(test_data, file="C:\\Users\\Mike\\Desktop\\test.csv", row.names=FALSE)

#train_hex <- h2o.importFile(
#  h2oServer, 
#  path = "C:\\Users\\Mike\\Desktop\\train.csv",
#  #parse_type="CSV", 
#  na.strings=T, 
#  header = F, 
#  parse=TRUE, 
#  sep = ',', 
#  destination_frame = 'train.hex')
#
#test_hex <- h2o.importFile(
#  h2oServer, 
#  parse=TRUE, 
#  path = "C:\\Users\\Mike\\Desktop\\test.csv",
#  #parse_type="CSV", 
#  na.strings=F, 
#  header = F, 
#  sep = ',', 
#  destination_frame = 'test.hex')


setwd("/slipstream/home/mmariani/hhv6_ann/output")

score_test_set=T  #disable if only interested in training throughput

run <- function(extra_params) {
  str(extra_params)
  print("Training.")
  model <- do.call(
    h2o.deeplearning, 
    modifyList(list(x=1:15, y=116,
    training_frame=train_hex, 
    model_id="dlmodel"), 
    extra_params))
  sampleshist <- model@model$scoring_history$samples
  samples <- sampleshist[length(sampleshist)]
  time <- model@model$run_time/1000
  print(paste0("training samples: ", samples))
  print(paste0("training time   : ", time, " seconds"))
  print(paste0("training speed  : ", samples/time, " samples/second"))
  
  if (score_test_set) {
    print("Scoring on test set.")
    ## Note: This scores full test set (10,000 rows) - can take time!
    p <- h2o.performance(model, test_hex)
    cm <- h2o.confusionMatrix(p)
    test_error <- cm$Error[length(cm$Error)]
    print(paste0("test set error  : ", test_error))
  } else {
    test_error <- 1.0
  }
  h2o.rm("dlmodel")
  c(paste(names(extra_params), extra_params, sep = "=", collapse=" "), 
    samples, sprintf("%.3f", time), 
    sprintf("%.3f", samples/time), sprintf("%.3f", test_error))
}

writecsv <- function(results, file, workdir) {
  table <- matrix(unlist(results), ncol = 5, byrow = TRUE)
  colnames(table) <- c("parameters", "training samples",
                       "training time", "training speed", "test set error")
  write.csv(table, 
            file.path(workdir,file), 
            #col.names = T, 
            row.names=F, 
            quote=T 
            #sep=","
            )
}


h2oServer <- h2o.init(ip="localhost", port=54321, max_mem_size="16g", nthreads=-1)

train_hex<-as.h2o(training_data)
test_hex<-as.h2o(test_data)

train_hex[,116] <- as.factor(train_hex[,116])
test_hex[,116] <- as.factor(test_hex[,116])

EPOCHS<-0.1

#args <- list(
#  list(hidden=c(64),             epochs=EPOCHS),
#  list(hidden=c(128),            epochs=EPOCHS),
#  list(hidden=c(256),            epochs=EPOCHS),
#  list(hidden=c(512),            epochs=EPOCHS),
#  list(hidden=c(1024),           epochs=EPOCHS),
#  list(hidden=c(64,64),          epochs=EPOCHS),
#  list(hidden=c(128,128),        epochs=EPOCHS),
#  list(hidden=c(256,256),        epochs=EPOCHS),
#  list(hidden=c(512,512),        epochs=EPOCHS),
#  list(hidden=c(1024,1024),      epochs=EPOCHS),
#  list(hidden=c(64,64,64),       epochs=EPOCHS),
#  list(hidden=c(128,128,128),    epochs=EPOCHS),
#  list(hidden=c(256,256,256),    epochs=EPOCHS),
#  list(hidden=c(512,512,512),    epochs=EPOCHS),
#  list(hidden=c(1024,1024,1024), epochs=EPOCHS)
#)
#writecsv(lapply(args, run), "network_topology.csv")

args <- list(
  list(hidden=c(1024, 1024), epochs=EPOCHS, 
       score_training_samples=100, adaptive_rate=T))

writecsv(lapply(args, run), "ann_slipstream_output.csv", "/slipstream/home/mmariani/hhv6_ann/output")


record_model <- h2o.deeplearning(
  x = 1:115, 
  y = 116, 
  #do_classification=T,  
  training_frame=train_hex, 
  validation_frame = test_hex,
  activation = "RectifierWithDropout", 
  hidden = c(1024,1024,2048),
  epochs = 8000, l1 = 1e-5, 
  input_dropout_ratio = 0.2,
  train_samples_per_iteration = -1, 
  classification_stop = -1)

#Look at model:
#record_model@model$validMetrics$cm

model_path <- h2o.saveModel(object=record_model, path=getwd(), force=TRUE)
