#!/usr/bin/env Rscript

library("h2o")
onp_files<- list.files(
  path="/slipstream/home/mmariani/hhv6_ann/training_data/input_not_empty_combined",
  pattern="onp",full.names=TRUE)
my_onp<- do.call(cbind, lapply(onp_files, function(x) read.csv(x, stringsAsFactors = FALSE)))

my_onp$category=rep("human",times=nrow(my_onp))

test_data<-my_onp

#There appears to be at least 1 duplicate column, so let's remove the duplicates:
test_data <- test_data[, !duplicated(colnames(test_data))]

h2oServer <- h2o.init(ip="localhost", port=54321, max_mem_size="16g", nthreads=-1)

test_hex<-as.h2o(test_data)

test_hex[,116] <- as.factor(test_hex[,116])

model_path<-"/slipstream/home/mmariani/hhv6_ann/output/DeepLearning_model_R_1535305455367_4"
saved_model <- h2o.loadModel(model_path)
predictions <- h2o.predict(saved_model, test_hex)

write.csv(as.data.frame(predictions), 
          file="/slipstream/home/mmariani/hhv6_ann/predictions.csv",
          row.names=FALSE)
