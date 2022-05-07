#!/usr/bin/env Rscript

#source("https://bioconductor.org/biocLite.R")
#biocLite("ShortRead")
library("rDNAse")
require(rDNAse)
library("ShortRead") #for converting fastq ONP to fasta

#onp_fasta = readFASTA(system.file(file.choose(),
#                               package = 'rDNAse'))
#neg_hs = readFASTA(system.file('dnaseq/non-hs.fasta',
#                               package = 'rDNAse'))

#onp_fasta<-readFASTA(file="/slipstream/home/mmariani/hhv6_ann/reads/basecalled.fasta")
#onp_fasta<-readFASTA(file="/slipstream/home/mmariani/hhv6_ann/reads/sim_telomere.fa")
#onp_fasta<-readFASTA(file="/slipstream/home/mmariani/hhv6_ann/reads/sim_hhv6.fa")

#onp_fasta<-writeFasta(
#  readFastq(
#    "/slipstream/home/mmariani/hhv6_ann/reads", 
#    "pass_merged.fastq"), 
#  "/slipstream/home/mmariani/hhv6_ann/reads/round_2_pass_merged.fasta"
#  )
onp_fasta<-readFASTA(file="/slipstream/home/mmariani/hhv6_ann/reads/round_2_pass_merged.fasta")

x1<-t(sapply(onp_fasta,kmer)) #- Basic kmer and Reverse compliment kmer
x2<-t(sapply(onp_fasta,make_idkmer_vec)) #- Increment of diversity (ID)
#Autocorrelation
x3<-t(sapply(onp_fasta,extrDAC)) #- Dinucleotide-based auto covariance
x4<-t(sapply(onp_fasta,extrDCC)) #- Dinucleotide-based cross covariance
x5<-t(sapply(onp_fasta,extrDACC)) #- Dinucleotide-based auto-cross covariance
x6<-t(sapply(onp_fasta,extrTAC)) #- Trinucleotide-based auto covariance
x7<-t(sapply(onp_fasta,extrTCC)) #- Trinucleotide-based cross covariance
x8<-t(sapply(onp_fasta,extrDACC)) #- Trinucleotide-based auto-cross covariance
#Pseudo nucleotide composition
x9<-t(sapply(onp_fasta,extrPseDNC)) #- Pseudo dinucleotide composition
x10<-t(sapply(onp_fasta,extrPseKNC)) #- Pseudo k-tupler nucleotide composition
x12<-t(sapply(onp_fasta,twoSeqSim))
x13<-t(sapply(onp_fasta,twoGOSim)) 

#output_dir<-"C:\\Users\\Mike\\Desktop\\training_data\\onp"
#output_dir<-"C:\\Users\\Mike\\Desktop\\training_data\\telomere"
#output_dir<-"C:\\Users\\Mike\\Desktop\\training_data\\hhv6"
output_dir<-"/slipstream/home/mmariani/hhv6_ann/output/hhv6_onp_2"

write.csv(x1, file=paste0(output_dir,"/","x1",".csv"), row.names = FALSE)
write.csv(x2, file=paste0(output_dir,"/","x2",".csv"), row.names = FALSE)
write.csv(x3, file=paste0(output_dir,"/","x3",".csv"), row.names = FALSE)
write.csv(x4, file=paste0(output_dir,"/","x4",".csv"), row.names = FALSE)
write.csv(x5, file=paste0(output_dir,"/","x5",".csv"), row.names = FALSE)
write.csv(x6, file=paste0(output_dir,"/","x6",".csv"), row.names = FALSE)
write.csv(x7, file=paste0(output_dir,"/","x7",".csv"), row.names = FALSE)
write.csv(x8, file=paste0(output_dir,"/","x8",".csv"), row.names = FALSE)
write.csv(x9, file=paste0(output_dir,"/","x9",".csv"), row.names = FALSE)
write.csv(x10, file=paste0(output_dir,"/","x10",".csv"), row.names = FALSE)
write.csv(x11, file=paste0(output_dir,"/","x11",".csv"), row.names = FALSE)
write.csv(x12, file=paste0(output_dir,"/","x12",".csv"), row.names = FALSE)
write.csv(x13, file=paste0(output_dir,"/","x13",".csv"), row.names = FALSE)
