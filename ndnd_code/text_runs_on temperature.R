# small test script to run 3 simulations with 3 temperatures

load("~/Documents/Work/NDND/NDND_in_R/NDND/ndnd_data/NDNDData_2019_04_09_10_22_33.RData") # my reference data file
save(NDNDData,file = './NDNDData.Rdata')
source('./ndnd_code/NDND_main.r')
NDNDData=TEOMR(NDNDData,c(1,1,1,1,1,1,0,0),1)
NDNDData$Data.tag=Sys.time()
NDNDData$comment="temperature increased by 1 degree from original file"
save(NDNDData,file = './NDNDData.Rdata')
source('./ndnd_code/NDND_main.r')
NDNDData=TEOMR(NDNDData,c(1,1,1,1,1,1,0,0),1)
NDNDData$Data.tag=Sys.time()
NDNDData$comment="temperature increased by 2 degree from original file"
save(NDNDData,file = './NDNDData.Rdata')
source('./ndnd_code/NDND_main.r')
