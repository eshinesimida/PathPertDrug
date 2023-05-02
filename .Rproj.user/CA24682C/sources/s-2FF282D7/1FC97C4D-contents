#Load data
DE_Colorectal <- get_data1()
ALL_Colorectal <- get_data2()

#get significant pathways
res <- get_pathway_disease(de=DE_Colorectal,all=ALL_Colorectal)

#get drug-induced gene expressionn
drug_exp <- get_drug_data1()
ALL <- get_drug_data2()

#get reverse score
score <- RS(res=res, drug_name = 'metformin', ALL = ALL)
