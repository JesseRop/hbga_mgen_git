library(SNPassoc)

library(snpStats)
library(dplyr)

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(haven)
library(data.table)
library(pander)
library(HardyWeinberg)
library(kableExtra)
library(stargazer)
library(emmeans)
library(tidyr)

setwd("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/r_analysis/abo_fut2_fut3_malariagen")

rm(list = ls())
mgen_bm_sample = read.table("D:/personal from lt0931/malariagen_analysis/Kenya_GWAS_bact_mal.sample", header = T, stringsAsFactors = F)

##Formating sample file appropriately for plink and R analysis
colnames(mgen_bm_sample)[4] = "SEX"

mgen_bm_sample = mutate(SEX = factor(SEX, levels = c("M","F","D")),
                        cc_status = factor(cc_status, levels = c("CONTROL", "BACT", "D", "MAL")),
                        PLATFORM = factor(PLATFORM),
                        mal_sub_status = factor(mal_sub_status, levels = c("CONTROL", "BACT","BOTH","CM", "D","OTHER","SMA")),
                        bact_sub_status = factor(bact_sub_status, levels = c("CONTROL", "ACI","BHS","D","ECOLI","HIB","MAL","NTS","PNEUMO","STAPH")))

mgen_bm_sample_no_dtypes = mgen_bm_sample[-1,]
mgen_bm_sample_no_dtypes[] = lapply(mgen_bm_sample_no_dtypes, function(x) if(is.factor(x)) factor(x) else x)
mgen_bm_sample_no_dtypes[] = lapply(mgen_bm_sample_no_dtypes, function(x) if(is.factor(x)) as.factor(as.numeric((x))) else x)
mgen_bm_plink = rbind(mgen_bm_sample[1,], mgen_bm_sample_no_dtypes)
mgen_bm_plink[] = lapply(mgen_bm_plink, function(x) if(is.factor(x)) factor(x) else x)

write.table(mgen_bm_plink, "D:/personal from lt0931/malariagen_analysis/Kenya_GWAS_bact_mal_plink.sample", quote = F, row.names = F)

##converting the .sample plink file to a phenotype plink file

mgen_bm_pheno_plink = mgen_bm_plink
colnames(mgen_bm_pheno_plink)[1] = "FID"
colnames(mgen_bm_pheno_plink)[2] = "IID"

mgen_bm_pheno_plink = mgen_bm_pheno_plink[-1,]

mgen_bm_pheno_plink = as.data.frame(lapply(mgen_bm_pheno_plink, function(x) if(is.factor(x)) factor(x) else x))
write.table(mgen_bm_pheno_plink, "D:/personal from lt0931/malariagen_analysis/Kenya_GWAS_bact_mal_pheno_plink.txt", quote = F, row.names = F)

##Creating dummy variables to enable R and plink analysis

mgen_bm_pheno_plink_dummy = as.data.frame(mgen_bm_pheno_plink)

##bactaeremia cases in the major group coded as 2 and controls as 1 with missing as NA under the vector cc_status_bact
mgen_bm_pheno_plink_dummy$cc_status_bact = replace(mgen_bm_pheno_plink_dummy$cc_status, mgen_bm_pheno_plink_dummy$cc_status == 3, NA)

##malaria cases in the major group coded as 2 and controls as 1 with missing as NA under the vector cc_status_mal
mgen_bm_pheno_plink_dummy$cc_status_mal = replace(mgen_bm_pheno_plink_dummy$cc_status, mgen_bm_pheno_plink_dummy$cc_status == 2, NA)
mgen_bm_pheno_plink_dummy$cc_status_mal = replace(mgen_bm_pheno_plink_dummy$cc_status_mal, mgen_bm_pheno_plink_dummy$cc_status_mal == 3, 2)

##severe malaria anemia cases in the malaria sub-group coded as 2 and controls as 1 with the rest as NA under the vector mal_sub_status_sma 
mgen_bm_pheno_plink_dummy$mal_sub_status_sma = replace(mgen_bm_pheno_plink_dummy$mal_sub_status, mgen_bm_pheno_plink_dummy$mal_sub_status %in% c(2,3,4,5), NA)
mgen_bm_pheno_plink_dummy$mal_sub_status_sma = replace(mgen_bm_pheno_plink_dummy$mal_sub_status_sma, mgen_bm_pheno_plink_dummy$mal_sub_status_sma == 6, 2)

##cerebral malaria cases in the malaria sub-group coded as 2 and controls as 1 with the rest as NA under the vector mal_sub_status_cm 
mgen_bm_pheno_plink_dummy$mal_sub_status_cm = replace(mgen_bm_pheno_plink_dummy$mal_sub_status, mgen_bm_pheno_plink_dummy$mal_sub_status %in% c(2,3,5,6), NA)
mgen_bm_pheno_plink_dummy$mal_sub_status_cm = replace(mgen_bm_pheno_plink_dummy$mal_sub_status_cm, mgen_bm_pheno_plink_dummy$mal_sub_status_cm == 4, 2)

##both cerebral malaria and severe malaria cases in the malaria sub-group coded as 2 and controls as 1 with the rest as NA under the vector mal_sub_status_both 
mgen_bm_pheno_plink_dummy$mal_sub_status_both = replace(mgen_bm_pheno_plink_dummy$mal_sub_status, mgen_bm_pheno_plink_dummy$mal_sub_status %in% c(2,4,5,6), NA)
mgen_bm_pheno_plink_dummy$mal_sub_status_both = replace(mgen_bm_pheno_plink_dummy$mal_sub_status_both, mgen_bm_pheno_plink_dummy$mal_sub_status_both == 3, 2)

##other combinations of severe malaria cases in the malaria sub-group coded as 2 and controls as 1 with the rest as NA under the vector mal_sub_status_other 
mgen_bm_pheno_plink_dummy$mal_sub_status_other = replace(mgen_bm_pheno_plink_dummy$mal_sub_status, mgen_bm_pheno_plink_dummy$mal_sub_status %in% c(2,3,4,6), NA)
mgen_bm_pheno_plink_dummy$mal_sub_status_other = replace(mgen_bm_pheno_plink_dummy$mal_sub_status_other, mgen_bm_pheno_plink_dummy$mal_sub_status_other == 5, 2)

##bactaeremia cases cases in the malaria sub-group coded as 2 and controls as 1 with the rest as NA under the vector mal_sub_status_bact 
mgen_bm_pheno_plink_dummy$mal_sub_status_bact = replace(mgen_bm_pheno_plink_dummy$mal_sub_status, mgen_bm_pheno_plink_dummy$mal_sub_status %in% c(3,4,5,6), NA)

##Acinetobacter cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_aci
mgen_bm_pheno_plink_dummy$bact_sub_status_aci = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(3,4,5,6,7,8,9), NA)

##Beta-haemolytic streptococci cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_bhs
mgen_bm_pheno_plink_dummy$bact_sub_status_bhs = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,4,5,6,7,8,9), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_bhs = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_bhs, mgen_bm_pheno_plink_dummy$bact_sub_status_bhs == 3, 2)

##Escherichia coli cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_ecoli
mgen_bm_pheno_plink_dummy$bact_sub_status_ecoli = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,3,5,6,7,8,9), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_ecoli = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_ecoli, mgen_bm_pheno_plink_dummy$bact_sub_status_ecoli == 4, 2)

##Escherichia coli cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_ecoli
mgen_bm_pheno_plink_dummy$bact_sub_status_hib = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,3,4,6,7,8,9), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_hib = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_hib, mgen_bm_pheno_plink_dummy$bact_sub_status_hib == 5, 2)

##Malaria cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_mal
mgen_bm_pheno_plink_dummy$bact_sub_status_mal = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,3,4,5,7,8,9), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_mal = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_mal, mgen_bm_pheno_plink_dummy$bact_sub_status_mal == 6, 2)

##Non typhoidal salmonella cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_nts
mgen_bm_pheno_plink_dummy$bact_sub_status_nts = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,3,4,5,6,8,9), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_nts = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_nts, mgen_bm_pheno_plink_dummy$bact_sub_status_nts == 7, 2)

##Pneumonia cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_pneumo
mgen_bm_pheno_plink_dummy$bact_sub_status_pneumo = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,3,4,5,6,7,9), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_pneumo = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_pneumo, mgen_bm_pheno_plink_dummy$bact_sub_status_pneumo ==8, 2)

##Staph aureus cases in the bactaeremia sub-group coded as 2 and controls as 1 with the rest as NA under the vector bact_sub_status_staph
mgen_bm_pheno_plink_dummy$bact_sub_status_staph = replace(mgen_bm_pheno_plink_dummy$bact_sub_status, mgen_bm_pheno_plink_dummy$bact_sub_status %in% c(2,3,4,5,6,7,8), NA)
mgen_bm_pheno_plink_dummy$bact_sub_status_staph = replace(mgen_bm_pheno_plink_dummy$bact_sub_status_staph, mgen_bm_pheno_plink_dummy$bact_sub_status_staph == 9, 2)

mgen_bm_pheno_plink_dummy[,15:29] = lapply(mgen_bm_pheno_plink_dummy[,15:29], function(x) ifelse(x == 1, 0, 1))

##returning phenotype labels
mgen_bm_pheno_plink_dummy$cc_status = factor(mgen_bm_pheno_plink_dummy$cc_status, labels = c("CONTROL", "BACT", "MAL"))
mgen_bm_pheno_plink_dummy$mal_sub_status = factor(mgen_bm_pheno_plink_dummy$mal_sub_status, labels = c("CONTROL", "BACT","BOTH","CM", "OTHER","SMA"))
mgen_bm_pheno_plink_dummy$bact_sub_status = factor(mgen_bm_pheno_plink_dummy$bact_sub_status, labels = c("CONTROL", "ACI","BHS","ECOLI","HIB","MAL","NTS","PNEUMO","STAPH"))


write.table(mgen_bm_pheno_plink_dummy, "D:/personal from lt0931/malariagen_analysis/Kenya_GWAS_bact_mal_dummy_pheno_plink.txt", quote = F, row.names = F)

##fORMATING mgen exclusion list to have an additional family ID column with similar ID as individual ID. 11/9/2018

mgen_exclusion_list = read.table("D:/personal from lt0931/malariagen_analysis/MGEN_sex_mismatch.excl", header = F, stringsAsFactors = F)
mgen_exclusion_list = mgen_exclusion_list[!(duplicated(mgen_exclusion_list$V1)),]
mgen_exclusion_list_plink = cbind(mgen_exclusion_list, mgen_exclusion_list)
colnames(mgen_exclusion_list_plink) [1] = 'FID'
colnames(mgen_exclusion_list_plink) [2] = 'IID'

write.table(mgen_exclusion_list_plink, "D:/personal from lt0931/malariagen_analysis/MGEN_sex_mismatch_plink", quote = F, row.names = F)

##reading in Sophie's and caro's datasets
su_sm = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/SM_Data.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>"))
su_dset = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/severe_malaria_final_APR2015_formated_to_remove_spaces.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>"))
su_dset = su_dset %>% mutate(serial_scode = coalesce(as.character(serial), as.character(source_code_old)))

cn_dset_raw = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/sm_bact_GWAS_combine06042019.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>"))
length(unique(cn_dset_raw$ID_2)) ##No duplicates for ID_2

cn_dset_raw$Kenyan_ID_raw = cn_dset_raw$Kenyan_ID
cn_dset_raw$Kenyan_ID = gsub("^.._", "", cn_dset_raw$Kenyan_ID)

##merging Jame's sample file with caro's dataset using genotyping IDs
bm_sample_jg_cn = merge(mgen_bm_pheno_plink_dummy, cn_dset_raw[,c(2,15:19)], by.x = "IID", by.y = "ID_2", all = T)

##removing those with missing kenyan IDs before using this as a key for merging
bm_mgen_su_dset = merge(bm_sample_jg_cn ,  su_dset[!is.na(su_dset[, "serial_scode"]),] , by.x = "Kenyan_ID"  ,  by.y = "serial_scode", all.x = T)

##merging in thal genotypes from Gideon
bm_mgen_gn_thal = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/gideon_thal_extra_clinical_data/abo_fut2_genopheno_missing_thal_GN.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>", "NO RESULTS") )
bm_mgen_gn_thal[, "thal_gene"] = ifelse(bm_mgen_gn_thal[, "thal_gene"] == "Hom", "Homo", bm_mgen_gn_thal[, "thal_gene"])
bm_mgen_su_dset = merge(bm_mgen_su_dset, bm_mgen_gn_thal[, c("iid","thal_gene")], by.x = "IID", by.y = "iid", all.x = T)

##merging in thal genotypes from Alex
alex_thal = read_dta("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/alex_thal/malariaGEN participants missing thalassaemia data.dta")
write.csv(alex_thal, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/alex_thal/alex_thal.csv", quote = F, row.names = F)
alex_thal = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/alex_thal/alex_thal.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>"))
bm_mgen_su_dset = merge(bm_mgen_su_dset, alex_thal[,c("kenyan_id", "thal_results")], by.x = "Kenyan_ID", by.y = "kenyan_id", all.x = T)

bm_mgen_su_dset = bm_mgen_su_dset %>% mutate(thall_type = coalesce(as.character(thall_type), as.character(thal_gene), as.character(thal_results)))

bm_mgen_su_dset$thall_type = ifelse(bm_mgen_su_dset$thall_type == "Noresults", NA, bm_mgen_su_dset$thall_type)
bm_mgen_su_dset$thall_type = factor(bm_mgen_su_dset$thall_type, levels = c("Norm", "Het", "Homo"))


#merging and counting strict syndromes
bm_mgen_gn_strict_case = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/gideon_thal_extra_clinical_data/abo_fut2_genopheno_missing_strict_case_GN.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>", "NO RESULTS") )
bm_mgen_su_dset = merge(bm_mgen_su_dset, bm_mgen_gn_strict_case[, c("iid", "bacterium1", "bacterium2", "bacterium3","gram_pos", "strep_pneumo","staph_aureus","grp_a_strep","grp_b_strep", "nts", "gp_other", "h_flu", "h_fluB", "h_flu_other", "acineto", "pseudomonas", "klebs", "gn_other", "anaemia", "resp_distress","resp_irregular", "resp_rate", "resp_deep", "severe_malaria", "u_malaria", "c_malaria", "ad_anaemia","ad_fever","ad_malaria", "dx1", "dx2", "dx3")], by.x = "IID", by.y = "iid", all.x = T)
bm_mgen_su_dset$resp_distress_gn = ifelse(bm_mgen_su_dset$resp_distress == "0", NA, bm_mgen_su_dset$resp_distress)
bm_mgen_su_dset = bm_mgen_su_dset %>% mutate(rd_su_gn = coalesce(rd, resp_distress_gn))
syndromes_count_mal_only = bm_mgen_su_dset %>% count(case, cm, sma, rd, rd_su_gn, resp_distress_gn, sma_not_cm_rd, not_cm_sma_rd, any_cm_sma_rd, cm_sma_rd, cm_not_sma_rd, sma_not_cm_rd, rd_not_cm_sma, cm_and_sma_not_rd, cm_and_rd_not_sma, sma_and_rd_not_cm, not_cm, not_sma, not_rd)

##dataset for James with individuals missing the different syndromes
bm_mgen_su_dset = merge(bm_mgen_su_dset, bm_mgen_su_dset %>% filter(case == "1" ) %>% select(c(2,57,59:65)) %>% distinct(Kenyan_ID, .keep_all = T) %>% gather(exclsv_syndromes, value, -Kenyan_ID) %>% na.omit() %>% select(-value), by = "Kenyan_ID", all.x = T)


write.csv(bm_mgen_su_dset, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_su_dset_including_those_in_excln_list.csv", quote = F, row.names = F)
#write.csv(bm_mgen_su_dset[,c(1,7:37, 52:68, 143, 144,71:74,129:139, 88,89,101,102, 140:142, 112:128 )], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_analysis_dset.csv", quote = F, row.names = F) 
write.csv(bm_mgen_su_dset[,c(1,34,2,6,14,15,32,52, 58, 143, 145,71:74,129:130,134,136, 101, 88,89,140:142, 112:128 )], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_pheno_7831.csv", quote = F, row.names = F) 

##writing out individuals with missing cc_status
write.csv(bm_mgen_su_dset[is.na(bm_mgen_su_dset[,"cc_status"]),c(1,34,2,6,14,15,32,52, 58, 143, 145,71:74,129:130,134,136, 101, 88,89,140:142, 112:128 )], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_pheno_230_missing_cc_status.csv", quote = F, row.names = F) 

##writing out controls for Gideon to issue date of birth/recruitment
write.csv(bm_mgen_su_dset[bm_mgen_su_dset[,"cc_status"] == "CONTROL" & !is.na(bm_mgen_su_dset[,"cc_status"]) ,c(1,2,6,32,52)] %>% distinct(Kenyan_ID, .keep_all = T), "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_pheno_controls_dob.csv", quote = F, row.names = F) 

##writing dataframes with Kenyan IDs with letters
bm_mgen_with_letters_in_kenyan_serial = bm_mgen_su_dset[grepl("[A-Z]", bm_mgen_su_dset$Kenyan_ID_raw), ]
write.csv(bm_mgen_with_letters_in_kenyan_serial, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_dset_with_letters_in_serial.csv", quote = F, row.names = F)

##Individuals in the malariagen JG dataset that are in the exclusion list
individuals_in_malariagen_sample_file_in_the_exclusion_list = bm_mgen_su_dset[bm_mgen_su_dset[,"IID"] %in% mgen_exclusion_list,]
write.csv(individuals_in_malariagen_sample_file_in_the_exclusion_list, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_individuals_in_the_exclusion_list.csv")

##Cases
mgen_pheno_3068 = read.table("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/JJG 18052019/Kenya_MGEN.CP1_only_b37.rs334.sample", stringsAsFactors = F, header = T, na.strings = c("NA", "", ".", "<NA>", "NO RESULTS"))
bm_mgen_su_dset_mgen_cases = merge(bm_mgen_su_dset, mgen_pheno_3068[!is.na(mgen_pheno_3068[,"popn_PC_1"]) & mgen_pheno_3068[,"popn_PC_1"] != "C", c("gtc_id", "nonstrict_caseorcontrol", "cm_sma_nonstrict")], by.x = "IID", by.y = "gtc_id", all = T)
bm_su_mgen_cases = bm_mgen_su_dset_mgen_cases[bm_mgen_su_dset_mgen_cases[,"nonstrict_caseorcontrol"] == "1" & !is.na(bm_mgen_su_dset_mgen_cases[,"nonstrict_caseorcontrol"]),]

bm_su_mgen_cases$exclsv_syndromes_rd_gn = ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"] == "BOTH" & bm_su_mgen_cases[,"resp_distress_gn"] == "1", "cm_sma_rd", ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"] == "CM" & bm_su_mgen_cases[,"resp_distress_gn"] == "1", "cm_and_rd_not_sma", ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"] == "SMA" & bm_su_mgen_cases[,"resp_distress_gn"] == "1", "cm_and_sma_not_rd", ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"] == "OTHER" & bm_su_mgen_cases[,"resp_distress_gn"] == "1", "rd_not_cm_sma", NA ))))
bm_su_mgen_cases = bm_su_mgen_cases %>% mutate(exclsv_syndromes_all = coalesce(exclsv_syndromes, exclsv_syndromes_rd_gn))

bm_su_mgen_cases$cm_mgen = ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"]== "BOTH" | bm_su_mgen_cases[,"cm_sma_nonstrict"]== "CM", "1", NA)
bm_su_mgen_cases$sma_mgen = ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"]== "BOTH" | bm_su_mgen_cases[,"cm_sma_nonstrict"]== "SMA", "1", NA)
bm_su_mgen_cases$other_mgen = ifelse(bm_su_mgen_cases[,"cm_sma_nonstrict"]== "OTHER", "1", NA)
bm_su_mgen_cases$cm_su = bm_su_mgen_cases$cm
bm_su_mgen_cases$sma_su = bm_su_mgen_cases$sma
bm_su_mgen_cases$other_su = bm_su_mgen_cases$not_cm_sma_rd
bm_su_mgen_cases$rd_su = bm_su_mgen_cases$rd
bm_su_mgen_cases$rd_gn = bm_su_mgen_cases$resp_distress_gn

bm_su_mgen_cases$missing_rd  = ifelse(is.na(bm_su_mgen_cases[,"exclsv_syndromes_all"]) , "yes", NA)
write.csv(bm_su_mgen_cases[,c(1,2,147,149:158)], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_su_mgen_cases_with_missing_rd.csv", quote = F, row.names = F)

##removing individuals in the exclusion list and writing out the clean dataset
bm_mgen_su_gwas = bm_mgen_su_dset[which(!(bm_mgen_su_dset[, "IID"] %in% mgen_exclusion_list)),]
#bm_mgen_su_gwas = merge(bm_mgen_su_gwas, bm_mgen_su_gwas %>% filter(case == "1") %>% select(c(2,57,59:65)) %>% gather(exclsv_syndromes, value, -Kenyan_ID) %>% na.omit() %>% select(-value), by = "Kenyan_ID", all.x = T)
write.csv(bm_mgen_su_gwas, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_su_dset_minus_those_in_excln_list.csv", quote = F, row.names = F)
write.csv(bm_mgen_su_gwas[,c(1,34,2,6,14,15,32,52, 58, 143, 145,71:74,129:130,134,136, 101, 88,89,140:142, 112:128 )], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_pheno_5456.csv", quote = F, row.names = F) 

##reading in bactaeremia cases
jg_additional_pheno_request = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/JG 17102019/additional_pheno_data.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>", "NO RESULTS"))
table(jg_additional_pheno_request$kenyan_id %in% bm_mgen_su_dset$Kenyan_ID) 

##reading the plink data into r for analysis
##using read.plink
ABO_FUT2_FUT3_plink_gtypes = read.plink("D:/personal from lt0931/malariagen_analysis/r_analysis/Kenya_GWAS_bact_mal_b37_ABO_FUT2_FUT3", na.strings = c("00", "-9", "NA") )
bm_mgen_gtypes_matrix = cbind(ABO_FUT2_FUT3_plink_gtypes$fam$member, ABO_FUT2_FUT3_plink_gtypes$genotypes@.Data)
colnames(bm_mgen_gtypes_matrix)[1] = "IID"

##merging genotypes to phenotypes
bm_mgen_su_gwas_genopheno = merge(bm_mgen_su_gwas, bm_mgen_gtypes_matrix, by = "IID", all.x = T)

bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno[,c(1:145,529,571,250,220)]

##incorporating sickle genotypes
bm_sickle_and_ethnicities = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/Rotavirus cases/JG controls ethnicity and rs334 data/bacteraemia_for_JR.sample.csv", header = T, stringsAsFactors = F)
bm_mgen_su_gwas_genopheno = merge(bm_mgen_su_gwas_genopheno, bm_sickle_and_ethnicities[,c("ID1","rs334","ethnicity")], by.x = "IID", by.y = "ID1", all.x = T)
bm_mgen_su_gwas_genopheno = merge(bm_mgen_su_gwas_genopheno, mgen_pheno_3068[,c("gtc_id","rsid.rs334.A.add", "compressed_ethnicity")], by.x = "IID", by.y = "gtc_id", all.x = T)
bm_mgen_su_gwas_genopheno$rsid.rs334.A.add = factor(round(as.numeric(bm_mgen_su_gwas_genopheno$rsid.rs334.A.add)), labels = c("AA","AT","TT"))
bm_mgen_su_gwas_genopheno$rs334 = factor(bm_mgen_su_gwas_genopheno$rs334, labels = c("AA","AT","TT"))
bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>% mutate(rs334_all = coalesce(rs334, rsid.rs334.A.add))

bm_mgen_su_gwas_genopheno$ethnicity.x = gsub("Kenya_", "", bm_mgen_su_gwas_genopheno$ethnicity.x)
bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>% mutate(compressed_ethnicity_all = coalesce(ethnicity.x, compressed_ethnicity))

##preparing variables for analysis by converting to proper datatypes and removing unnecessary NAs
bm_mgen_su_gwas_genopheno[,c("rs480133:C", "rs601338:A", "rs8176719:TC", "rs8176746:T")] = lapply(bm_mgen_su_gwas_genopheno[,c("rs480133:C", "rs601338:A", "rs8176719:TC", "rs8176746:T")], function(x) ifelse(x == "00", NA, as.character(x)))

bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg = factor(ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "03" & bm_mgen_su_gwas_genopheno$`rs8176746:T` == "03", "B", ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "03" & bm_mgen_su_gwas_genopheno$`rs8176746:T` == "02", "AB", ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "03" & bm_mgen_su_gwas_genopheno$`rs8176746:T` == "01", "A", ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "02" & bm_mgen_su_gwas_genopheno$`rs8176746:T` == "03", "B" , ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "02" & bm_mgen_su_gwas_genopheno$`rs8176746:T` == "02", "B" , ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "02" & bm_mgen_su_gwas_genopheno$`rs8176746:T` == "01", "A", ifelse(bm_mgen_su_gwas_genopheno$`rs8176719:TC` == "01" , "OO", NA ))))))), levels = c("OO", "A", "B", "AB"))
bm_mgen_su_gwas_genopheno[!is.na(bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg), "non_OvsO"] = factor(ifelse( bm_mgen_su_gwas_genopheno[!is.na(bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg),"phenotypic_bldgrp_4bg"] == "B" | bm_mgen_su_gwas_genopheno[!is.na(bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg),"phenotypic_bldgrp_4bg"] == "A"   | bm_mgen_su_gwas_genopheno[!is.na(bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg), "phenotypic_bldgrp_4bg"] == "AB", "non_O", "OO"), levels = c("OO", "non_O"))

bm_mgen_su_gwas_genopheno$secretor_status = factor(ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03", "nonsecretor", "secretor"), levels = c("secretor", "nonsecretor"))
bm_mgen_su_gwas_genopheno$OvsnonO_secretor_combi = factor(ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` != "03" & bm_mgen_su_gwas_genopheno$non_OvsO == "OO", "O_secretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` != "03" & bm_mgen_su_gwas_genopheno$non_OvsO == "non_O", "nonO_secretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03" & bm_mgen_su_gwas_genopheno$non_OvsO == "OO", "O_nonsecretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03" & bm_mgen_su_gwas_genopheno$non_OvsO == "non_O", "nonO_nonsecretor", NA)))), levels = c("O_secretor", "O_nonsecretor", "nonO_secretor",  "nonO_nonsecretor"))
bm_mgen_su_gwas_genopheno$ABO_secretor_combi = factor(ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` != "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "OO", "O_secretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` != "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "A", "A_secretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` != "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "B", "B_secretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` != "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "AB", "AB_secretor",ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "OO", "O_nonsecretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "A", "A_nonsecretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "B", "B_nonsecretor", ifelse(bm_mgen_su_gwas_genopheno$`rs601338:A` == "03" & bm_mgen_su_gwas_genopheno$phenotypic_bldgrp_4bg == "AB", "AB_nonsecretor", NA)))))))), levels = c("O_secretor", "O_nonsecretor", "A_secretor",  "A_nonsecretor", "B_secretor",  "B_nonsecretor", "AB_secretor",  "AB_nonsecretor"))

##Writing out final file
write.csv(bm_mgen_su_gwas_genopheno, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_su_gwas_genopheno_final_dset_for_analysis.csv", quote = F, row.names = F)

##analysis variables
bm_genopheno = bm_mgen_su_gwas_genopheno[,c(1,2,5,7:30,32,36:40,47,53:68,75,76,83,88,89,90,92,94,95,101,140:142,145:149,154:160)]
bm_genopheno$IID = as.character(bm_genopheno$IID)
bm_genopheno[,c(4:9)] = lapply(bm_genopheno[,c(4:9)], function(x) as.numeric(x))
bm_genopheno[,c(13:28, 35:50, 54:55, 61:75)] = lapply(bm_genopheno[,c(13:28, 35:50, 54:55, 61:75)], function(x) as.factor(x))
bm_genopheno = bm_genopheno[bm_genopheno[,"rs334_all"] != "TT",]

##rs601338
##distribution of secretors in cases and controls
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$secretor_status ))
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$phenotypic_bldgrp_4bg ))
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$non_OvsO ))
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$ABO_secretor_combi ))
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$OvsnonO_secretor_combi ))
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$rs334_all))


prop.table(table(bm_genopheno$bact_sub_status, bm_genopheno$secretor_status), 1)
bm_genopheno %>% filter(!is.na(`rs601338:A`)) %>% count(`rs601338:A`)
bm_genopheno %>% filter(!is.na(`rs601338:A`)) %>% count(`rs601338:A`) %>% with(HWChisq(n))
bm_genopheno %>% filter(!is.na(`rs601338:A`) & cc_status == "CONTROL") %>% count(`rs601338:A`)
bm_genopheno %>% filter(!is.na(`rs601338:A`) & cc_status == "CONTROL") %>% count(`rs601338:A`) %>% with(HWChisq(n))

##rs480133
##distribution of secretors in cases and controls
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$`rs480133:C`))
prop.table(table(bm_genopheno$bact_sub_status, bm_genopheno$`rs480133:C`), 1)
bm_genopheno %>% filter(!is.na(`rs480133:C`)) %>% count(`rs480133:C`)
bm_genopheno %>% filter(!is.na(`rs480133:C`)) %>% count(`rs480133:C`) %>% with(HWChisq(n))
bm_genopheno %>% filter(!is.na(`rs480133:C`) & cc_status == "CONTROL") %>% count(`rs480133:C`)
bm_genopheno %>% filter(!is.na(`rs480133:C`) & cc_status == "CONTROL") %>% count(`rs480133:C`) %>% with(HWChisq(n))

##rs8176719
##distribution of secretors in cases and controls
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$`rs8176719:TC`))
prop.table(table(bm_genopheno$bact_sub_status, bm_genopheno$`rs8176719:TC`), 1)
bm_genopheno %>% filter(!is.na(`rs8176719:TC`)) %>% count(`rs8176719:TC`)
bm_genopheno %>% filter(!is.na(`rs8176719:TC`)) %>% count(`rs8176719:TC`) %>% with(HWChisq(n))
bm_genopheno %>% filter(!is.na(`rs8176719:TC`) & cc_status == "CONTROL") %>% count(`rs8176719:TC`)
bm_genopheno %>% filter(!is.na(`rs8176719:TC`) & cc_status == "CONTROL") %>% count(`rs8176719:TC`) %>% with(HWChisq(n))

##rs8176746
##distribution of secretors in cases and controls
addmargins(table(bm_genopheno$bact_sub_status, bm_genopheno$`rs8176746:T`))
prop.table(table(bm_genopheno$bact_sub_status, bm_genopheno$`rs8176746:T`), 1)
bm_genopheno %>% filter(!is.na(`rs8176746:T`)) %>% count(`rs8176746:T`)
bm_genopheno %>% filter(!is.na(`rs8176746:T`)) %>% count(`rs8176746:T`) %>% with(HWChisq(n))
bm_genopheno %>% filter(!is.na(`rs8176746:T`) & cc_status == "CONTROL") %>% count(`rs8176746:T`)
bm_genopheno %>% filter(!is.na(`rs8176746:T`) & cc_status == "CONTROL") %>% count(`rs8176746:T`) %>% with(HWChisq(n))


prop.table(table(bm_genopheno[!is.na(bm_genopheno[, "rs480133:C"]), c("phenotypic_bldgrp_4bg", "rs480133:C")]), margin = 1)


##converting rs480133:C to a factor
bm_genopheno$`rs480133:C` = factor(bm_genopheno$`rs480133:C`, labels = c("TT", "TC", "CC"))

a_b_ab_oo = c("A", "B", "AB", "OO")
#bact_syndromes = c("cc_status_mal" , "case","cm","cm_not_sma_rd","cm_and_rd_not_sma", "cm_and_sma_not_rd","cm_sma_rd","any_cm_sma_rd",   "not_cm_sma_rd","rd", "rd_not_cm_sma", "sma", "sma_not_cm_rd",  "sma_and_rd_not_cm", "not_cm", "not_sma", "not_rd", "sm_sujg", "cm_sujg", "sma_sujg", "cm_and_sma_sujg" )
##using JG syndromes
#bact_syndromes = c("cc_status_mal" , "mal_sub_status_cm","mal_sub_status_all_cm", "cm_not_sma_rd", "mal_sub_status_sma", "mal_sub_status_all_sma", "mal_sub_status_both", "mal_sub_status_other" )

##writing out final d_set
write.csv(bm_genopheno, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/raw_processed/bm_genopheno_final.csv" , quote = F, row.names = F)

bact_syndromes = c("cc_status_bact", "mal_sub_status_bact", "bact_sub_status_aci", "bact_sub_status_bhs", "bact_sub_status_ecoli", "bact_sub_status_hib", "bact_sub_status_mal", "bact_sub_status_nts", "bact_sub_status_pneumo", "bact_sub_status_staph")

#"cc_status_mal", "mal_sub_status_sma", "mal_sub_status_all_sma", "mal_sub_status_cm", "mal_sub_status_all_cm", "mal_sub_status_both", "mal_sub_status_other",
snpid = c("secretor_status", "`rs480133:C`")

##13052019 Number of cases and controls in individuals included in logistic regression model i.e with complete case, sex, ethnicity, hbb_rs334 and secretor status info
for (i in 1:length(bact_syndromes)) {
  print(bact_syndromes[i])
  print(addmargins(table(bm_genopheno[!is.na(bm_genopheno[,"PC1"]) & !is.na(bm_genopheno[,"rs334_all"]) & !is.na(bm_genopheno[, "phenotypic_bldgrp_4bg"]),bact_syndromes[i]])))
}

for (i in 1:length(bact_syndromes)) {
  print(bact_syndromes[i])
  print(addmargins(table(bm_genopheno[!is.na(bm_genopheno[,"PC1"])  & !is.na(bm_genopheno[, "rs334_all"]) & !is.na(bm_genopheno[, "secretor_status"]) ,bact_syndromes[i]])))
}

for (i in 1:length(bact_syndromes)) {
  print(bact_syndromes[i])
  print(addmargins(table(bm_genopheno[!is.na(bm_genopheno[,"PC1"])  & !is.na(bm_genopheno[, "rs334_all"]) & !is.na(bm_genopheno[, "ABO_secretor_combi"]) ,bact_syndromes[i]])))
}

for (i in 1:length(bact_syndromes)) {
  print(bact_syndromes[i])
  print(addmargins(table(bm_genopheno[!is.na(bm_genopheno[,"PC1"])  & !is.na(bm_genopheno[, "rs334_all"]) & !is.na(bm_genopheno[, "rs480133:C"]) ,bact_syndromes[i]])))
}

for (i in 1:length(bact_syndromes)) {
  print(bact_syndromes[i])
  print(addmargins(table(bm_genopheno[!is.na(bm_genopheno[,"PC1"]) &  !is.na(bm_genopheno[, "rs334_all"]) & !is.na(bm_genopheno[, "rs480133:C"]) & !is.na(bm_genopheno[, "phenotypic_bldgrp_4bg"]) ,bact_syndromes[i]])))
}

##counts of cases individuals with rs334, blood group and not missing either rs601338 or rs480133
sm_sub_syndromes_rs344_abo_su_jcr_4bg_pc_counts_jg = bm_genopheno %>%  filter(!is.na(cc_status_bact) & !is.na(phenotypic_bldgrp_4bg) & !is.na(rs334_all) & !is.na(thall_type) ) %>% count( bact_sub_status, phenotypic_bldgrp_4bg, rs334_all, !is.na(PC1))
write.csv(sm_sub_syndromes_rs344_abo_su_jcr_4bg_pc_counts_jg, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/bact_syndromes_rs344_abo_su_jcr_4bg_pc_counts_jg.csv", quote = F, row.names = F)
sm_sub_syndromes_rs344_abo_nonovso_jg = bm_genopheno %>%  filter(!is.na(cc_status_bact) & !is.na(non_OvsO_su_jcr) & !is.na(rs334_all) & !is.na(thall_type) ) %>% count(bact_sub_status,non_OvsO_su_jcr, rs334_all, !is.na(PC1))
write.csv(sm_sub_syndromes_rs344_abo_nonovso_jg, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/bact_syndromes_rs344_abo_nonovso_jg.csv", quote = F, row.names = F)

#bact_syndromes_count_jg = KG_bact_mal_adtnl_ptypes.ABO_FUT2_FUT3_dosages.sample[which(!(KG_bact_mal_adtnl_ptypes.ABO_FUT2_FUT3_dosages.sample  [, "ID_2"] %in% malariagen_exclusion_list)) ,] %>% filter(!is.na(phenotypic_bldgrp) & !is.na(rs334_all) & cc_status != "D"  ) %>%count(cc_status, mal_sub_status_all_cm, cm, sma, rd, sma_not_cm_rd, not_cm_sma_rd, any_cm_sma_rd, cm_sma_rd, cm_not_sma_rd, sma_not_cm_rd, rd_not_cm_sma, cm_and_sma_not_rd, cm_and_rd_not_sma, sma_and_rd_not_cm, not_cm, not_sma, not_rd)
#write.csv(bact_syndromes_count_jg, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_results/bact_syndromes_count_jg.csv", quote = F, row.names = F)

##checking whether all of the subsyndromes are in the case vector and outsheeting a file with severe malaria cases that are lacking in su dset but occur in jg dset
for (i in 1:length(bact_syndromes)){
  print(table(bm_genopheno[,bact_syndromes[i]]))
  print(table(bm_genopheno[!is.na(bm_genopheno[,bact_syndromes[i]]) & bm_genopheno[,bact_syndromes[i]] == 1 ,bact_syndromes[i]] == bm_genopheno[!is.na(bm_genopheno[,bact_syndromes[i]]) & bm_genopheno[,bact_syndromes[i]] == 1, "cc_status_mal"]))
  
}



##ANALYSIS
##Univariate analysis
sink("possible_covariates_univariate.txt")
posbl_predctr = c( "rs334_all", "thall_type")
for (i in 1:length(bact_syndromes)){
  for (k in 1:length(posbl_predctr)){
    print(paste(">>>>>>>>", bact_syndromes[i], "<<<<<<<<", "/n"))
    print(paste("WITH MALARIA----------------UNIVARIATE"))
    abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  posbl_predctr[k], "+" ,  paste(colnames(bm_genopheno)[c(4:9)], collapse = "+") )),data = bm_genopheno , family = "binomial")
    print(summary(abo_fut_model))
    print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
    print(nobs(abo_fut_model))
    #print(paste("MINUS MALARIA----------------UNIVARIATE"))
    #abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  posbl_predctr[k], "+" ,  paste(colnames(bm_genopheno)[c(4:9)], collapse = "+") )),data = bm_genopheno[bm_genopheno[,"malaria"] !="YES"| is.na(bm_genopheno[,"malaria"]),] , family = "binomial")
    #print(summary(abo_fut_model))
    #print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
    #print(nobs(abo_fut_model))
    
  }
}
sink()

##counts with sickle
counts_with_sickle = bm_genopheno %>% filter(!is.na(cc_status) & !is.na(ethnicity))  %>% count(cc_status, bact_sub_status, sickle)
write.csv(counts_with_sickle, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_sickle.csv", quote = F, row.names = F)

#within everyone
snp_log_rgrsn_everyone = function (predctr) {
  for (i in 1:length(bact_syndromes)){
    
    for (k in 1:length(predctr)){
      print(paste(">>>>>>>>", bact_syndromes[i], "<<<<<<<<"))
      ##including accounting for ethnicity instead of PCs
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k], "+" ,  paste(colnames(bm_genopheno)[c(4:9)], collapse = "+") )),data = bm_genopheno, family = "binomial")
      print(summary(abo_fut_model))
      print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
      print(nobs(abo_fut_model))
      ##10062019 unadjusted
      print(paste("----------------unadjusted----------------"))
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k] )),data = bm_genopheno, family = "binomial")
      print(summary(abo_fut_model))
      print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
      print(nobs(abo_fut_model))
    }
    
  }
}

snp_forest_everyone = function (predctr) {
  for (i in 1:length(bact_syndromes)){
    
    for (k in 1:length(predctr)){
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k], "+" ,  paste(colnames(bm_genopheno)[c(4:9)], collapse = "+") )),data = bm_genopheno, family = "binomial")
      print(plot_model(abo_fut_model, show.values = TRUE, value.offset = .2, value.size = 6, show.p = T, colors = "bw", vline.color = "red",transform = "exp", title = paste(bact_syndromes[i]), rm.terms = c("sickleAT","sickleTT",  "PC_1", "popn_PC_2", "popn_PC_3", "popn_PC_4", "popn_PC_5", "popn_PC_6", "thall_typeHet", "thall_typeHomo"), order.terms =  c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) ) + theme_sjplot(base_size = 18))
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k] )),data = bm_genopheno, family = "binomial")
      print(plot_model(abo_fut_model, show.values = TRUE, value.offset = .2, value.size = 6, show.p = T, colors = "bw", vline.color = "red",transform = "exp", title = paste(bact_syndromes[i], "_u"), order.terms =  c(1,2,3,4,5,6,7,8,9,10,11) ) + theme_sjplot(base_size = 18))
    }
    
  }
}

##OvsnonO
sink("bact_logistic_regression_OvsnonO_su_jcr_everyone.txt")
snp_log_rgrsn_everyone(predctr = "non_OvsO")
sink()

pdf("bact_logistic_regression_OvsnonO_su_jcr_everyone.pdf")
snp_forest_everyone(predctr = "non_OvsO")
dev.off()

##counts with o vs nono
counts_with_O_nonO = bm_genopheno %>% filter(!is.na(cc_status) & !is.na(non_OvsO_jcr))  %>% count(cc_status, P4_exclusive, P6_exclusive, P8_exclusive,non_OvsO_jcr, !is.na(ethnicity) )
write.csv(counts_with_O_nonO, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_O_nonO.csv", quote = F, row.names = F)

##rs601338
#sm_sujg_rs344_secretor_counts = bm_genopheno %>% filter(!is.na(cm_sma_nonstrict) & !is.na(sickle) &!is.na(secretor_status) & !is.na(PC1) & !is.na(thall_type)) %>% count(cm_sma_nonstrict, mal_sub_status_all_cm, cm_not_sma_rd, mal_sub_status_cm, mal_sub_status_all_sma,  mal_sub_status_sma,   mal_sub_status_both, mal_sub_status_other, cm_not_sma_rd)
#write.csv(sm_sujg_rs344_secretor_counts, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/sm_sujg_rs344_secretor_counts.csv", quote = F, row.names = F)

sink("bact_logistic_regression_rs601338_everyone.txt")
snp_log_rgrsn_everyone(predctr = "secretor_status")
sink()

pdf("bact_logistic_regression_rs601338_everyone.pdf")
snp_forest_everyone(predctr = "secretor_status")
dev.off()

counts_with_secretor = bm_genopheno %>% filter(!is.na(cc_status)  & !is.na(secretor_status))  %>% count(cc_status, P4_exclusive, P6_exclusive, P8_exclusive,secretor_status, !is.na(ethnicity))
write.csv(counts_with_secretor, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_secretor.csv", quote = F, row.names = F)

sink("bact_logistic_ABO_secretor_combi.txt")
snp_log_rgrsn_everyone(predctr = "ABO_secretor_combi")
sink()

pdf("bact_logistic_ABO_secretor_combi.pdf")
snp_forest_everyone(predctr = "ABO_secretor_combi")
dev.off()

counts_with_ABO_secretor_combi = bm_genopheno %>% filter(!is.na(cc_status) & !is.na(ABO_secretor_combi) & !is.na(P8_exclusive))  %>% count(P8_exclusive, ABO_secretor_combi, !is.na(ethnicity))
write.csv(counts_with_ABO_secretor_combi, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_ABO_secretor_combi.csv", quote = F, row.names = F)

sink("bact_logistic_OvsnonO_secretor_combi.txt")
snp_log_rgrsn_everyone(predctr = "OvsnonO_secretor_combi")
sink()

pdf("bact_logistic_OvsnonO_secretor_combi.pdf")
snp_forest_everyone(predctr = "OvsnonO_secretor_combi")
dev.off()

counts_with_OvsnonO_secretor_combi = bm_genopheno %>% filter(!is.na(cc_status) & !is.na(P8_exclusive) & !is.na(OvsnonO_secretor_combi))  %>% count( P8_exclusive, OvsnonO_secretor_combi, !is.na(ethnicity))
write.csv(counts_with_OvsnonO_secretor_combi, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_OvsnonO_secretor_combi.csv", quote = F, row.names = F)

sink("bact_logistic_rs480133.txt")
snp_log_rgrsn_everyone(predctr = "`rs480133:C`")
sink()


