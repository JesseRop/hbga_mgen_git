library(SNPassoc)

library(snpStats)
library(tidyverse)

library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(ggplot2)
library(data.table)
library(pander)
library(HardyWeinberg)
library(kableExtra)
library(stargazer)
library(emmeans)
library(haven)
library(lubridate)
library(gtsummary)
library(flextable)
library(officer)
library(labelled)

setwd("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/r_analysis/abo_fut2_fut3_malariagen")

rm(list = ls())
mgen_bm_sample = read.table("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/Kenya_GWAS_bact_mal.sample", header = T, stringsAsFactors = F)

##Formating sample file appropriately for plink and R analysis
colnames(mgen_bm_sample)[4] = "SEX"

mgen_bm_sample = mgen_bm_sample %>% mutate(
  SEX = factor(SEX, levels = c("M","F","D")),
  cc_status = factor(cc_status, levels = c("CONTROL", "BACT", "D", "MAL")),
  PLATFORM = factor(PLATFORM),
  mal_sub_status = factor(mal_sub_status, levels = c("CONTROL", "BACT","BOTH","CM", "D","OTHER","SMA")),
  bact_sub_status = factor(bact_sub_status, levels = c("CONTROL", "ACI","BHS","D","ECOLI","HIB","MAL","NTS","PNEUMO","STAPH")))

##Writing a sample file for snptest analysis input
mgen_bm_sample[-1,] %>%
  lapply(., function(x) if(is.factor(x)) as.factor(as.numeric((x))) else x) %>% 
  as.data.frame() %>%
  rbind(mgen_bm_sample[1,], .) %>%
  #lapply(., function(x) if(is.factor(x)) factor(x) else x) %>%
  write.table(., "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/Kenya_GWAS_bact_mal_plink.sample", quote = F, row.names = F)

##Writing a sample file for plink analysis input
mgen_bm_sample[-1,] %>%
  lapply(., function(x) if(is.factor(x)) as.factor(as.numeric((x))) else x) %>% 
  as.data.frame() %>%
  rename("FID" = "ID_1", "IID" = "ID_2") %>%
  write.table(., "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/Kenya_GWAS_bact_mal_pheno_plink.txt", quote = F, row.names = F)

##Reading in mgen exclusion list and removing duplicates
mgen_exclusion_list = read.table("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/MGEN_sex_mismatch.excl", header = F, stringsAsFactors = F) %>% 
  filter(!(duplicated(V1)))

##Mgen exclusion list for plink
mgen_exclusion_list_plink = mgen_exclusion_list %>%
  mutate(V2 = V1) %>%
  rename(FID = V1, IID = V2) %>% 
  write.table(., "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/MGEN_sex_mismatch_plink", quote = F, row.names = F)


##Creating dummy variables to enable R and plink analysis
mgen_bm_pheno_plink_dummy = mgen_bm_sample[-1,] %>%
  as.data.frame() %>%
  mutate(
    cc_status_bact = case_when(cc_status == "CONTROL" ~ 0, cc_status == "BACT" ~ 1),
    cc_status_mal = case_when(cc_status == "CONTROL" ~ 0, cc_status == "MAL" ~ 1),
    mal_sub_status_sma = case_when(mal_sub_status == "CONTROL" ~ 0, mal_sub_status == "SMA" ~ 1),
    mal_sub_status_cm = case_when(mal_sub_status == "CONTROL" ~ 0, mal_sub_status == "CM" ~ 1),
    mal_sub_status_both = case_when(mal_sub_status == "CONTROL" ~ 0, mal_sub_status == "BOTH" ~ 1),
    mal_sub_status_other = case_when(mal_sub_status == "CONTROL" ~ 0, mal_sub_status == "OTHER" ~ 1),
    mal_sub_status_bact = case_when(mal_sub_status == "CONTROL" ~ 0, mal_sub_status == "BACT" ~ 1),
    bact_sub_status_aci = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "ACI" ~ 1),
    bact_sub_status_bhs = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "BHS" ~ 1),
    bact_sub_status_ecoli = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "ECOLI" ~ 1),
    bact_sub_status_hib = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "HIB" ~ 1),
    bact_sub_status_mal = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "MAL" ~ 1),
    bact_sub_status_nts = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "NTS" ~ 1),
    bact_sub_status_pneumo = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "PNEUMO" ~ 1),
    bact_sub_status_staph = case_when(bact_sub_status == "CONTROL" ~ 0, bact_sub_status == "STAPH" ~ 1)
  )

##Writing out mgen file with dummy variables
write.table(mgen_bm_pheno_plink_dummy, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/Kenya_GWAS_bact_mal_dummy_pheno_plink.txt", quote = F, row.names = F)


##reading in Sophie's and caro's datasets
su_sm = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/SM_Data.csv", na = c("NA", "", ".", "<NA>"))

su_dset = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/severe_malaria_final_APR2015_formated_to_remove_spaces.csv", na = c("NA", "", ".", "<NA>")) %>% 
  mutate(serial_scode = coalesce(as.character(serial), as.character(source_code_old)))

cn_dset_raw = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/sm_bact_GWAS_combine06042019.csv", na = c("NA", "", ".", "<NA>"), col_types = cols(Kenyan_ID = col_character())) %>% 
  mutate(Kenyan_ID_raw = Kenyan_ID,
         Kenyan_ID = gsub("^.._", "", Kenyan_ID))

length(unique(cn_dset_raw$ID_2)) ##No duplicates for ID_2

##merging Jame's sample file with caro's and Sophie's dataset (after removing those with missing Kenyan IDs) using genotyping IDs
bm_mgen_su_dset = cn_dset_raw %>% 
  select(c("ID_2", "Oxford_ID", "Kenyan_ID", "Study", "missing_kenyanid", "Kenyan_ID_raw")) %>%
  left_join(mgen_bm_pheno_plink_dummy, ., by = "ID_2") %>% 
  rename(IID = ID_2) %>%
  left_join(., su_dset %>% filter(!is.na(serial_scode)), by = c("Kenyan_ID" = "serial_scode" ))

##merging in thal genotypes from Gideon
bm_mgen_su_dset = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/gideon_thal_extra_clinical_data/abo_fut2_genopheno_missing_thal_GN.csv", na = c("NA", "", ".", "<NA>", "NO RESULTS")) %>%
  mutate(thal_gene = case_when(thal_gene == "Hom" ~ "Homo", TRUE ~ thal_gene)) %>%
  select(c("iid","thal_gene")) %>%
  left_join(bm_mgen_su_dset, ., by = c("IID" = "iid"))


##merging in thal genotypes from Alex
bm_mgen_su_dset = read_dta("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/alex_thal/malariaGEN participants missing thalassaemia data.dta") %>%
  mutate_all(~na_if(.,  "")) %>%
  mutate_all(~na_if(.,  "NA")) %>% 
  select(c("kenyan_id", "thal_results")) %>%
  mutate(kenyan_id = as.character(kenyan_id)) %>%
  left_join(bm_mgen_su_dset, ., by = c("Kenyan_ID" = "kenyan_id"))

##Coalescing and formatting Thalassaemia genotype 
bm_mgen_su_dset = bm_mgen_su_dset %>% 
  mutate(thall_type = coalesce(as.character(thall_type), as.character(thal_gene), as.character(thal_results)),
         thall_type = case_when(thall_type == "Homo" ~ "Hom", TRUE ~ as.character(thall_type)),
         thall_type = na_if(thall_type, "Noresults"),
         thall_type = factor(thall_type, levels = c("Norm", "Het", "Hom"))
  )

##Coalescing function from https://alistaire.rbind.io/blog/coalescing-joins/
coalesce_join <- function(x, y, 
                          by = NULL, suffix = c(".x", ".y"), 
                          join = dplyr::full_join, ...) {
  joined <- join(x, y, by = by, suffix = suffix, ...)
  # names of desired output
  cols <- union(names(x), names(y))
  
  to_coalesce <- names(joined)[!names(joined) %in% cols]
  suffix_used <- suffix[ifelse(endsWith(to_coalesce, suffix[1]), 1, 2)]
  # remove suffixes and deduplicate
  to_coalesce <- unique(substr(
    to_coalesce, 
    1, 
    nchar(to_coalesce) - nchar(suffix_used)
  ))
  
  coalesced <- purrr::map_dfc(to_coalesce, ~dplyr::coalesce(
    joined[[paste0(.x, suffix[1])]], 
    joined[[paste0(.x, suffix[2])]]
  ))
  names(coalesced) <- to_coalesce
  
  dplyr::bind_cols(joined, coalesced)[cols]
}

##formatting the data types appropriately bfore merging
bm_mgen_su_dset = bm_mgen_su_dset %>%
  mutate(dod = dmy(dod),
         dob = dmy(dob),
         doa = dmy(doa),
         dset1 = "Yes")

##Reading in control dates from gideon and coalescing to main dataset
bm_mgen_su_dset = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/Rotavirus cases/Gideon extra data for cases and controls/rota_controls_extra_data_gn.csv", na =c("NA", "", ".", "<NA>")) %>% 
  select(kenyan_id, date_entry, dob) %>%
  mutate(doa = dmy(date_entry),
         dob = dmy(dob),
         Kenyan_ID = as.character(kenyan_id)) %>%
  select(-c(date_entry, kenyan_id)) %>%
  coalesce_join(bm_mgen_su_dset, ., by = 'Kenyan_ID', join = dplyr::left_join)

#merging and counting strict syndromes
bm_mgen_su_dset = read_dta("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/gideon_thal_extra_clinical_data/1995_2020_severe_malaria.dta")  %>%
  mutate(serialno = as.character(serialno),
         dss = as.character(dss),
         sex = as.character(sex),
         ethnic = as.character(ethnic),
         dset2 = "Yes") %>%
  mutate_if(is.character, ~na_if(.,"")) %>%
  rename(Kenyan_ID = serialno) %>%
  coalesce_join(bm_mgen_su_dset, ., by = 'Kenyan_ID', join = dplyr::left_join)

##dataset for James with individuals missing the different syndromes
bm_mgen_su_dset = bm_mgen_su_dset %>% 
  mutate (not_cm_sma_rd = case_when(cc_status != "CONTROL" & is.na(cm) & is.na(sma) & is.na(rd) ~ "1", cc_status == "CONTROL" ~ "0"), 
          cm_not_sma_rd = case_when(cm == "1" & is.na(sma) & is.na(rd) ~ "1", cc_status == "CONTROL" ~ "0"), 
          sma_not_cm_rd = case_when(sma == "1" & is.na(cm) & is.na(rd) ~ "1", cc_status == "CONTROL" ~ "0"), 
          cm_and_sma_not_rd = case_when(cm == "1" & sma == "1" & is.na(rd) ~ "1", cc_status == "CONTROL" ~ "0"), 
          cm_and_rd_not_sma = case_when(cm == "1" & rd == "1" & is.na(sma) ~ "1", cc_status == "CONTROL" ~ "0"), 
          sma_and_rd_not_cm = case_when(sma == "1" & rd == "1" & is.na(cm) ~ "1", cc_status == "CONTROL" ~ "0"),
          rd_not_cm_sma = case_when(rd == "1" & is.na(cm) & is.na(sma)  ~ "1", cc_status == "CONTROL" ~ "0"),
          any_cm_sma_rd = case_when(rd == "1" | cm == "1" | sma == "1" ~ "1", cc_status == "CONTROL" ~ "0")) %>%
  mutate_at(vars(not_cm_sma_rd, cm_not_sma_rd, sma_not_cm_rd, cm_and_sma_not_rd, cm_and_rd_not_sma, sma_and_rd_not_cm, rd_not_cm_sma, any_cm_sma_rd), ~factor(., labels = c("CONTROL", "CASE")))
  
bm_mgen_su_dset = bm_mgen_su_dset %>%
  filter(cc_status != "CONTROL") %>%
  select(c("Kenyan_ID", "not_cm_sma_rd", "cm_sma_rd", "cm_not_sma_rd", "sma_not_cm_rd", "rd_not_cm_sma", "cm_and_sma_not_rd", "cm_and_rd_not_sma", "sma_and_rd_not_cm")) %>% 
  distinct(Kenyan_ID, .keep_all = T) %>% 
  gather(exclsv_syndromes, value, -Kenyan_ID) %>% 
  na.omit() %>% 
  select(-value) %>% 
  left_join(bm_mgen_su_dset, ., by = "Kenyan_ID")


write.csv(bm_mgen_su_dset[,c("Kenyan_ID", "Kenyan_ID_raw", "IID", "cc_status", "mal_sub_status","bact_sub_status", "Study", "case", "any_cm_sma_rd", "exclsv_syndromes", "prostration", "severe_malaria_anaemia", "cerebral_malaria", "respiratory_distress", "thall_type", "diag1", "diag2")], "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_pheno_7831.csv", quote = F, row.names = F) 

##writing dataframes with Kenyan IDs with letters
bm_mgen_su_dset[grepl("[A-Z]", bm_mgen_su_dset$Kenyan_ID_raw), ] %>% write_csv(., "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_dset_with_letters_in_serial.csv")

##Individuals in the malariagen JG dataset that are in the exclusion list
bm_mgen_su_dset[bm_mgen_su_dset[,"IID"] %in% mgen_exclusion_list,] %>% write_csv(., "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_individuals_in_the_exclusion_list.csv")

##Reading 3068 mgen participants and merging with mgen_su 
bm_su_mgen_cases = read.table("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/JJG 18052019/Kenya_MGEN.CP1_only_b37.rs334.sample", stringsAsFactors = F, header = T, na.strings = c("NA", "", ".", "<NA>", "NO RESULTS")) %>% 
  filter(popn_PC_1 != "C") %>% 
  select(c("gtc_id", "caseorcontrol", "cm_sma_nonstrict","rsid.rs334.A.add", "compressed_ethnicity"), contains("popn_")) %>% 
  left_join(bm_mgen_su_dset, ., by = c("IID" = "gtc_id")) #%>% 
  # filter(nonstrict_caseorcontrol == "1" & !is.na(nonstrict_caseorcontrol))

##Formatting malaria syndromes appropriately
bm_su_mgen_cases = bm_su_mgen_cases %>%
  mutate(exclsv_syndromes_rd = case_when(cm_sma_nonstrict == "BOTH" & rd == "1" ~ "cm_sma_rd", 
                                            cm_sma_nonstrict == "CM" & rd == "1" ~ "cm_and_rd_not_sma", 
                                            cm_sma_nonstrict == "SMA" & rd == "1" ~ "sma_and_rd_not_cm", 
                                            cm_sma_nonstrict == "OTHER" & rd == "1" ~ "rd_not_cm_sma"),
         exclsv_syndromes_all = coalesce(exclsv_syndromes, exclsv_syndromes_rd),
         cm_mgen = case_when(cm_sma_nonstrict== "BOTH" | cm_sma_nonstrict== "CM" ~ "1"),
         sma_mgen = case_when(cm_sma_nonstrict== "BOTH" | cm_sma_nonstrict== "SMA" ~ "1"),
         other_mgen = case_when(cm_sma_nonstrict== "OTHER" ~ "1"), 
         missing_rd  = case_when(is.na(exclsv_syndromes_all) ~ "yes")
  )


##removing individuals in the exclusion list and writing out the clean dataset
# bm_mgen_su_gwas = bm_mgen_su_dset[which(!(bm_mgen_su_dset[, "IID"] %in% mgen_exclusion_list)),]
# 

##reading in bactaeremia cases
jg_additional_pheno_request = read.csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/JG 17102019/additional_pheno_data.csv", stringsAsFactors = F, na.strings = c("NA", "", ".", "<NA>", "NO RESULTS"))
table(jg_additional_pheno_request$kenyan_id %in% bm_mgen_su_dset$Kenyan_ID) 

##reading the plink data into r for analysis
##using read.plink
ABO_FUT2_FUT3_plink_gtypes = read.plink("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/from local F/malariagen_analysis/r_analysis/Kenya_GWAS_bact_mal_b37_ABO_FUT2_FUT3", na.strings = c("00", "-9", "NA") )

bm_mgen_gtypes_matrix = cbind(ABO_FUT2_FUT3_plink_gtypes$fam$member, ABO_FUT2_FUT3_plink_gtypes$genotypes@.Data)
colnames(bm_mgen_gtypes_matrix)[1] = "IID"

##merging genotypes to phenotypes
bm_mgen_su_gwas_genopheno = merge(bm_su_mgen_cases, bm_mgen_gtypes_matrix, by = "IID", all.x = T)

bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>% select(c("IID", "Kenyan_ID", "Missing", "SEX", "cc_status", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PLATFORM", "mal_sub_status", "bact_sub_status", "cc_status_bact", "cc_status_mal", "mal_sub_status_sma", "mal_sub_status_cm", "mal_sub_status_both", "mal_sub_status_other", "mal_sub_status_bact", "bact_sub_status_aci", "bact_sub_status_bhs", "bact_sub_status_ecoli", "bact_sub_status_hib", "bact_sub_status_mal", "bact_sub_status_nts", "bact_sub_status_pneumo", "bact_sub_status_staph", "Oxford_ID", "Study", "missing_kenyanid", "Kenyan_ID_raw", "serial", "doa", "wbc", "rbc", "mps_100_wbc", "mps_500_rbc", "parasite_gn", "sample_code", "source_code_old", "agemths", "ageyr", "sex", "dob", "ethnicity", "ethnic", "loc_dss", "dss", "case", "died", "cm", "sma", "rd", "not_cm_sma_rd", "any_cm_sma_rd", "cm_sma_rd", "cm_not_sma_rd", "sma_not_cm_rd", "rd_not_cm_sma", "cm_and_sma_not_rd", "cm_and_rd_not_sma", "sma_and_rd_not_cm", "not_cm", "not_sma", "not_rd", "slide_pos", "bcstot", "prostration", "severe_malaria_anaemia", "cerebral_malaria", "respiratory_distress", "bpd", "bps", "pul_hrate", "crefill", "tempaxil", "oxysat", "fontanelle", "spleen", "muac", "height", "weight", "glucose", "dod", "diag1", "diag2", "hb", "hco3", "hct", "creat", "platelet", "mcv", "baseexc", "k", "na", "logpf_mcl", "phenotypic_bldgrp", "thall_type", "hbb_rs334", "g6pd_rs1050828", "g6pd_rs1050829", "inDSS", "cr1_rs17047660",  "cr1_rs17047661", "abo_rs8176719", "abo_rs8176746", "thal_gene", "thal_results", "exclsv_syndromes", "rs480133:C", "rs601338:A", "rs8176719:TC", "rs8176746:T", "caseorcontrol", "cm_sma_nonstrict"),contains("popn_"))

##incorporating sickle genotypes
bm_mgen_su_gwas_genopheno = read_csv("C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/JRWT Study data/Rotavirus cases/JG controls ethnicity and rs334 data/bacteraemia_for_JR.sample.csv") %>%
  select(c("ID1","rs334","ethnicity")) %>%
  left_join(bm_mgen_su_gwas_genopheno, ., by = c("IID" = "ID1")) %>%
  left_join(., bm_su_mgen_cases %>% select(c("IID","rsid.rs334.A.add", "compressed_ethnicity")), by = c("IID")) %>% 
  mutate(rsid.rs334.A.add = factor(round(as.numeric(rsid.rs334.A.add)), labels = c("Norm","Het","Hom")),
         rs334 = factor(rs334, labels = c("Norm","Het","Hom")),
         rs334_all = coalesce(rs334, rsid.rs334.A.add)
    
  )

##Formatting ethnicity
bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>%
  mutate(ethnicity.x = gsub("Kenya_", "", ethnicity.x),
         compressed_ethnicity_all = str_to_sentence(coalesce(ethnicity.x, compressed_ethnicity))
  )

##preparing variables for analysis by converting to proper datatypes and removing unnecessary NAs
bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>% 
  mutate_at(vars(c(`rs480133:C`, `rs601338:A`, `rs8176719:TC`, `rs8176746:T`)), ~na_if(., "00")) %>%
  mutate(phenotypic_bldgrp_6bg = factor(case_when(`rs8176719:TC` == "03" & `rs8176746:T` == "03" ~ "BB", 
                                                  `rs8176719:TC` == "03" & `rs8176746:T` == "02" ~ "AB", 
                                                  `rs8176719:TC` == "03" & `rs8176746:T` == "01" ~ "AA", 
                                                  `rs8176719:TC` == "02" & `rs8176746:T` == "03" ~ "BO" , 
                                                  `rs8176719:TC` == "02" & `rs8176746:T` == "02" ~ "BO" , 
                                                  `rs8176719:TC` == "02" & `rs8176746:T` == "01" ~ "AO", 
                                                  `rs8176719:TC` == "01" ~ "OO"), levels = c("OO", "AO", "AA", "BO", "BB", "AB")),
         phenotypic_bldgrp_4bg = factor(case_when(phenotypic_bldgrp_6bg == "BB" | phenotypic_bldgrp_6bg == "BB" ~ "B",
                                                  phenotypic_bldgrp_6bg == "AA" | phenotypic_bldgrp_6bg == "AO" ~ "A",
                                                  TRUE ~ as.character(phenotypic_bldgrp_6bg)), levels = c("OO", "A", "B", "AB")))

bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>% 
  mutate(
  non_OvsO = factor(
    case_when(phenotypic_bldgrp_4bg == "B" | phenotypic_bldgrp_4bg == "A" | phenotypic_bldgrp_4bg == "AB" ~ "non_O", 
              TRUE ~ as.character(phenotypic_bldgrp_4bg)), 
    levels = c("OO", "non_O")
  ),
  secretor_status = factor(
    case_when(`rs601338:A` == "03" ~ "nonsecretor", 
              `rs601338:A` == "01" | `rs601338:A` == "02" ~ "secretor"), 
    levels = c("secretor", "nonsecretor")))

##HBGA combinations
bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>%
  mutate(
    OvsnonO_secretor_combi = factor(
      case_when(`rs601338:A` != "03" & non_OvsO == "OO" ~ "O_secretor", 
                `rs601338:A` != "03" & non_OvsO == "non_O" ~ "nonO_secretor", 
                `rs601338:A` == "03" & non_OvsO == "OO" ~ "O_nonsecretor", 
                `rs601338:A` == "03" & non_OvsO == "non_O" ~ "nonO_nonsecretor"), levels = c("O_secretor", "O_nonsecretor", "nonO_secretor",  "nonO_nonsecretor")),
    ABO_secretor_combi = factor(
      case_when(`rs601338:A` != "03" & phenotypic_bldgrp_4bg == "OO" ~ "O_secretor",
                `rs601338:A` != "03" & phenotypic_bldgrp_4bg == "A" ~ "A_secretor",
                `rs601338:A` != "03" & phenotypic_bldgrp_4bg == "B" ~ "B_secretor",
                `rs601338:A` != "03" & phenotypic_bldgrp_4bg == "AB" ~ "AB_secretor",
                `rs601338:A` == "03" & phenotypic_bldgrp_4bg == "OO" ~ "O_nonsecretor", 
                `rs601338:A` == "03" & phenotypic_bldgrp_4bg == "A" ~ "A_nonsecretor",
                `rs601338:A` == "03" & phenotypic_bldgrp_4bg == "B" ~ "B_nonsecretor",
                `rs601338:A` == "03" & phenotypic_bldgrp_4bg == "AB" ~ "AB_nonsecretor"), 
      levels = c("O_secretor", "O_nonsecretor", "A_secretor",  "A_nonsecretor", "B_secretor",  "B_nonsecretor", "AB_secretor",  "AB_nonsecretor")
      )
  )

##hospitalization days calculation
bm_mgen_su_gwas_genopheno = bm_mgen_su_gwas_genopheno %>%
  mutate(hospitalization_days = as.numeric(as.duration(interval(doa, dod)), "days"),
         admission_year = factor(lubridate::year(doa)),
         age_years = as.numeric(as.duration(interval(dob, doa)), "years"),
         age_months = as.numeric(as.duration(interval(dob, doa)), "months")
  )

##Writing out final file
write.csv(bm_mgen_su_gwas_genopheno, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_mgen_su_gwas_genopheno_final_dset_for_analysis.csv", quote = F, row.names = F)

##analysis variables
bm_genopheno = bm_mgen_su_gwas_genopheno %>% select(c("IID", "Kenyan_ID", "SEX", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PLATFORM", "mal_sub_status", "bact_sub_status","cc_status", "cc_status_bact", "cc_status_mal", "mal_sub_status_sma", "mal_sub_status_cm", "mal_sub_status_both", "mal_sub_status_other", "mal_sub_status_bact", "bact_sub_status_aci", "bact_sub_status_bhs", "bact_sub_status_ecoli", "bact_sub_status_hib", "bact_sub_status_mal", "bact_sub_status_nts", "bact_sub_status_pneumo", "bact_sub_status_staph", "Study", "doa", "wbc", "rbc", "mps_100_wbc", "mps_500_rbc", "dob", "died", "cm", "sma", "rd", "not_cm_sma_rd", "any_cm_sma_rd", "cm_sma_rd", "cm_not_sma_rd", "sma_not_cm_rd", "rd_not_cm_sma", "cm_and_sma_not_rd", "cm_and_rd_not_sma", "sma_and_rd_not_cm", "not_cm", "not_sma", "not_rd", "bpd", "bps", "muac", "weight", "diag1", "diag2", "hb", "hct", "platelet", "mcv", "thall_type", "exclsv_syndromes", "rs480133:C", "rs601338:A", "rs8176719:TC", "rs8176746:T", "rs334_all", "compressed_ethnicity_all",  "phenotypic_bldgrp_4bg", "phenotypic_bldgrp_6bg", "non_OvsO", "secretor_status", "OvsnonO_secretor_combi", "ABO_secretor_combi", "hospitalization_days", "admission_year", "age_years", "age_months", "caseorcontrol", "cm_sma_nonstrict"), contains("popn_"))

bm_genopheno$IID = as.character(bm_genopheno$IID)

bm_genopheno = bm_genopheno %>% 
  mutate_at(vars(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6")), as.numeric) %>%
  mutate_at(vars(c("cc_status_bact", "cc_status_mal", "mal_sub_status_sma", "mal_sub_status_cm", "mal_sub_status_both", "mal_sub_status_other", "mal_sub_status_bact", "bact_sub_status_aci", "bact_sub_status_bhs", "bact_sub_status_ecoli", "bact_sub_status_hib", "bact_sub_status_mal", "bact_sub_status_nts", "bact_sub_status_pneumo", "bact_sub_status_staph", "Study", "cm", "sma", "rd", "not_cm_sma_rd", "any_cm_sma_rd", "cm_sma_rd", "cm_not_sma_rd", "sma_not_cm_rd", "rd_not_cm_sma", "cm_and_sma_not_rd", "cm_and_rd_not_sma", "sma_and_rd_not_cm", "not_cm", "not_sma", "not_rd", "diag1", "diag2", "exclsv_syndromes", "rs480133:C", "rs601338:A", "rs8176719:TC", "rs8176746:T", "rs334_all", "phenotypic_bldgrp_4bg", "non_OvsO", "secretor_status", "OvsnonO_secretor_combi", "ABO_secretor_combi","caseorcontrol", "cm_sma_nonstrict")), as.factor) %>%
  mutate(`rs480133:C` = factor(`rs480133:C`, labels = c("TT", "TC", "CC")),
         died = factor(died, labels = c("No", "Yes")),
         SEX = factor(SEX, labels = c("Male", "Female"))) 

##Removing participants
##Removing individuals in the exclusion list for JG bact dataset
bm_genopheno_jg_bact = bm_genopheno[which(!(bm_genopheno[, "IID"] %in% mgen_exclusion_list$V1)),]

##Retaining 3068 mgen malaria participants from 7831 total individuals from James Gilchrist (JG)
bm_genopheno_3068 = bm_genopheno %>% 
  filter(!is.na(popn_PC_1) & caseorcontrol != "PARENT")
  
##Removing sicklers
#bm_genopheno_3068 = bm_genopheno_3068[bm_genopheno_3068[,"rs334_all"] != "TT",]

##rs601338 HWE
bm_genopheno_3068 %>% filter(!is.na(`rs601338:A`) & cc_status == "CONTROL") %>% count(`rs601338:A`) %>% with(HWChisq(n))

##rs480133 HWE
##distribution of secretors in cases and controls
bm_genopheno_3068 %>% filter(!is.na(`rs480133:C`) & cc_status == "CONTROL") %>% count(`rs480133:C`) %>% with(HWChisq(n))

##rs8176719 HWE
bm_genopheno_3068 %>% filter(!is.na(`rs8176719:TC`) & cc_status == "CONTROL") %>% count(`rs8176719:TC`) %>% with(HWChisq(n))

##rs8176746 HWE
bm_genopheno_3068 %>% filter(!is.na(`rs8176746:T`) & cc_status == "CONTROL") %>% count(`rs8176746:T`) %>% with(HWChisq(n))


##writing out final d_set
write.csv(bm_genopheno_3068, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/malariagen_data_raw/raw_processed/bm_genopheno_3068_final.csv" , quote = F, row.names = F)

##Declaring variables
a_b_ab_oo = c("A", "B", "AB", "OO")

bact_syndromes = c("cc_status_bact", "mal_sub_status_bact", "bact_sub_status_aci", "bact_sub_status_bhs", "bact_sub_status_ecoli", "bact_sub_status_hib", "bact_sub_status_mal", "bact_sub_status_nts", "bact_sub_status_pneumo", "bact_sub_status_staph")

#"cc_status_mal", "mal_sub_status_sma", "mal_sub_status_all_sma", "mal_sub_status_cm", "mal_sub_status_all_cm", "mal_sub_status_both", "mal_sub_status_other",
snpid = c("secretor_status", "`rs480133:C`")

##Labelling variables appropriately
# rcc_unit = expression(paste("Red cell count (10"^"12","/L)"))
bm_genopheno_3068 = bm_genopheno_3068 %>%
  set_variable_labels(secretor_status = "Secretor status", phenotypic_bldgrp_4bg = "ABO 4", phenotypic_bldgrp_6bg = "ABO 6", hospitalization_days = "Hospitalization duration (Days)", weight = "Weight (Kg)", muac = "Mid upper arm circumference (cm)", compressed_ethnicity_all = "Ethnicity", exclsv_syndromes = "Exclusive malaria syndromes", age_years = "Age (years)", SEX = "Sex", rs334_all = "Sickle", thall_type = "Alpha thalassaemia", `rs480133:C` = "rs480133 genotype", hb = "Haemoglobin (g/dL)", hct = "Haematocrit (%)", mcv = "Mean cell volume (fL)", platelet = "Platelets (10^9/L)", died = "Died during admission", wbc = "White cell count (10^9/L)", rbc = "Red cell count (10^12/L)", mps_100_wbc = "Parasites in 100 WBCs", mps_500_rbc = "Parasites in 500 WBCs", cm = "Cerebral malaria", sma = "Severe malaria anaemia", rd = "Respiratory distress", non_OvsO = "O VS Non-O")

##28072019 Number of cases and controls in individuals included in logistic regression model i.e with complete case, sex, ethnicity, hbb_rs334 and secretor status info
##Table 1 summary
tbl_cases = bm_genopheno_3068 %>% 
  filter(cc_status == "MAL") %>% 
  dplyr::select(age_years, SEX, compressed_ethnicity_all, phenotypic_bldgrp_4bg, rs334_all, thall_type, `rs480133:C`, hospitalization_days, exclsv_syndromes, hb, hospitalization_days, weight, muac, secretor_status, mps_500_rbc, died, cm, sma, rd) %>%
  mutate_at(vars(cm, sma, rd), droplevels) %>%
  tbl_summary(., by = secretor_status, missing = "ifany", type = list(all_dichotomous() ~ "categorical")) %>% 
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>% 
  add_stat_label()%>%
  bold_p() #%>% 
# add_p_footnotes()

tbl_ctrls = bm_genopheno_3068 %>% 
  dplyr::filter(cc_status == "CONTROL") %>% 
  dplyr::select(age_years, SEX, compressed_ethnicity_all, phenotypic_bldgrp_4bg, rs334_all, thall_type, `rs480133:C`, secretor_status) %>%
  tbl_summary(., by = secretor_status, missing = "ifany") %>% 
  add_p(pvalue_fun = function(x) style_pvalue(x, digits = 2)) %>% 
  add_stat_label()%>%
  bold_p()

tbl_demog = tbl_merge(list(tbl_cases, tbl_ctrls),tab_spanner = c("Cases", "Controls")) %>%  
  gtsummary::as_flextable() %>% 
  padding(padding=0) %>% 
  set_table_properties(layout = "autofit") %>% 
  flextable::fontsize(size = 8, part = "all") %>% 
  flextable::border(i = ~ label %in% c("Age (years), median (IQR)", "Sex, n (%)", "Hospitalization duration (Days), median (IQR)", "Weight (Kg), median (IQR)", "Mid upper arm circumference (cm), median (IQR)", "Ethnicity, n (%)", "Haemoglobin (g/dL), median (IQR)", "O VS Non-O, n (%)", "ABO 4, n (%)", "Mixed & non-typable P genotypes, n (%)", "Sickle, n (%)", "Alpha thalassaemia, n (%)", "Died during admission, n (%)", "Parasites in 500 WBCs, median (IQR)", "cm, n (%)", "sma, n (%)", "rd, n (%)", "rs480133 genotype, n (%)", "Exclusive malaria syndromes, n (%)"), border.top = fp_border(color="gray70", width = 1), part = "body") %>%
  bold(i = ~ label %in% c("Age (years), median (IQR)", "Sex, n (%)", "Hospitalization duration (Days), median (IQR)", "Weight (Kg), median (IQR)", "Mid upper arm circumference (cm), median (IQR)", "Ethnicity, n (%)", "Haemoglobin (g/dL), median (IQR)", "O VS Non-O, n (%)", "ABO 4, n (%)", "Mixed & non-typable P genotypes, n (%)", "Sickle, n (%)", "Alpha thalassaemia, n (%)", "Died during admission, n (%)", "Parasites in 500 WBCs, median (IQR)", "cm, n (%)", "sma, n (%)", "rd, n (%)", "rs480133 genotype, n (%)", "Exclusive malaria syndromes, n (%)"), j=1, bold = TRUE, part = "body") %>%
  bold(bold = TRUE, part = "header")


##ANALYSIS
##Univariate analysis
sink("possible_covariates_univariate.txt")
posbl_predctr = c( "rs334_all", "thall_type")
for (i in 1:length(bact_syndromes)){
  for (k in 1:length(posbl_predctr)){
    print(paste(">>>>>>>>", bact_syndromes[i], "<<<<<<<<", "/n"))
    print(paste("WITH MALARIA----------------UNIVARIATE"))
    abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  posbl_predctr[k], "+" ,  paste(colnames(bm_genopheno_3068)[c(4:9)], collapse = "+") )),data = bm_genopheno_3068 , family = "binomial")
    print(summary(abo_fut_model))
    print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
    print(nobs(abo_fut_model))
    #print(paste("MINUS MALARIA----------------UNIVARIATE"))
    #abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  posbl_predctr[k], "+" ,  paste(colnames(bm_genopheno_3068)[c(4:9)], collapse = "+") )),data = bm_genopheno_3068[bm_genopheno_3068[,"malaria"] !="YES"| is.na(bm_genopheno_3068[,"malaria"]),] , family = "binomial")
    #print(summary(abo_fut_model))
    #print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
    #print(nobs(abo_fut_model))
    
  }
}
sink()

##counts with sickle
counts_with_sickle = bm_genopheno_3068 %>% filter(!is.na(cc_status) & !is.na(ethnicity))  %>% count(cc_status, bact_sub_status, sickle)
write.csv(counts_with_sickle, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_sickle.csv", quote = F, row.names = F)

#within everyone
snp_log_rgrsn_everyone = function (predctr) {
  for (i in 1:length(bact_syndromes)){
    
    for (k in 1:length(predctr)){
      print(paste(">>>>>>>>", bact_syndromes[i], "<<<<<<<<"))
      ##including accounting for ethnicity instead of PCs
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k], "+" ,  paste(colnames(bm_genopheno_3068)[c(4:9)], collapse = "+") )),data = bm_genopheno_3068, family = "binomial")
      print(summary(abo_fut_model))
      print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
      print(nobs(abo_fut_model))
      ##10062019 unadjusted
      print(paste("----------------unadjusted----------------"))
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k] )),data = bm_genopheno_3068, family = "binomial")
      print(summary(abo_fut_model))
      print(exp(cbind(OR = coef(abo_fut_model), confint(abo_fut_model))))
      print(nobs(abo_fut_model))
    }
    
  }
}

snp_forest_everyone = function (predctr) {
  for (i in 1:length(bact_syndromes)){
    
    for (k in 1:length(predctr)){
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k], "+" ,  paste(colnames(bm_genopheno_3068)[c(4:9)], collapse = "+") )),data = bm_genopheno_3068, family = "binomial")
      print(plot_model(abo_fut_model, show.values = TRUE, value.offset = .2, value.size = 6, show.p = T, colors = "bw", vline.color = "red",transform = "exp", title = paste(bact_syndromes[i]), rm.terms = c("sickleAT","sickleTT",  "PC_1", "popn_PC_2", "popn_PC_3", "popn_PC_4", "popn_PC_5", "popn_PC_6", "thall_typeHet", "thall_typeHomo"), order.terms =  c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17) ) + theme_sjplot(base_size = 18))
      abo_fut_model = glm(as.formula(paste(bact_syndromes[i], "~",  predctr[k] )),data = bm_genopheno_3068, family = "binomial")
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
counts_with_O_nonO = bm_genopheno_3068 %>% filter(!is.na(cc_status) & !is.na(non_OvsO_jcr))  %>% count(cc_status, P4_exclusive, P6_exclusive, P8_exclusive,non_OvsO_jcr, !is.na(ethnicity) )
write.csv(counts_with_O_nonO, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_O_nonO.csv", quote = F, row.names = F)

##rs601338

sink("bact_logistic_regression_rs601338_everyone.txt")
snp_log_rgrsn_everyone(predctr = "secretor_status")
sink()

pdf("bact_logistic_regression_rs601338_everyone.pdf")
snp_forest_everyone(predctr = "secretor_status")
dev.off()

counts_with_secretor = bm_genopheno_3068 %>% filter(!is.na(cc_status)  & !is.na(secretor_status))  %>% count(cc_status, P4_exclusive, P6_exclusive, P8_exclusive,secretor_status, !is.na(ethnicity))
write.csv(counts_with_secretor, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_secretor.csv", quote = F, row.names = F)

sink("bact_logistic_ABO_secretor_combi.txt")
snp_log_rgrsn_everyone(predctr = "ABO_secretor_combi")
sink()

pdf("bact_logistic_ABO_secretor_combi.pdf")
snp_forest_everyone(predctr = "ABO_secretor_combi")
dev.off()

counts_with_ABO_secretor_combi = bm_genopheno_3068 %>% filter(!is.na(cc_status) & !is.na(ABO_secretor_combi) & !is.na(P8_exclusive))  %>% count(P8_exclusive, ABO_secretor_combi, !is.na(ethnicity))
write.csv(counts_with_ABO_secretor_combi, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_ABO_secretor_combi.csv", quote = F, row.names = F)

sink("bact_logistic_OvsnonO_secretor_combi.txt")
snp_log_rgrsn_everyone(predctr = "OvsnonO_secretor_combi")
sink()

pdf("bact_logistic_OvsnonO_secretor_combi.pdf")
snp_forest_everyone(predctr = "OvsnonO_secretor_combi")
dev.off()

counts_with_OvsnonO_secretor_combi = bm_genopheno_3068 %>% filter(!is.na(cc_status) & !is.na(P8_exclusive) & !is.na(OvsnonO_secretor_combi))  %>% count( P8_exclusive, OvsnonO_secretor_combi, !is.na(ethnicity))
write.csv(counts_with_OvsnonO_secretor_combi, "C:/Users/Jesse Rop/OneDrive - Kemri Wellcome Trust/WT 18 month project/malariagen_analysis/bact_mal/bm_counts/counts_with_OvsnonO_secretor_combi.csv", quote = F, row.names = F)

sink("bact_logistic_rs480133.txt")
snp_log_rgrsn_everyone(predctr = "`rs480133:C`")
sink()


