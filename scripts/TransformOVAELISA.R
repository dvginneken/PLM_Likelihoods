#Binder Properties

#Clean the data
df_b <- read.csv("OneDrive - UMC Utrecht/Documenten/Students/Karlis/Binders_May31_fixed.csv", sep = ";")
df_b <- df_b[,c("HC_TRANSLATED_AA", "LC_TRANSLATED_AA","Bind..ELISA.signal.0.2.", "octet.affinity..nM.",
             "Clone_id", "Mouse_clone_HC", "Mouse_clone_LC", "HC_CDR3","LC_CDR3","VDJ_Vgene","VDJ_Jgene",
             "VDJ_Cgene","VJ_Cgene","VJ_Vgene","VJ_Jgene","VDJ_CDR3_aa","VJ_CDR3_aa")]
df_nb <- read.csv("OneDrive - UMC Utrecht/Documenten/Students/Karlis/Nonbinders_May31_fixed.csv", sep = ";")
df_nb <- df_nb[,c("HC_TRANSLATED_AA", "LC_TRANSLATED_AA","Bind..ELISA.signal.0.2.", "octet.affinity..nM.",
                "Clone_id", "Mouse_clone_HC", "Mouse_clone_LC", "HC_CDR3","LC_CDR3","VDJ_Vgene","VDJ_Jgene",
                "VDJ_Cgene","VJ_Cgene","VJ_Vgene","VJ_Jgene","VDJ_CDR3_aa","VJ_CDR3_aa")]
df <- rbind(df_b, df_nb)

#Clean up
df$Bind..ELISA.signal.0.2. <- case_match(df$Bind..ELISA.signal.0.2.,
                                         "yes (own ELISA)" ~ "yes",
                                         "yes" ~ "yes",
                                         "no" ~ "no")
df$full_sequence <- df$HC_TRANSLATED_AA
df$full_sequence
write.csv(df, file = "OneDrive - UMC Utrecht/Documenten/Project_EvoLikelihood/PLM-likelihoods/data/OVA_V7/Binder_properties_cleaned.csv",
          row.names = F)
