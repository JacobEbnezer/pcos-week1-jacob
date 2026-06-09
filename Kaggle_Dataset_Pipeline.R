rm(list = ls())

# Setting directory
setwd("C:/Users/jacob/Desktop/pcos-week1-jacob") 

# 2. Loading the Kaggle dataset
clinical_data <- read.csv("PCOS_infertility.csv")

# 3. Printing the size of the dataset and the column names
print("Dataset Dimensions (Rows, Columns):")
print(dim(clinical_data))

print("Column Names:")
print(colnames(clinical_data))

# The amh text problem and convert to numeric
clinical_data$AMH.ng.mL. <- as.numeric(clinical_data$AMH.ng.mL.)

# Remove the row with missing data
clean_clinical <- na.omit(clinical_data)

# Average amh levels for pcos vs control
aggregate(AMH.ng.mL. ~ PCOS..Y.N., data = clean_clinical, FUN = median)

# Average beta-hcg levels for pcos vs control
aggregate(I...beta.HCG.mIU.mL. ~ PCOS..Y.N., data = clean_clinical, FUN = median)

# finding if the amh difference is statistically significant
wilcox.test(AMH.ng.mL. ~ PCOS..Y.N., data = clean_clinical)

# plot
boxplot(AMH.ng.mL. ~ PCOS..Y.N., data = clean_clinical,
        main = "AMH Levels: Control vs PCOS",
        ylab = "AMH (ng/mL)", 
        xlab = "Group (0 = Control, 1 = PCOS)",
        col = c("yellow", "tomato"), 
        outline = FALSE) # hiding extreme outliers so the boxes are readable
