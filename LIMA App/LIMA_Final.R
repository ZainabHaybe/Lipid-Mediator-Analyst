
# ------> LOAD PACKAGES:
if (!require('shiny')) install.packages('shiny'); library('shiny')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('plotly')) install.packages('plotly'); library('plotly')
if (!require('corrplot')) install.packages('plotly'); library('corrplot')
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require('mixOmics')) BiocManager::install('mixOmics'); library('mixOmics')
if (!require('ropls')) BiocManager::install('ropls'); library('ropls')
if (!require('DT')) install.packages('DT'); library('DT')
if (!require('MVN')) install.packages('MVN'); library('MVN')
if (!require('randomForest')) install.packages('randomForest'); library('randomForest')
if (!require('gridExtra')) install.packages('gridExtra'); library('gridExtra')
if (!require('snowfall')) install.packages('snowfall')
if (!require('neldermead')) install.packages('neldermead')
if (!require('optimbase')) install.packages('optimbase')
if (!require('classyfire')) install.packages("classyfire_0.1-2.tar.gz",
                                             repos = NULL,
                                             type = "source"); library('classyfire')
if (!require('glmnet')) install.packages('glmnet'); library('glmnet')
if (!require('caret')) install.packages('caret'); library('caret')
if (!require('arm')) install.packages('arm'); library('arm')
if (!require('xgboost')) install.packages('xgboost'); library('xgboost')
if (!require('DiagrammeR')) install.packages('DiagrammeR'); library('DiagrammeR')
if (!require('Matrix')) install.packages('Matrix'); library('Matrix')
if (!require('Ckmeans.1d.dp')) install.packages('Ckmeans.1d.dp'); library('Ckmeans.1d.dp')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('tibble')) install.packages('tibble'); library('tibble')
if (!require('pROC')) install.packages('pROC'); library('pROC')
#if (!require('bslib')) install.packages('bslib');library('bslib')
if (!require('shinythemes')) install.packages('shinythemes');library('shinythemes')


#set working directory to current location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



###Open rename and fatty acid file needed throughout multiple pages 
Rename_Met_Table<-read.csv(file = "coefficients/Rename_Met2.csv", header = TRUE, sep = ",") #open file
colnames(Rename_Met_Table)<- c("Final_Name", "Other_Name") #ensure column names are standardized 

#Create a dictionary with values as the old names and the rename as the names
Rename_Met<-setNames(Rename_Met_Table$Final_Name, Rename_Met_Table$Other_Name)

# #dictionary for SPMS
SPM_FA_Table<- read.csv(file = "coefficients/SPM_FA.csv", header = TRUE, sep = ",")
colnames(SPM_FA_Table)<-c("Metabolite", "Fatty_Acid") #ensure column names are standardized

#create a dictionary with the metabolites as values and the names as fatty acids 
SPMs_FA<-setNames(SPM_FA_Table$Fatty_Acid, SPM_FA_Table$Metabolite)



####UI = the user interface script with the inputs and outputs for each page 
ui <- navbarPage(title = "LIMA",
                 #theme for the application 
                 theme = shinytheme("flatly"),
                 
                 #menu page for the data preparation 
                 navbarMenu("Data Preparation",
                            
                            #first page with the data progessing
                            tabPanel(
                              
                              #name of page 
                              title ="Data Processing",
                              
                              #ID is how it is referenced in other code
                              id = "data",
                              
                              #Page title 
                              titlePanel("Data Processing"),
                              
                              #sidebar layout with input and output definitions 
                              sidebarLayout(
                                
                                #Sidebar panel for inputs
                                sidebarPanel(
                                  #size of sidebar
                                  width = 3,
                                  
                                  #input: Select a file using fileInput()
                                  fileInput("Data_File", "Choose Raw MS file",
                                            multiple = FALSE,
                                            #accept is list of file types that are acceptable
                                            accept = c("text/csv", 
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")),
                                  
                                  #checkbox input to see if table has headers 
                                  checkboxInput("Header_Data", "Header", TRUE),
                                  
                                  #Input: select separator 
                                  radioButtons("sep_Data", "Separator",
                                               choices = c(Comma = ",",
                                                           Tab = "\t"),
                                               selected = ","),
                                  #line separator
                                  tags$hr(),
                                  
                                  #option for user to use their own coefficients 
                                  checkboxInput("Coef_Upload", "Do you want to use your own coefficeints?", FALSE),
                                  
                                  #conditional panel that appears if user wants to use set coefficients (Checkbox input = FALSE)
                                  conditionalPanel(
                                    condition = "input.Coef_Upload == false",
                                    
                                    #Button for user to select which Mass Spec machine used which coordinates with coefficient values 
                                    radioButtons("MS_Machine", "Please select Mass Spec Machine",
                                                 choices = c("Mass Spec 2" = "MS2",
                                                             "Mass Spec 3" = "MS3",
                                                             "Mass Spec 4" = "MS4"),
                                                 selected = "MS2")),
                                  
                                  #conditional panel that appears if user wants to use their own coefficients (Checkbox input = TRUE)
                                  conditionalPanel(
                                    condition = "input.Coef_Upload == true",
                                    
                                    fileInput("CrossRef_File", "Choose Cross Reference Coefficent File",
                                              multiple = FALSE,
                                              #accept is list of file types that are acceptable
                                              accept = c("text/csv", 
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    #checkbox input to see if table has headers 
                                    checkboxInput("Header_CrossRef_File", "Header", TRUE),
                                    
                                    #Input: select separator 
                                    radioButtons("sep_CrossRef_File", "Separator",
                                                 choices = c(Comma = ",",
                                                             Tab = "\t"),
                                                 selected = ","),
                                    #download example cross reference coeficient file 
                                    downloadButton("DownloadCrossRef_File", "Example Cross Reference File"),
                                    
                                    #line separator
                                    tags$hr(),
                                    
                                    fileInput("ConvCoef_File", "Choose  Standard Curv Coefficient file",
                                              multiple = FALSE,
                                              #accept is list of file types that are acceptable
                                              accept = c("text/csv", 
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    #checkbox input to see if table has headers 
                                    checkboxInput("Header_ConvCoef_File", "Header", TRUE),
                                    
                                    #Input: select separator 
                                    radioButtons("sep_ConvCoef_File", "Separator",
                                                 choices = c(Comma = ",",
                                                             Tab = "\t"),
                                                 selected = ","),
                                    #download sample list file
                                    downloadButton("DownloadConvCoef_File", "Example Standard Curve Coefficient File")
                                  ),
                                  
                                  
                                  #line separator
                                  tags$hr(),
                                  
                                  #if sample has cross refence 
                                  checkboxInput("CrossRef_Incl", "Does samples include Cross Reference?", TRUE),
                                  
                                  conditionalPanel(
                                    
                                    condition = "input.CrossRef_Incl == true",
                                    
                                    #Input - Cross Reference Name 
                                    selectInput('Cross_Name', 'Select Name of Cross Reference Sample', "")
                                  ),
                                  
                                  #select excess columns that are not sampels 
                                  selectInput('not_Samples', 'Select all that are not samples', "", multiple=TRUE, selectize=TRUE),
                                  
                                  #Input - Max cokumn name name 
                                  selectInput('MaxCol_Name', 'Select Name of 100% Sample', ""),
                                  
                                  #Input: Filter threshold for Area 
                                  numericInput("Area_Filter", "Select Area to filter to zero ", 0, min = 0, max = NA),
                                  
                                  
                                  #line separator
                                  tags$hr(),
                                  
                                  checkboxInput("SameVol_Data", "Are your samples all the same volumne/weight?", TRUE),
                                  conditionalPanel(
                                    condition = "input.SameVol_Data == true",
                                    #input volume of samples 
                                    numericInput("SameData_Amount", "Select the volume/weight of your sample", 1, min = 0, max = NA),
                                    
                                    numericInput("SameData_Standard", "Select the volume/weight you want to standardize the measurement to.", 1, min = 0, max = NA),
                                    ),
                                  
                                  checkboxInput("DiffVol_Data", "Are your samples different volume/Weight?", FALSE),
                                  conditionalPanel(
                                    condition = "input.DiffVol_Data == true",
                                    #download sample list file
                                    downloadButton("DownloadSampleList", "Download Sample List"),
                                    fileInput("DiffVol_File", "Choose Samples Volume/Mass File",
                                              multiple = FALSE,
                                              #accept is list of file types that are acceptable
                                              accept = c("text/csv", 
                                                         "text/comma-separated-values,text/plain",
                                                         ".csv")),
                                    #checkbox input to see if table has headers 
                                    checkboxInput("Header_DiffVol", "Header", TRUE),
                                    
                                    #Input: select separator 
                                    radioButtons("sep_DiffVol", "Separator",
                                                 choices = c(Comma = ",",
                                                             Tab = "\t"),
                                                 selected = ","),
                                    numericInput("DiffData_Standard", "Select the volume/weight you want to standardize the measurement to.", 1, min = 0, max = NA)
                                  ),
                                  
                                  #line separator
                                  tags$hr(),
                                  
                                  
                                  #calculate button 
                                  actionButton("Go_Data", "Process")
                                  
                                  
                                  
                                ),
                                mainPanel(
                                  
                                  #ensures plots fit into tabset 
                                  width = 9,
                                  height = 12,
                                  tabsetPanel(type="tabs",
                                              tabPanel("Final Table", 
                                                       #Error Message text output
                                                       htmlOutput("ConcTable_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#ConcTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       dataTableOutput("ConcTable"),
                                                       downloadButton("DownloadConcTable", "Download Table")),
                                              tabPanel("Raw Data", 
                                                       #Error Message text output
                                                       htmlOutput("Raw_Table_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#Raw_Table_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       #datatable output if there is no error
                                                       dataTableOutput("Raw_Table"),
                                                       downloadButton("DownloadRaw_Table", "Download Table")), 
                                              tabPanel("Cross Reference Coefficent", 
                                                       #Error Message text output
                                                       htmlOutput("CrossRefTable_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#CrossRefTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       dataTableOutput("CrossRefTable"),
                                                       downloadButton("DownloadCrossRefTable", "Download Table")),
                                              tabPanel("Amount in pg", 
                                                       #Error Message text output
                                                       htmlOutput("AmountTable_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#AmountTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       dataTableOutput("AmountTable"),
                                                       downloadButton("DownloadAmountTable", "Download Table")),
                                              tabPanel("% Recovery",
                                                       #Error Message text output
                                                       htmlOutput("PerRecovery_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#PerRecovery_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       dataTableOutput("PerRecovery"),
                                                       downloadButton("DownloadPerRecoveryTable", "Download Table")),
                                              tabPanel("Normalized Amount", 
                                                       #Error Message text output
                                                       htmlOutput("NormTable_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#NormTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       dataTableOutput("NormTable"),
                                                       downloadButton("DownloadNormTable", "Download Table")),
                                              tabPanel("Coefficients Info", 
                                                       textOutput("Cross Reference Coefficients"),
                                                       dataTableOutput("CrossRefCoef_Table"),
                                                       tags$hr(),
                                                       textOutput("Conversion Coefficients"),
                                                       dataTableOutput("ConvCoef_Table")))
                                  ))),
                            
                            tabPanel(
                              #name of page 
                              title ="Data Formatting",
                              
                              #ID is how it is referenced in other code
                              id = "data2",
                              
                              #Page title 
                              titlePanel("Data Formating"),
                              
                              #sidebar layout with input and output definitions 
                              sidebarLayout(
                                
                                #Sidebar panel for inputs
                                sidebarPanel(
                                  #size of sidebar
                                  width = 3,
                                  
                                  #input: Select a file using fileInput()
                                  fileInput("Data_File2", "Choose Concentration Table file",
                                            multiple = FALSE,
                                            #accept is list of file types that are acceptable
                                            accept = c("text/csv", 
                                                       "text/comma-separated-values,text/plain",
                                                       ".csv")),
                                  
                                  #checkbox input to see if table has headers 
                                  checkboxInput("Header_Data2", "Header", TRUE),
                                  
                                  #Input: select separator 
                                  radioButtons("sep_Data2", "Separator",
                                               choices = c(Comma = ",",
                                                           Tab = "\t"),
                                               selected = ","),
                                  
                                  #line separator
                                  tags$hr(),
                                  
                                  #if sample has cross refence 
                                  checkboxInput("Filtering_Data", "Would you like to filter which mediators to include?", FALSE),
                                  
                                  conditionalPanel(
                                    
                                    condition = "input.Filtering_Data == true",
                                    
                                    #Input - Cross Reference Name 
                                    selectInput('Met_DataList2', 'Select all Mediators of Interest', "Please upload file", multiple=TRUE, selectize=TRUE)),
                                  
                                  #calculate button 
                                  actionButton("Go_Data_Format", "Format Data")),
                                mainPanel(
                                  
                                  #ensures plots fit into tabset 
                                  width = 9,
                                  height = 12,
                                  tabsetPanel(type="tabs",
                                              tabPanel("Statistics Table", 
                                                       #Error Message text output
                                                       textOutput("StatsInputTable_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#StatsInputTable_Error{color: red; font-size: 25px; font-style: italic;}")),
                                                       dataTableOutput("Stats_Table"),
                                                       downloadButton("DownloadStats_Table", "Download Table")),
                                              tabPanel("Machine Learning Table", 
                                                       textOutput("MLInputTable_Error"),
                                                       #styling text output to be red and bigger
                                                       tags$head(tags$style("#MLInputTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                       dataTableOutput("ML_Input_Table"),
                                                       downloadButton("DownloadML_Input_Table", "Download Table"))))
                                ))),
                 
                 navbarMenu("Multivariate Statistics",
                            tabPanel(title ="PCA/PLS-DA",
                                     id = "PCA/PLS-DA",
                                     titlePanel("PCA/PLS-DA Analysis"),
                                     sidebarLayout(sidebarPanel(
                                       #size of sidebar
                                       width = 3,
                                       
                                       #File upload code 
                                       fileInput("PCAPLSDA_file", "Choose upload file",
                                                 multiple = FALSE,
                                                 accept = c("text/csv", 
                                                            "text/comma-separated-values,text/plain",
                                                            ".csv")),
                                       
                                       #Check box to see if file has headers 
                                       checkboxInput("header_PCA", "Header", TRUE),
                                       
                                       #input--Separator for the file 
                                       radioButtons("sep_PCA", "Separator", 
                                                    choices = c(Comma = ",",
                                                                Tab = "\t")),
                                       
                                       #input---Location of samples 
                                       radioButtons("samp_PCA", "Sample location",
                                                    choices = c(Row = "row_PCA", 
                                                                Column = "col_PCA"),
                                                    selected = "row_PCA"),
                                       
                                       #input -- where the sample names column is
                                       #numericInput("SampCol_PCA", "Sample name Column/Row Number", 1, min = 1, max = NA),
                                       selectInput("SampleCol_PCA", "Select Sample Column Name", "Please upload file"),
                                       helpText("If samples were the rownames, please select x."),
                                       #input -- where the sample group info column is
                                       selectInput("Group_PCA", "Select Group Name", "Please upload file"),
                                       
                                       
                                       
                                       #line separator 
                                       tags$hr(),
                                       
                                       #Input --- Selecting the test
                                       radioButtons("test_PCAPLSDA", "Select test",
                                                    choices = c(PCA = "PCA",
                                                                PLSDA = "PLS-DA"),
                                                    selected = "PCA"),
                                       
                                       #line separator
                                       tags$hr(),
                                       
                                       #calculate button 
                                       actionButton("Go_PCAPLSDA", "Calculate")
                                     ),
                                     mainPanel(
                                       #ensures plots fit into tabset 
                                       width = 9,
                                       height = 12,
                                       
                                       tabsetPanel(type="tabs", id ="PCA/PLSDA_tabs",
                                                   tabPanel("2D Scores Plot",
                                                            #Error Message text output
                                                            htmlOutput("Xscoresplot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Xscoresplot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            plotlyOutput("Xscoresplot", height = "600px"),
                                                            #Check box to show samples 
                                                            checkboxInput("ViewSamples_PCA", "View Samples", TRUE),
                                                            #checkbox to select condidence interval
                                                            numericInput("PCA_Confidence", "Select The condience Interval ", 0.95, min = 0, max = 1),
                                                            #Update information 
                                                            actionButton("Go_PCAPLSDA2", "Update Plot"),
                                                            downloadButton("DownloadScore", "Download X-scores Data"),
                                                            tags$hr(),
                                                            uiOutput("R2Q2_Title"),
                                                            dataTableOutput("R2Q2_Table"),
                                                            plotlyOutput("R2Q2_Plot")),
                                                   tabPanel("Loading Plot", 
                                                            #Error Message text output
                                                            htmlOutput("Loadingplot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Loadingplot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            plotOutput("Loadingplot", height = "800px"),
                                                            downloadButton("downloadLoading_Plot", "Download Loading Plot"),
                                                            downloadButton("DownloadLoading", "Download Loading Score Data")),
                                                   tabPanel("Variance Plot", 
                                                            #Error Message text output
                                                            htmlOutput("PerVarplot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#PerVarplot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            plotlyOutput("PerVarplot", height = "600px"),
                                                            downloadButton("DownloadVar", "Download Variance Data")),
                                                   tabPanel(title = "VIP Scores", 
                                                            #Error Message text output
                                                            htmlOutput("VIPplot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#VIPplot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            id = "VIPScore",
                                                            textOutput("VIPtext"),
                                                            plotOutput("VIPplot", height = "800px", width = "800px"), 
                                                            tags$hr(),
                                                            #dataTableOutput("VIPtable"),
                                                            downloadButton("downloadPLSDA_VIP_Plot", "Download VIP plot"),
                                                            downloadButton("DownloadVIP", "Download VIP Table"),
                                                            tags$hr(),
                                                            uiOutput("VIP_Pathway_Title"),
                                                            plotlyOutput("VIP_PathwayPlot")
                                                            ))
                                     ))),
                            tabPanel(title ="Differential Analysis",
                                     id = "Diff",
                                     titlePanel("Differential Analysis"),
                                     sidebarLayout(sidebarPanel(
                                       #size of sidebar
                                       width = 3,
                                       
                                       #File upload code 
                                       fileInput("lm_profiles", "Choose file to upload",
                                                 multiple = FALSE,
                                                 accept = c("text/csv", 
                                                            "text/comma-separated-values,text/plain",
                                                            ".csv")),
                                       
                                       #Check box to see if file has headers 
                                       checkboxInput("header_Diff", "Header", TRUE),
                                       
                                       #input--Separator for the file 
                                       radioButtons("sep_Diff", "Separator", 
                                                    choices = c(Comma = ",",
                                                                Tab = "\t")),
                                       
                                       #input---Location of samples 
                                       radioButtons("samp_Diff", "Sample location",
                                                    choices = c(Row = "row_Diff", 
                                                                Column = "col_Diff"),
                                                    selected = "row_Diff"),
                                       
                                      
                                       #line separator 
                                       tags$hr(),
                                       
                                       # Input: Select separator ----
                                       radioButtons("mvn", "Multivariate Normality Test",
                                                    choices = c(Mardia = "mardia",
                                                                `Henze-Zirkler` = "hz",
                                                                Royston = "royston"),
                                                    selected = "royston"),
                                       
                                       #line separator
                                       tags$hr(),
                                       
                                       # Group selection ----
                                       selectInput("group_A", "Group A",""),
                                       selectInput("group_B", "Group B",""),
                                       
                                       
                                       #calculate button 
                                       actionButton("Go_Diff", "Calculate")
                                       
                                       # #refresh button
                                       # actionButton("reset_diff", "Reset inputs")
                                     ),
                                     
                                     mainPanel(
                                       #ensures plots fit into tabset 
                                       width = 9,
                                       height = 12,
                                       
                                       tabsetPanel(type="tabs", id ="Diff_tabs",
                                                   tabPanel("Results Table",
                                                            #Error Message text output
                                                            htmlOutput("DiffTable_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#DiffTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            dataTableOutput("DiffTable"),
                                                            downloadButton("DownloadDiffTable", "Download Data")),
                                                   tabPanel("Volcano Plots", 
                                                            #Error Message text output
                                                            htmlOutput("DiffPlot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#DiffPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            numericInput("Sig_DiffPlot", "Adjusted P-Value cut-off", 0.05, min = 0, max = NA),
                                                            numericInput("FC_DiffPlot", "Fold Change cut-off", 1, min = 0, max = NA),
                                                            numericInput("TopN_DiffPlot", "Number of Top points to Label", 10, min = 0, max = NA),
                                                            actionButton("Go_DiffPlot", "Update Plot"),
                                                            plotOutput("DiffPlot", height = "600px"),
                                                            downloadButton("downloadDiffPlot", "Download Volcano Plot")),
                                                   tabPanel("Pathways",
                                                            #Error Message text output
                                                            htmlOutput("Diff_UpRegPath_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Diff_UpRegPath_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("Diff_UpRegPath_Title"),
                                                            plotlyOutput("Diff_UpRegPath_Plot"),
                                                            tags$hr(),
                                                            uiOutput("Diff_DownRegPath_Title"),
                                                            plotlyOutput("Diff_DownRegPath_Plot"),
                                                            tags$hr(),
                                                            uiOutput("Diff_PValPath_Title"),
                                                            plotlyOutput("Diff_PValPath_Plot")))
                                     ))),
                            tabPanel(title ="Correlation Analysis",
                                     id = "Corr",
                                     titlePanel("Correlation Analysis"),
                                     sidebarLayout(
                                       sidebarPanel(
                                         #size of sidebar
                                         width = 3,
                                         checkboxInput("Corr_OneTable", "Is all your data in one table", TRUE),
                                         #File upload code
                                         conditionalPanel(
                                           condition = "input.Corr_OneTable == true",
                                           fileInput("Corr_file", "Choose file to upload",
                                                     multiple = FALSE,
                                                     accept = c("text/csv", 
                                                              "text/comma-separated-values,text/plain",
                                                              ".csv")),
                                           #Check box to see if file has headers 
                                           checkboxInput("header_Corr", "Header", TRUE),
                                           
                                           #input--Separator for the file 
                                           radioButtons("sep_Corr", "Separator",
                                                        choices = c(Comma = ",",
                                                                    Tab = "\t")),
                                           #input---Location of samples
                                           radioButtons("samp_Corr", "Sample location",
                                                        choices = c(Row = "row_Corr",
                                                                    Column = "col_Corr"),
                                                        selected = "row_Corr"),
                                           
                                           #line separator 
                                           tags$hr(),
                                           
                                           #input---condition information column
                                           selectInput("Sample_Corr", "Select Sample Column Name", ""),
                                           
                                           #input---condition information column
                                           selectInput("Cond_Corr", "Select Condition Info Name", ""),
                                           
                                           #line separator
                                           tags$hr(),
                                           
                                           #Group A start column 
                                           selectInput("GrpAstart_Corr", "Row/Column Name Group A Starts", ""),
                                           
                                           #Group A start column 
                                           selectInput("GrpAend_Corr", "Row/Column Name Group A Ends", ""),
                                           
                                           #Group A start column 
                                           selectInput("GrpBstart_Corr", "Row/Column Name Group B Starts", ""),
                                           
                                           #Group A start column 
                                           selectInput("GrpBend_Corr", "Row/Column Name Group B Ends", "")), 
                                           
                                          
                                           #line separator
                                           tags$hr(),
                                           
                                         
                                         conditionalPanel(
                                           condition = "input.Corr_OneTable == false",
                                           
                                           h3("Group A File"),
                                           
                                           #Upload for Group A file 
                                           fileInput("Corr_fileA", "Upload File for Group A",
                                                     multiple = FALSE,
                                                     accept = c("text/csv", 
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           #Check box to see if Group A file has headers 
                                           checkboxInput("header_CorrA", "Header for Group A", TRUE),
                                           
                                           #input--Separator for Group A file 
                                           radioButtons("sep_CorrA", "Separator for Group A",
                                                        choices = c(Comma = ",",
                                                                    Tab = "\t")),
                                           #input---Location of samples for Group A
                                           radioButtons("samp_CorrA", "Sample location for Group A",
                                                        choices = c(Row = "row_CorrA",
                                                                    Column = "col_CorrA"),
                                                        selected = "row_CorrA"),
                                           
                                           #line separator 
                                           tags$hr(),
                                           
                                           #input---condition information column number fir Group A 
                                           selectInput("Cond_CorrA", "Select Condition Info Name", ""),
                                           
                                           #line separator 
                                           tags$hr(),
                                           
                                           h3("Group B File"),
                                           
                                           #Upload for Group B file 
                                           fileInput("Corr_fileB", "Upload File for Group B",
                                                     multiple = FALSE,
                                                     accept = c("text/csv", 
                                                                "text/comma-separated-values,text/plain",
                                                                ".csv")),
                                           #Check box to see if Group B file has headers 
                                           checkboxInput("header_CorrB", "Header for Group B", TRUE),
                                           
                                           #input--Separator for Group B file 
                                           radioButtons("sep_CorrB", "Separator for Group B",
                                                        choices = c(Comma = ",",
                                                                    Tab = "\t")),
                                           #input---Location of samples for Group B
                                           radioButtons("samp_CorrB", "Sample location for Group B",
                                                        choices = c(Row = "row_CorrB",
                                                                    Column = "col_CorrB"),
                                                        selected = "row_CorrB"),
                                           
                                           #line separator 
                                           tags$hr(),
                                           
                                           #input---condition information column number fir Group A 
                                           selectInput("Cond_CorrB", "Select Condition Info Name", "")
                                           ),
                                         
                                         #calculate button
                                         actionButton("Go_Corr", "Calculate")
                                       ),
                                       mainPanel(
                                         width = 9,
                                         height = 12,
                                         tabsetPanel(type="tabs", id ="Corr_tabs",
                                                     tabPanel("Correlation Table",
                                                              #Error Message text output
                                                              htmlOutput("CorrTable_Error"),
                                                              #styling text output to be red and bigger
                                                              tags$head(tags$style("#CorrTable_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                              selectInput("Corr_Table", "Select Correlation Table", ""),
                                                              actionButton("Go_CorrTable", "View Table"),
                                                              dataTableOutput("CorrTable"),
                                                              downloadButton("DownloadCorrTable", "Download Data")),
                                                     tabPanel("Correlation Plots", 
                                                              #Error Message text output
                                                              htmlOutput("CorrPlot_Error"),
                                                              #styling text output to be red and bigger
                                                              tags$head(tags$style("#CorrPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                              selectInput("Corr_Plot", "Select Correlation plot", ""),
                                                              numericInput("Sig_CorrPlot", "Significance cut-off", 0.1, min = 1, max = NA),
                                                              actionButton("Go_CorrPlot", "View Plot"),
                                                              plotOutput("CorrPlot", height = '650px', width = '700px'),
                                                              downloadButton("DownloadCorrPlot", "Download Plot"))))
                                                     ))
                                     ),
                 navbarMenu("ML Models", 
                            tabPanel(title ="Machine Learning ",
                                     id = "ML",
                                     titlePanel("Machine Learning"),
                                     
                                     # Sidebar panel for inputs ----
                                     sidebarPanel(
                                       
                                       width = 3,
                                       
                                       # Input: Select a file ----
                                       
                                       h3("Lipid mediator data"),
                                       
                                       fileInput("LM_ML_File", "Lipid Mediator profiling file:",
                                                 multiple = FALSE,
                                                 accept = c("text/csv",
                                                            "text/comma-separated-values,text/plain",
                                                            ".csv", ".tsv", ".txt")),
                                       
                                       # Input: Select separator ----
                                       radioButtons("sep_LM_ML", "Separator",
                                                    choices = c(Comma = ",",
                                                                Semicolon = ";",
                                                                Tab = "\t"),
                                                    selected = ","),
                                       
                                       #input---Location of samples
                                       radioButtons("samp_ML", "Sample location",
                                                    choices = c(Row = "row_ML",
                                                                Column = "col_ML"),
                                                    selected = "row_ML"),
                                       
                                       #input -- where the sample group info column is
                                       selectInput("Group_ML", "Select Group Column Name", "Please upload file"),
                                       
                                       # Help text ----
                                       helpText(style="text-align: justify;", paste("The LM profiling file consist in a table with samples per row and each 
              lipid mediator as columns. An extra column is added called "), HTML(paste0("<b>","groups","</b>")), 
                                                " specifying to which group every sample belongs to; and an extra row, added after the lipid mediators 
             names, meaning the", HTML(paste0("<b>","second row","</b>")), ", which contains the fatty acid substrates 
             to which every lipid mediator comes from.", sep = ""),
                                       helpText("You can see the format dowloading the example file."),
                                       
                                       
                                       # Download button for the example file ----
                                       downloadButton("lm_example.tsv", "Example LM file"),
                                       
                                       # Horizontal line ----
                                       tags$hr(),
                                       
                                       # Input: Checkbox if you want to run machine learning in clinical data ----
                                       
                                       
                                       
                                       # Horizontal line ----
                                       tags$hr(),
                                       
                                       checkboxGroupInput("Model_ML", "Select Machine Learning Model",
                                                          c("Random Forests" = "RF_ML", 
                                                            "Support Vector Machine" = "SVM_ML",
                                                            "Bayseian Classifier" = "BC_ML",
                                                            "Elastic Net Regresion" = "LA_ML", 
                                                            "Extreme Gradient Boosting" = "XGB_ML")),
                                       
                                       # Action bottom to create and run the ML models ---
                                       actionButton("Go_ML","Create Machine Learning Models")
                                     ),
                                     # Main panel for displaying outputs ----
                                     mainPanel(
                                       width = 9,
                                       height = 12,
                                       tabsetPanel(type = "tabs", id = "ML_tabs",
                                                   tabPanel("Accuracy",
                                                            #Error Message text output
                                                            htmlOutput("Accuracy_ML_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Accuracy_ML_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("AccPlot_Title"),
                                                            downloadButton("downloadAcc_ML_Plot", "Download Accuracy Plot"),
                                                            downloadButton("downloadAcc_ML_Table", "Download Accuracy Table"),
                                                            tags$hr(),
                                                            plotOutput("Accuracy_ML"),
                                                            tags$hr(),
                                                            uiOutput("AccTable_Title"),
                                                            tags$hr(),
                                                            dataTableOutput("ML_Table")),
                                                   tabPanel("Random Forest",
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Parameters",
                                                                                 #Error Message text output
                                                                                 htmlOutput("RF_Plot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#RF_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("RF_Plot_Title"),
                                                                                 tags$hr(),
                                                                                 plotOutput("RF_Plot", height = "1500px"),
                                                                                 selectInput("RF_Mod_Num1", "Select which model you want to download", ""),
                                                                                 downloadButton("download_RF_Plot", "Download Random Forest Plot"),
                                                                                 downloadButton("download_RF_Mod1", "Download Random Forest Model"),
                                                                                 tags$hr()),
                                                                        tabPanel("Importance",
                                                                                 #Error Message text output
                                                                                 htmlOutput("RF_VIPPlot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#RF_VIPPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("RF_VIPPlot_Title"),
                                                                                 tags$hr(),
                                                                                 plotOutput("RF_VIPPlot", height = "2000px"),
                                                                                 selectInput("RF_Mod_Num", "Select which plot/model you want to download", ""),
                                                                                 downloadButton("download_RF_VIPPlot", "Download Random Forest VIP Plot"),
                                                                                 downloadButton("download_RF_Mod", "Download Random Forest Model"),
                                                                                 tags$hr()))
                                                   ),
                                                            
                                                   tabPanel("Extreme Gradient Boosting",
                                                            tabsetPanel(
                                                              type = "tabs",
                                                              tabPanel("Parameters", 
                                                                       #Error Message text output
                                                                       htmlOutput("XGB_Models_Error"),
                                                                       #styling text output to be red and bigger
                                                                       tags$head(tags$style("#XGB_Models_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                       uiOutput("XGB_Plot_Title"),
                                                                       tags$hr(), dataTableOutput("XGB_Models"),
                                                                       tags$hr(), plotlyOutput("XGB_Plot1"),
                                                                       tags$hr(), plotlyOutput("XGB_Plot2"),
                                                                       tags$hr(), plotlyOutput("XGB_Plot3"),
                                                                       tags$hr(), plotlyOutput("XGB_Plot4"),
                                                                       tags$hr(), plotlyOutput("XGB_Plot5"),
                                                                       selectInput("XGB_Mod_Num1", "Select which Model you want to download", ""),
                                                                       downloadButton("download_XGB_Table", "Download XGBoost Table"),
                                                                       downloadButton("download_XGB_Mod1", "Download Extreme Gradien Boosting Model"),
                                                                       tags$hr()),
                                                              tabPanel("Importance Plots", 
                                                                       #Error Message text output
                                                                       htmlOutput("XGB_VIPPlot_Error"),
                                                                       #styling text output to be red and bigger
                                                                       tags$head(tags$style("#XGB_VIPPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                       uiOutput("XGB_VIPPlot_Title"),
                                                                       plotOutput("XGB_VIPPlot", height = "2000px"),
                                                                       tags$hr(),
                                                                       selectInput("XGB_Mod_Num", "Select which plot/Model you want to download", ""),
                                                                       downloadButton("download_XGB_VIPPlot", "Download XGBoost VIP Plot"),
                                                                       downloadButton("download_XGB_Mod", "Download Extreme Gradien Boosting Model"),
                                                                       tags$hr()))
                                                              ),
                                                   tabPanel("SVM",
                                                            #Error Message text output
                                                            htmlOutput("SVM_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#SVM_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("SVM_Plot_Title"),
                                                            tags$hr(),
                                                            plotOutput("SVM_Plot", height = "1000px"),
                                                            selectInput("SVM_Mod_Num", "Select which Model you want to download", ""),
                                                            downloadButton("download_SVM_Plot", "Download SVM Plot"),
                                                            downloadButton("download_SVM_Mod", "Download SVM Model"),
                                                            tags$hr()),
                                                   tabPanel("Elastic Net Regression",
                                                            #Error Message text output
                                                            htmlOutput("LA_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#LA_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("LA_Plot_Title"),
                                                            tags$hr(),
                                                            plotOutput("LA_Plot", height = "1000px"),
                                                            selectInput("LA_Mod_Num", "Select which Model you want to download", ""),
                                                            downloadButton("download_LA_Plot", "Download Elastic Net Plot"),
                                                            downloadButton("download_LA_Mod", "Download Elastic Net Model"),
                                                            tags$hr()),
                                                   tabPanel("Bayesian Linear Model",
                                                            #Error Message text output
                                                            htmlOutput("BC_Mod_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#BC_Mod_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            selectInput("BC_Mod_Num", "Select which Model you want to download", ""),
                                                            downloadButton("download_BC_Mod", "Download Bayes GLM Model"))
                                                  
                                       ))),
                            
                            tabPanel(title ="Optimize ML Model",
                                     id = "ML2",
                                     titlePanel("Optimize Machine Learning Model"),
                                     sidebarPanel(
                                       
                                       width = 3,
                                       
                                       # Input: Select a file ----
                                       
                                       h3("Lipid mediator data"),
                                       
                                       fileInput("LM_ML_File2", "Lipid Mediator profiling file:",
                                                 multiple = FALSE,
                                                 accept = c("text/csv",
                                                            "text/comma-separated-values,text/plain",
                                                            ".csv", ".tsv", ".txt")),
                                       
                                       # Input: Select separator ----
                                       radioButtons("sep_LM_ML2", "Separator",
                                                    choices = c(Comma = ",",
                                                                Semicolon = ";",
                                                                Tab = "\t"),
                                                    selected = ","),
                                       #input---Location of samples
                                       radioButtons("samp_ML2", "Sample location",
                                                    choices = c(Row = "row_ML2",
                                                                Column = "col_ML2"),
                                                    selected = "row_ML2"),
                                       
                                       #input -- where the sample group info column is
                                       selectInput("Group_ML2", "Select Group Column Name", "Please upload file"),
                                       
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a random forest model 
                                       checkboxInput("RF_Build", "Do you want to build a Random Forest Model?", FALSE),
                                       
                                       conditionalPanel(
                                         condition = "input.RF_Build == true",
                                         #select number of trees
                                         numericInput("RF_ntrees", "Select Number of Trees", 10000, min = 1, max = NA)),
                                       
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a XGBoost model 
                                       checkboxInput("XGB_Build", "Do you want to build a Extreme Gradient Boosting Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.XGB_Build == true",
                                         #select number of rounds 
                                         numericInput("XGB_nrounds", "Select Number of Rounds", 10000, min = 1, max = NA),
                                         
                                         #select model for All LM
                                         numericInput("XGB_ALL_LM_Model", "Select Model for All LM.", 1, min = 1, max = 5),
                                         
                                         #select model for DHA
                                         checkboxInput("XGB_Build_DHA", "Does Data contain DHA SPMs?", FALSE),
                                         conditionalPanel( 
                                           condition = "input.XGB_Build_DHA == true",
                                           numericInput("XGB_DHA_Model", "Select Model for DHA.", 1, min = 1, max = 5)),
                                         
                                         #Select model for n3-DPA
                                         checkboxInput("XGB_Build_n3DPA", "Does Data contain n3-DPA SPMs?", FALSE),
                                         conditionalPanel( 
                                           condition = "input.XGB_Build_n3DPA == true",
                                           numericInput("XGB_n3DPA_Model", "Select Model for n3-DPA.", 1, min = 1, max = 5)),
                                         
                                         #Select model for EPA
                                         checkboxInput("XGB_Build_EPA", "Does Data contain EPA SPMs?", FALSE),
                                         conditionalPanel( 
                                           condition = "input.XGB_Build_EPA == true",
                                           numericInput("XGB_EPA_Model", "Select Model for EPA.", 1, min = 1, max = 5)),
                                         
                                         #Select model for n3-DPA
                                         checkboxInput("XGB_Build_AA", "Does Data contain AA SPMs?", FALSE),
                                         conditionalPanel( 
                                           condition = "input.XGB_Build_AA == true",
                                           numericInput("XGB_AA_Model", "Select Model for AA.", 1, min = 1, max = 5))
                                       ),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a support vector mechanism model 
                                       checkboxInput("SVM_Build", "Do you want to build a Support Vector Machine Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.SVM_Build == true",
                                         #select number of trees
                                         numericInput("SVM_Ensemble", "Select Number of Ensembles", 70, min = 1, max = NA),
                                         numericInput("SVM_BootNum", "Select Number of Boostraps", 70, min = 1, max = NA)),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a Lasso model 
                                       checkboxInput("LA_Build", "Do you want to build a Elastic Net Regressin Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.LA_Build == true",
                                         #select number of trees
                                         numericInput("LA_BootNum", "Select Number of Bootstraps", 70, min = 1, max = NA)),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a Bayseian model 
                                       checkboxInput("BC_Build", "Do you want to build a Baysian Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.BC_Build == true",
                                         #select number of trees
                                         numericInput("BC_BootNum", "Select Number of Bootstraps", 70, min = 1, max = NA)),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       # Action bottom to create and run the ML models ---
                                       actionButton("Build_ML","Build Machine Learning Models")
                                     ),
                                     mainPanel(
                                       width = 9,
                                       tabsetPanel(type ="tabs", id = "build_tabs",
                                                   tabPanel("Accuracy",
                                                            #Error Message text output
                                                            htmlOutput("Build_Accuracy_ML_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Build_Accuracy_ML_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("Build_AccPlot_Title"),
                                                            downloadButton("downloadBuild_Acc_ML_Plot", "Download Accuracy Plot"),
                                                            downloadButton("downloadBuild_Acc_ML_Table", "Download Accuracy Table"),
                                                            tags$hr(),
                                                            plotOutput("Build_Accuracy_ML"),
                                                            tags$hr(),
                                                            uiOutput("Build_AccTable_Title"),
                                                            tags$hr(),
                                                            dataTableOutput("Build_ML_Table")),
                                                   tabPanel("Random Forest", 
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Parameters",
                                                                                 #Error Message text output
                                                                                 htmlOutput("Build_RF_Plot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#Build_RF_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("Build_RF_Plot_Title"),
                                                                                 tags$hr(),
                                                                                 selectInput("Build_RF_Mod_Num1", "Select which model you want to download", ""),
                                                                                 downloadButton("downloadBuild_RF_Plot", "Download Random Forest Plot"),
                                                                                 downloadButton("downloadBuild_RF_Mod1", "Download Random Forest Model"),
                                                                                 plotOutput("Build_RF_Plot", height = "1500px")),
                                                                        tabPanel("Importance",
                                                                                 #Error Message text output
                                                                                 htmlOutput("Build_RF_VIPPlot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#Build_RF_VIPPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("Build_RF_VIPPlot_Title"),
                                                                                 tags$hr(),
                                                                                 selectInput("Build_RF_Mod_Num", "Select which model you want to download", ""),
                                                                                 downloadButton("downloadBuild_RF_Mod", "Download Random Forest Model"),
                                                                                 downloadButton("downloadBuild_RF_VIPPlot", "Download Random Forest VIP Plot"),
                                                                                 plotOutput("Build_RF_VIPPlot", height = "2000px")))),
                                                   tabPanel("Extreme Gradient Boosting", 
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Parameters",
                                                                                 #Error Message text output
                                                                                 htmlOutput("Build_XGB_Models_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#Build_XGB_Models_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("Build_XGB_Plot_Title"),
                                                                                 tags$hr(),
                                                                                 selectInput("Build_XGB_Mod_Num1", "Select which model you want to download", ""),
                                                                                 downloadButton("downloadBuild_XGB_Table", "Download XGBoost Table"),
                                                                                 downloadButton("downloadBuild_XGB_Mod1", "Download XGBoost Model"),
                                                                                 dataTableOutput("XGB_Build_Models"),
                                                                                 plotlyOutput("Build_XGB_Plot1"),
                                                                                 plotlyOutput("Build_XGB_Plot2"),
                                                                                 plotlyOutput("Build_XGB_Plot3"),
                                                                                 plotlyOutput("Build_XGB_Plot4"),
                                                                                 plotlyOutput("Build_XGB_Plot5")
                                                                        ),
                                                                        tabPanel("Importance",
                                                                                 #Error Message text output
                                                                                 htmlOutput("Build_XGB_VIPPlot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#Build_XGB_VIPPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("Build_XGB_VIPPlot_Title"),
                                                                                 selectInput("Build_XGB_Mod_Num", "Select which model you want to download", ""),
                                                                                 downloadButton("downloadBuild_XGB_VIPPlot", "Download XGBoost VIP Plot"),
                                                                                 downloadButton("downloadBuild_XGB_Mod", "Download XGBoost Model"),
                                                                                 tags$hr(),
                                                                                 plotOutput("Build_XGB_VIPPlot", height = "2000px")
                                                                        ))),
                                                   tabPanel("SVM",
                                                            #Error Message text output
                                                            htmlOutput("Build_SVM_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Build_SVM_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("Build_SVM_Plot_Title"),
                                                            tags$hr(),
                                                            selectInput("Build_SVM_Mod_Num", "Select which model you want to download", ""),
                                                            downloadButton("downloadBuild_SVM_Plot", "Download SVM Plot"),
                                                            downloadButton("downloadBuild_SVM_Mod", "Download SVM Model"),
                                                            plotOutput("Build_SVM_Plot", height = "1500px")),
                                                   tabPanel("Elastic Net Regression",
                                                            #Error Message text output
                                                            htmlOutput("Build_LA_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Build_LA_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("Build_LA_Plot_Title"),
                                                            tags$hr(),
                                                            selectInput("Build_LA_Mod_Num", "Select which model you want to download", ""),
                                                            downloadButton("downloadBuild_LA_Plot", "Download Elastic Net Regression Plot"),
                                                            downloadButton("downloadBuild_LA_Mod", "Download Elastic Net Regression Model"),
                                                            tags$hr(),
                                                            plotOutput("Build_LA_Plot", height = "1500px")),
                                                   tabPanel("Bayes",
                                                            #Error Message text output
                                                            htmlOutput("Build_BC_Mod_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Build_BC_Mod_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            selectInput("Build_BC_Mod_Num", "Select which model you want to download", ""),
                                                            downloadButton("downloadBuild_BC_Mod", "Download Bayes GLM Model"))
                                                   ))),
                            
                            tabPanel(title ="Build ML Model",
                                     id = "ML2_5",
                                     titlePanel("Build Machine Learning Model"),
                                     sidebarPanel(
                                       
                                       width = 3,
                                       
                                       # Input: Select a file ----
                                       
                                       h3("Lipid mediator data"),
                                       
                                       fileInput("LM_ML_File2_5", "Lipid Mediator profiling file:",
                                                 multiple = FALSE,
                                                 accept = c("text/csv",
                                                            "text/comma-separated-values,text/plain",
                                                            ".csv", ".tsv", ".txt")),
                                       
                                       # Input: Select separator ----
                                       radioButtons("sep_LM_ML2_5", "Separator",
                                                    choices = c(Comma = ",",
                                                                Semicolon = ";",
                                                                Tab = "\t"),
                                                    selected = ","),
                                       
                                       #input---Location of samples
                                       radioButtons("samp_ML2_5", "Sample location",
                                                    choices = c(Row = "row_ML2_5",
                                                                Column = "col_ML2_5"),
                                                    selected = "row_ML2_5"),
                                       
                                       #input -- where the sample group info column is
                                       selectInput("Group_ML2_5", "Select Group Column Name", "Please upload file"),
                                       
                                       #selecting mediators of interest 
                                       selectInput("MetName_ML", "Select all Mediators of intrest", "Please upload file", multiple = TRUE, selectize = TRUE),
                                       
                                       tags$hr(),
                                       
                                       #option to build a random forest model 
                                       checkboxInput("RF_BuildMet", "Do you want to build a Random Forest Model?", FALSE),
                                       
                                       #option to build a XGBoost model 
                                       checkboxInput("XGB_BuildMet", "Do you want to build a Extreme Gradient Boosting Model?", FALSE),
                                      
                                       #option to build a support vector mechanism model 
                                       checkboxInput("SVM_BuildMet", "Do you want to build a Support Vector Machine Model?", FALSE),
                                       
                                       #option to build a Lasso model 
                                       checkboxInput("LA_BuildMet", "Do you want to build a Elastic Net Regressin Model?", FALSE),
                                      
                                       #option to build a Bayseian model 
                                       checkboxInput("BC_BuildMet", "Do you want to build a Baysian Model?", FALSE),
                                       
                                       tags$hr(),
                                       
                                       # Action bottom to create and run the ML models ---
                                       actionButton("Build_ML2_5","Build Machine Learning Models")
                                     ),
                                     mainPanel(
                                       width = 9,
                                       tabsetPanel(type ="tabs", id = "BuildMet_tabs",
                                                   tabPanel("Accuracy",
                                                            #Error Message text output
                                                            htmlOutput("BuildMet_Accuracy_ML_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#BuildMet_Accuracy_ML_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("BuildMet_AccPlot_Title"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_AccPlot", "Download Accuracy Plot"),
                                                            downloadButton("downloadBuildMet_Acc_ML_Table", "Download Accuracy Table"),
                                                            plotOutput("BuildMet_AccPlot"),
                                                            dataTableOutput("BuildMet_ML_Table")),
                                                   tabPanel("Random Forest", 
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Parameters",
                                                                                 #Error Message text output
                                                                                 htmlOutput("BuildMet_RF_Plot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#BuildMet_RF_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("BuildMet_RF_Plot_Title"),
                                                                                 tags$hr(),
                                                                                 downloadButton("downloadBuildMet_RF_Plot", "Download Random Forest Plot"),
                                                                                 downloadButton("downloadBuildMet_RF_Mod", "Download Random Forest Model"),
                                                                                 plotOutput("BuildMet_RF_Plot", height = "500px")),
                                                                        tabPanel("Importance",
                                                                                 #Error Message text output
                                                                                 htmlOutput("BuildMet_RF_VIPPlot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#BuildMet_RF_VIPPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("BuildMet_RF_VIPPlot_Title"),
                                                                                 tags$hr(),
                                                                                 downloadButton("downloadBuildMet_RF_VIPPlot", "Download Random Forest VIP Plot"),
                                                                                 downloadButton("downloadBuildMet_RF_Mod1", "Download Random Forest Model"),
                                                                                 plotOutput("BuildMet_RF_VIPPlot", height = "500px")))),
                                                   tabPanel("Extreme Gradient Boosting", 
                                                            tabsetPanel(type = "tabs",
                                                                        tabPanel("Parameters",
                                                                                 #Error Message text output
                                                                                 htmlOutput("BuildMet_XGB_Models_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#BuildMet_XGB_Models_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("BuildMet_XGB_Plot_Title"),
                                                                                 tags$hr(),
                                                                                 downloadButton("downloadBuildMet_XGB_Mod", "Download Extreme Gradient Boosting Model"),
                                                                                 plotlyOutput("BuildMet_XGB_Plot")),
                                                                        tabPanel("Importance",
                                                                                 #Error Message text output
                                                                                 htmlOutput("BuildMet_XGB_VIPPlot_Error"),
                                                                                 #styling text output to be red and bigger
                                                                                 tags$head(tags$style("#BuildMet_XGB_VIPPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                                                 uiOutput("BuildMet_XGB_VIPPlot_Title"),
                                                                                 downloadButton("downloadBuildMet_XGB_VIPPlot", "Download XGBoost VIP Plot"),
                                                                                 downloadButton("downloadBuildMet_XGB_Mod1", "Download Extreme Gradient Boosting Model"),
                                                                                 tags$hr(),
                                                                                 plotOutput("BuildMet_XGB_VIPPlot", height = "500px")
                                                                        ))),
                                                   tabPanel("SVM",
                                                            #Error Message text output
                                                            htmlOutput("BuildMet_SVM_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#BuildMet_SVM_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("BuildMet_SVM_Plot_Title"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_SVM_Plot", "Download SVM Plot"),
                                                            downloadButton("downloadBuildMet_SVM_Mod", "Download SVM Model"),
                                                            plotOutput("BuildMet_SVM_Plot", height = "500px")),
                                                   tabPanel("Elastic Net Regression",
                                                            #Error Message text output
                                                            htmlOutput("BuildMet_LA_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#BuildMet_LA_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            uiOutput("BuildMet_LA_Plot_Title"),
                                                            tags$hr(),
                                                            downloadButton("downloadBuildMet_LA_Plot", "Download Elastic Net Regression Plot"),
                                                            downloadButton("downloadBuildMet_LA_Mod", "Download Elastic Net Model"),
                                                            tags$hr(),
                                                            plotOutput("BuildMet_LA_Plot", height = "500px")),
                                                   tabPanel("Bayes GLM",
                                                            #Error Message text output
                                                            htmlOutput("BuildMet_BC_Mod_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#BuildMet_BC_Mod_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            downloadButton("downloadBuildMet_BC_Mod", "Download Bayes GLM Model"))
                                       ))),
                            
                            tabPanel(title ="Run ML Model",
                                     id = "ML3",
                                     titlePanel("Run Machine Learning Model"),
                                     sidebarPanel(
                                       
                                       width = 3,
                                       
                                       # Input: Select a file ----
                                       
                                       h3("Lipid mediator data"),

                                       fileInput("LM_ML_File3", "Lipid Mediator profiling file:",
                                                 multiple = FALSE,
                                                 accept = c("text/csv",
                                                            "text/comma-separated-values,text/plain",
                                                            ".csv", ".tsv", ".txt")),
                                       
                                       # actionButton("filechoose=",label = "Pick a file")
                                       
                                       # Input: Select separator ----
                                       radioButtons("sep_LM_ML3", "Separator",
                                                    choices = c(Comma = ",",
                                                                Semicolon = ";",
                                                                Tab = "\t"),
                                                    selected = ","),
                                       
                                       #input---Location of samples
                                       radioButtons("samp_ML3", "Sample location",
                                                    choices = c(Row = "row_ML3",
                                                                Column = "col_ML3"),
                                                    selected = "row_ML3"),
                                       
                                       #input -- where the sample group info column is
                                       selectInput("Group_ML3", "Select Group Column Name", "Please upload file"),
                                       
                                       
                                       #checkbox input
                                       checkboxInput("Choose_MetName3", "Do you want to filter the mediators included?", FALSE),
                                       
                                       #conditional panel
                                       conditionalPanel(
                                         condition = "input.Choose_MetName3 == true",
                                         #selecting metabolites of interest 
                                         selectInput("MetName_ML3", "Select all Mediators of intrest", "Please upload file", multiple = TRUE, selectize = TRUE)
                                         
                                       ),
                                       
                                       
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a random forest model 
                                       checkboxInput("RF_Run", "Do you want to build a Random Forest Model?", FALSE),
                                       
                                       #conditional panel for random forest model 
                                       conditionalPanel(
                                         condition = "input.RF_Run == true",
                                         numericInput('numRF_Model', 'select number of random Forest models', 1, min = 1, max = 5),  
                                         fileInput("RF_Mod_File1", "Random Forest Model:",
                                                   multiple = FALSE,
                                                   accept = ".rds"),
                                        
                                        #conditional panels for each additional model 
                                        conditionalPanel(
                                          condition = "input.numRF_Model > '1'",
                                          fileInput("RF_Mod_File2", "Random Forest Model:",
                                                    multiple = FALSE,
                                                    accept = ".rds")),
                                        conditionalPanel(
                                          condition = "input.numRF_Model > '2'",
                                          fileInput("RF_Mod_File3", "Random Forest Model:",
                                                    multiple = FALSE,
                                                    accept = ".rds")),
                                        conditionalPanel(
                                          condition = "input.numRF_Model > '3'",
                                          fileInput("RF_Mod_File4", "Random Forest Model:",
                                                    multiple = FALSE,
                                                    accept = ".rds")),
                                        conditionalPanel(
                                          condition = "input.numRF_Model > '4'",
                                          fileInput("RF_Mod_File5", "Random Forest Model:",
                                                    multiple = FALSE,
                                                    accept = ".rds"))
                                        ),
                                       
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a XGBoost model 
                                       checkboxInput("XGB_Run", "Do you want to Run an Extreme Gradient Boosting Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.XGB_Run == true",
                                         
                                         #input to see how many models of XGboost to test
                                         numericInput('numXGB_Model', 'select number of XGBoost models', 1, min = 1, max = 5),
                                         
                                         #input for the first model 
                                         fileInput("XGB_Mod_File1", "XGBoost Model:",
                                                   multiple = FALSE,
                                                   accept = ".rds"),
                                         
                                         #conditional panels for each additional model 
                                         conditionalPanel(
                                           condition = "input.numXGB_Model > '1'",
                                           fileInput("XGB_Mod_File2", "XGBoost Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numXGB_Model > '2'",
                                           fileInput("XGB_Mod_File3", "XGBoost Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numXGB_Model > '3'",
                                           fileInput("XGB_Mod_File4", "XGBoost Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         
                                         conditionalPanel(
                                           condition = "input.numXGB_Model > '4'",
                                           fileInput("XGB_Mod_File5", "XGBoost Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds"))
                                         ),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a support vector mechanism model 
                                       checkboxInput("SVM_Run", "Do you want to Run a Support Vector Machine Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.SVM_Run == true",
                                         #number of SVM models to run 
                                         numericInput('numSVM_Model', 'select number of SVM models', 1, min = 1, max = 5),
                                         fileInput("SVM_Mod_File1", "SVM Model:",
                                                   multiple = FALSE,
                                                   accept = ".rds"),
                                         
                                         #conditional panels for each additional model 
                                         conditionalPanel(
                                           condition = "input.numSVM_Model > '1'",
                                           fileInput("SVM_Mod_File2", "SVM Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numSVM_Model > '2'",
                                           fileInput("SVM_Mod_File3", "SVM Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numSVM_Model > '3'",
                                           fileInput("SVM_Mod_File4", "SVM Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numSVM_Model > '4'",
                                           fileInput("SVM_Mod_File5", "SVM Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds"))
                                         ),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to Run elastic Net regresssion Model  
                                       checkboxInput("LA_Run", "Do you want to Run a Elastic Net Regressin Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.LA_Run == true",
                                         #input to select how many Elastic Net Models to Run
                                         numericInput('numLA_Model', 'select number of Elastic Net Regression models', 1, min = 1, max = 5),
                                         fileInput("LA_Mod_File1", "Elastic Net Regression Model:",
                                         multiple = FALSE,
                                         accept = ".rds"),
                                         
                                         #conditional panels for each additional model 
                                         conditionalPanel(
                                           condition = "input.numLA_Model > '1'",
                                           fileInput("LA_Mod_File2", "Elastic Net Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numLA_Model > '2'",
                                           fileInput("LA_Mod_File3", "Elastic Net Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numLA_Model > '3'",
                                           fileInput("LA_Mod_File4", "Elastic Net Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numLA_Model > '4'",
                                           fileInput("LA_Mod_File5", "Elastic Net Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds"))
                                         ),
                                       
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       #option to build a Bayseian model 
                                       checkboxInput("BC_Run", "Do you want to run a Baysian Model?", FALSE),
                                       conditionalPanel(
                                         condition = "input.BC_Run == true",
                                         
                                         #input to select how many Bayes Regression Models to Run
                                         numericInput('numBC_Model', 'select number of Bayes GLM models', 1, min = 1, max = 5),
                                         
                                         #upload BC model 
                                         fileInput("BC_Mod_File1", "Bayes Regression Model:",
                                                   multiple = FALSE,
                                                   accept = ".rds"),
                                         
                                         #conditional panels for each additional model 
                                         conditionalPanel(
                                           condition = "input.numBC_Model > '1'",
                                           fileInput("BC_Mod_File2", "Bayes Regression Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numBC_Model > '2'",
                                           fileInput("BC_Mod_File3", "Bayes Regression Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numBC_Model > '3'",
                                           fileInput("BC_Mod_File4", "Bayes Regression Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds")),
                                         conditionalPanel(
                                           condition = "input.numBC_Model > '4'",
                                           fileInput("BC_Mod_File5", "Bayes Regression Model:",
                                                     multiple = FALSE,
                                                     accept = ".rds"))
                                         
                                         ),
                                       # Horizontal line ----
                                       tags$hr(style = "border-top: 1px solid #000000;"),
                                       
                                       # Action bottom to create and run the ML models ---
                                       actionButton("Run_ML","Run Machine Learning Models")
                                     ),
                                     mainPanel(
                                       width = 9,
                                       tabsetPanel(type ="tabs", id = "Run_tabs",
                                                   tabPanel("Accuracy",
                                                            htmlOutput("ROC_Plot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#ROC_Plot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            plotOutput("ROC_Plot", height = "600px"),
                                                            dataTableOutput("ROC_Table"),
                                                            downloadButton("downloadRun_ROCPlot", "Download Accuracy Plot"),
                                                            downloadButton("downloadRun_ROCTable", "Download Accuracy Table")),
                                                   tabPanel("Random Forest ROC",
                                                            htmlOutput("Run_RF_ROCPlot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Run_RF_ROCPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            downloadButton("downloadRun_RF_ROCPlot", "Download Random Forest ROCPlot"),
                                                            plotOutput("Run_RF_ROCPlot", height = "600px")),
                                                   tabPanel("Extreme Gradient Boosting ROC",
                                                            htmlOutput("Run_XGB_ROCPlot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Run_XGB_ROCPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            downloadButton("downloadRun_XGB_ROCPlot", "Download XGboost ROC Plot"),
                                                            tags$hr(),
                                                            plotOutput("Run_XGB_ROCPlot", height = "600px")),
                                                   tabPanel("SVM ROC",
                                                            htmlOutput("Run_SVM_ROCPlot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Run_SVM_ROCPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            downloadButton("downloadRun_SVM_ROCPlot", "Download SVM ROC Plot"),
                                                            tags$hr(),
                                                            plotOutput("Run_SVM_ROCPlot", height = "600px")),
                                                   tabPanel("Elastic Net Regression ROC",
                                                            htmlOutput("Run_LA_ROCPlot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Run_LA_ROCPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            downloadButton("downloadRun_LA_ROCPlot", "Download Elastic Net ROC Plot"),
                                                            tags$hr(),
                                                            plotOutput("Run_LA_ROCPlot", height = "600px")),
                                                   tabPanel("Bayesian Regression ROC",
                                                            htmlOutput("Run_BC_ROCPlot_Error"),
                                                            #styling text output to be red and bigger
                                                            tags$head(tags$style("#Run_BC_ROCPlot_Error{color: red; font-size: 20px; font-style: italic;}")),
                                                            downloadButton("downloadRun_BC_ROCPlot", "Download Bayes Classifier ROC Plot"),
                                                            tags$hr(),
                                                            plotOutput("Run_BC_ROCPlot", height = "600px"))
                                       )))
                            
                            )
                 )
                 
                 #tabPanel("Network Analysis")
                 

server <- function(input, output, session) {
  
  options(shiny.maxRequestSize=100*1024^2)
  
  observe({
    req(input$Data_File)
    
    options(digits = 3)
    
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {MS_Table<- read.csv(input$Data_File$datapath,
                            header = input$Header_Data,
                            sep = input$sep_Data)
        
        #get list of sample names as options 
        Sample_Names<-unique(as.character(unname(unlist(MS_Table$Sample.Name))))
        
        #if list to update inputs is empty, input has error message
        if (Sample_Names == ''){
          updateSelectInput(session,"Cross_Name", choices = "Error!!Please check file and Separator")
          updateSelectInput(session,"not_Samples", choices = "Error!!Please check file and Separator")
          updateSelectInput(session,"MaxCol_Name", choices = "Error!!Please check file and Separator")
        }
        
        #update inputs 
        updateSelectInput(session, "Cross_Name", choices = Sample_Names) #update cross reference selection
        updateSelectInput(session, "not_Samples", choices = Sample_Names) #update not samples selection
        updateSelectInput(session, 'MaxCol_Name', choices = Sample_Names) #update Max column selection
        
        Samples <-dplyr::select(MS_Table, Sample.Name, Component.Name)
        
        #exclude standards, and cross references 
        for (x in c("std", "at", "iso", "bio", "cross ref")){
          Samples <- Samples[Samples["Sample.Name"] != x, ]
        }
        #exclude internal standards from the metabolites names 
        IS_Names<-c('d85SHETE 116',	'd4LTB4 197',	'd5MaR1 177',	'd5MaR2 177',	'd4PGE2 193',	'd5LXA4 115',	'd5RvD3 147',	'd5RvD2 141',	
                    'd4RvE1 197',	'd517RRvD1 141', 'd5LTE4 194', 'd5LTD4 194', 'd5LTC4 194')
        for (x in IS_Names){
          Samples <- Samples[Samples["Component.Name"] != x, ]
        }
        
        #download sample list 
        Sample_List<-as.data.frame(unique(as.character(unname(unlist(Samples$Sample.Name)))))
        colnames(Sample_List)[1]<-"Sample_Name"
        Sample_List["Weight"]<-1:length(Sample_List$Sample_Name)
        
        output$DownloadSampleList <- downloadHandler(filename = function(){"Sample_List.csv"}, 
                                                     content = function(fname){
                                                       write.csv(Sample_List, fname, row.names = FALSE)})
      
      }, 
      
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check all inputs are correct','',type = "error")
        updateSelectInput(session,"Cross_Name", choices = "Error!!Please check file and Separator")
        updateSelectInput(session,"not_Samples", choices = "Error!!Please check file and Separator")
        updateSelectInput(session,"MaxCol_Name", choices = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
  })
  
  observe({ 
    ##Example file downloads 
    #read file with crossRef their own cross reference 
    MS2_CrossRef<-read.csv(file = "coefficients/MS2_CrossRef.csv",
                           header = TRUE,
                           sep = ",")
    
    #download table 
    output$DownloadCrossRef_File <- downloadHandler(filename = function(){"Example_CrossRef.csv"}, 
                                                    content = function(fname){
                                                      write.csv(MS2_CrossRef, fname, row.names = TRUE)})
    
    #read file with crossRef their own cross reference 
    MS2_ConvCoef<-read.csv(file = "coefficients/MS2_ConvCoef.csv",
                           header = TRUE,
                           sep = ",")
    
    #download table 
    output$DownloadConvCoef_File <- downloadHandler(filename = function(){"Example_ConvCoef.csv"}, 
                                                    content = function(fname){
                                                      write.csv(MS2_ConvCoef, fname, row.names = TRUE)})})
  

  observeEvent(input$Go_Data|input$Filt_Table, {
    req(input$Data_File)
    
    tryCatch(
      ##try to run the code to process the table 
      {
        #form dataframe with file 
        MS_Table<- read.csv(input$Data_File$datapath,
                            header = input$Header_Data,
                            sep = input$sep_Data)
        
        #Step 1: Obtain the Cross-Reference Coefficient ----------------
        #filter to only keep rows of interest 
        MS_Table_Filtered <-dplyr::select(MS_Table, Sample.Name,  Sample.Index, Component.Type, Component.Name,  IS.Name, Area)
        
        #remove hyphens in the names 
        MS_Table_Filtered$Component.Name <- gsub("\\-", "", MS_Table_Filtered$Component.Name)
        MS_Table_Filtered$IS.Name <- gsub("\\-", "", MS_Table_Filtered$IS.Name)
        
        #get samples of interest
        Sample_List <- MS_Table_Filtered
        
        #exclude standards, and cross references 
        for (x in input$not_Samples){
          Sample_List <- Sample_List[Sample_List["Sample.Name"] != x, ]
        }
        
        Sample_Order<-unique(as.character(unname(unlist(Sample_List$Sample.Name))))
        
        #only if the sample has a cross reference coefficient 
        if(input$CrossRef_Incl == TRUE){
          #filter to keep cross reference values for internal standards 
          Cross_Ref <- MS_Table_Filtered[MS_Table_Filtered["Sample.Name"] == input$Cross_Name, ]
          Cross_Ref <- Cross_Ref[Cross_Ref["Component.Type"] == "Internal Standards", ]
          
          #Seperate tables for each of three cross reference values 
          index <- unique(Cross_Ref$Sample.Index)
          Cross_Ref1 <- Cross_Ref[Cross_Ref["Sample.Index"] == index[1], ]
          Cross_Ref2 <- Cross_Ref[Cross_Ref["Sample.Index"] == index[2], ]
          Cross_Ref3 <- Cross_Ref[Cross_Ref["Sample.Index"] == index[3], ]
          
          #only keep columns of interest 
          Cross_Ref1<-dplyr::select(Cross_Ref1,  Component.Name, Area)
          Cross_Ref2<-dplyr::select(Cross_Ref2,  Component.Name, Area)
          Cross_Ref3<-dplyr::select(Cross_Ref3,  Component.Name, Area)
          
          #merge the three tables 
          Cross_Ref_Table <- merge(Cross_Ref1, Cross_Ref2, by= 'Component.Name', all = TRUE)
          Cross_Ref_Table <- merge(Cross_Ref_Table, Cross_Ref3, by= 'Component.Name', all = TRUE)
          
          #get the average for the cross reference value
          Cross_Ref_Table[,-1]<- as.data.frame(sapply(Cross_Ref_Table[,-1], function (x) as.numeric(as.character(x))))
          Mean_CrossRef<-as.data.frame(rowMeans(Cross_Ref_Table[,-1]))
          
          Mean_CrossRef<-cbind(Cross_Ref_Table$Component.Name, Mean_CrossRef)
          colnames(Mean_CrossRef)<-c("Component.Name","Mean_CrossRef")
        }
        
        #if the user wants to upload their own file
        if (input$Coef_Upload == TRUE){
          MS_CrossRef<- read.csv(input$CrossRef_File$datapath,
                                 header = input$Header_CrossRef_File,
                                 sep = input$sep_CrossRef_File)
          
        }
        
        #if user wants to use the set standards 
        if (input$Coef_Upload == FALSE){
          #Merge both tables and get the cross-reference coefficient 
          if (input$MS_Machine == "MS2"){
            
            #read file with crossRef their own cross reference 
            MS_CrossRef<-read.csv(file = "coefficients/MS2_CrossRef.csv",
                                  header = TRUE,
                                  sep = ",")
          }
          
          if (input$MS_Machine == "MS3"){
            
            MS_CrossRef<-read.csv(file = "coefficients/MS3_CrossRef.csv",
                                  header = TRUE,
                                  sep = ",")
          }
          
          if (input$MS_Machine == "MS4"){
            MS_CrossRef<-read.csv(file = "coefficients/MS4_CrossRef.csv",
                                  header = TRUE,
                                  sep = ",")
            
          }
          
        }
        
        colnames(MS_CrossRef)<-c("Component.Name", "MS_CrossRef") #ensure columns names are standized
        
        #output table on the coefficients info page
        output$CrossRefCoef_Table<-renderDataTable(MS_CrossRef,  extensions = c('FixedColumns',"FixedHeader"),
                                                   options = list(paging = TRUE, pagelength = 20, scrollY = TRUE, 
                                                                  autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                                   selection = "single", rownames = TRUE)
        if(input$CrossRef_Incl == TRUE){
          Mean_CrossRef <- merge(Mean_CrossRef, MS_CrossRef, by= 'Component.Name', all.y = TRUE)
          
          Mean_CrossRef["CrossRef_Coef"]<-Mean_CrossRef["Mean_CrossRef"]/Mean_CrossRef["MS_CrossRef"]
          Cross_Ref_Coef<-dplyr::select(Mean_CrossRef, Component.Name, CrossRef_Coef)
        }
        
        if(input$CrossRef_Incl == FALSE){
          #make table with component name and crossRef value of 1
          Cross_Ref_Coef<-dplyr::select(MS_CrossRef, Component.Name)
          CrossRef_Coef<-rep(1,times=length(Cross_Ref_Coef$Component.Name))
          Cross_Ref_Coef<-cbind(Cross_Ref_Coef, CrossRef_Coef)
        }
        
        
        
        
        
        
        ###Step 2: Filter low area value to zero----
        
        #get samples of interest
        Samples <- MS_Table_Filtered
        
        #select are the names 
        #exclude standards, and cross references 
        for (x in input$not_Samples){
          Samples <- Samples[Samples["Sample.Name"] != x, ]
        }
        
        #obtain the Raw data in a table 
        Raw<- dplyr::select(Samples, Sample.Name, Component.Name, Area)
        Raw$Area <-as.numeric(Raw$Area) 
        Raw$Area <-round(Raw$Area, 3) #round to 3 digits
        Raw_Table<- pivot_wider(Raw, names_from = Component.Name, values_from = Area)
        Raw_Table<-tibble::column_to_rownames(Raw_Table, "Sample.Name")
        
        #reorder samples names to be in order samples were run 
        Raw_Table <- Raw_Table[Sample_Order, ]
        
        output$Raw_Table<-renderDataTable(Raw_Table, extensions = c('FixedColumns',"FixedHeader"),
                                          options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, 
                                                         autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                          selection = "single", rownames = TRUE)
        
        #download table 
        output$DownloadRaw_Table <- downloadHandler(filename = function(){"Raw_Table.csv"}, 
                                                    content = function(fname){
                                                      write.csv(Raw_Table, fname, row.names = TRUE)})
        
        
        
        
        #filter area to change area below threshold to zero 
        Samples$Area <- as.numeric(Samples$Area)
        Samples$Area[Samples$Area <input$Area_Filter] <- 0
        
        #Step 3: Apply the Cross-Reference Coefficient to Samples------------------
        
        #add the rename values 
        Rename<-read.csv(file = "coefficients/Rename_Comp.csv", header = TRUE, sep = ",")
        colnames(Rename)<-c("Comp_Name", "Rename")
        #Create a dictionary with values as the old names and the rename as the values
        Rename_Comp<-setNames(Rename$Rename, Rename$Comp_Name)
        
        
        ##put coefficents into file!! 
        
        Samples["Component.Rename"]<-Samples$Component.Name
        
        
        for (x in 1:length(Samples$Component.Name)){
          if (Samples$Component.Rename[x] %in% Rename$Comp_Name){
            Samples$Component.Rename[x] <-Rename_Comp[[Samples$Component.Rename[x]]]
          }else{
            Samples$Component.Rename[x]<-Samples$Component.Rename[x]
          }
        }
        
        
        #Now get the metabolites based on the renamed 
        
        #add the rename values 
        Met_IntStd<-read.csv(file = "coefficients/Met_IntStd.csv", header = TRUE, sep = ",")
        colnames(Met_IntStd)<-c("Metabolite", "Internal_Standard")
        
        #Create a dictionary with values as the old names and the rename as the values
        
        IntStd<-setNames(Met_IntStd$Internal_Standard, Met_IntStd$Metabolite)
        
        ##put coefficents into file!! 
        Samples["IS_Name"]<-Samples$IS.Name
        
        NoIS<-c("N/A", "(No IS)")
        for (x in 1:length(Samples$IS.Name)){
          if (Samples$IS_Name[x] %in% NoIS ){
            if (Samples$Component.Rename[x] %in% Met_IntStd$Metabolite){
              Samples$IS_Name[x] <-IntStd[[Samples$Component.Rename[x]]]
            }else{
              Samples$IS_Name[x]<-Samples$IS.Name[x]
            }
            
          }else{
            Samples$IS_Name[x]<-Samples$IS.Name[x]
          }
        }
        
        Samples<-within(Samples, rm(IS.Name))
        names(Samples)[names(Samples) == 'IS_Name'] <- 'IS.Name'
        
        
        #add cross-reference coefficent 
        Samples_Coef <- merge(Samples, Cross_Ref_Coef, by.x ="IS.Name", by.y="Component.Name", all.x = TRUE)
        
        #replace NAN with 1 so mediators without internal standard have area that isnt corrected
        Samples_Coef$CrossRef_Coef<- replace_na(Samples_Coef$CrossRef_Coef,1)
        
        #divide area by cross-reference coefficient 
        Samples_Coef$Area<-as.numeric(Samples_Coef$Area)#convert area values to numeric
        Samples_Coef["CrossRef_Area"]<-Samples_Coef["Area"]/Samples_Coef["CrossRef_Coef"] 
        
        #obtain Area after cross reference applied table 
        Area_CrossCoef<-dplyr::select(Samples_Coef, Sample.Name, Component.Name, CrossRef_Area)
        Area_CrossCoef$CrossRef_Area<-round(Area_CrossCoef$CrossRef_Area, 3) #round to 3 digits
        Area_CrossCoef_Table<- pivot_wider(Area_CrossCoef, names_from = Component.Name, values_from = CrossRef_Area)
        Area_CrossCoef_Table<-tibble::column_to_rownames( Area_CrossCoef_Table, "Sample.Name") #make row names 
        #reorder samples names to be in order samples were run 
        Area_CrossCoef_Table<-Area_CrossCoef_Table[Sample_Order, ]
        
        output$CrossRefTable<-renderDataTable(Area_CrossCoef_Table,  extensions = c('FixedColumns',"FixedHeader"),
                                              options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, 
                                                             autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                              selection = "single", rownames = TRUE)
        
        #download table 
        output$DownloadCrossRefTable <- downloadHandler(filename = function(){"CrossRef_Area_Table.csv"}, 
                                                        content = function(fname){
                                                          write.csv(Area_CrossCoef_Table, fname, row.names = TRUE)})
        
        
        
        
        ###Step 4: Convert area to pg------
        
        #if the user wants to upload their own conversion coefficient 
        if (input$Coef_Upload == TRUE){
          Conv_Coef_Table<- read.csv(input$ConvCoef_File$datapath,
                                     header = input$Header_ConvCoef_File,
                                     sep = input$sep_ConvCoef_File)
          
        } 
        if (input$Coef_Upload == FALSE){
          #get the conversion coef for MS machine
          if (input$MS_Machine == "MS2"){
            Conv_Coef_Table <-read.csv(file = "coefficients/MS2_ConvCoef.csv",
                                       header = TRUE,
                                       sep = ",")
          }
          
          if (input$MS_Machine == "MS3"){
            
            Conv_Coef_Table <-read.csv(file = "coefficients/MS3_ConvCoef.csv",
                                       header = TRUE,
                                       sep = ",")
          }
          
          if (input$MS_Machine == "MS4"){
            Conv_Coef_Table <-read.csv(file = "coefficients/MS4_ConvCoef.csv",
                                       header = TRUE,
                                       sep = ",")
          }
          
        }
        
        colnames(Conv_Coef_Table)<-c("Mediators", "Conv_Coef")
        
        #merge the conversion coefficient table with the mediators 
        Samples_Conv<-merge(Samples_Coef, Conv_Coef_Table, by.x = "Component.Rename", by.y = "Mediators", all.x = TRUE)
        
        #output table on the coefficients info page
        output$ConvCoef_Table<-renderDataTable(Conv_Coef_Table, extensions = c('FixedColumns',"FixedHeader"),
                                               options = list(paging = TRUE, pagelength = 20, scrollY = TRUE, 
                                                              autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                               selection = "single", rownames = TRUE)
        
        
        
        Samples_Conv["Amount_pg"]<-Samples_Conv["CrossRef_Area"]/Samples_Conv["Conv_Coef"]
        
        #obtain amount after aonversion coefficent applied applied table 
        Amount_convCoef<-dplyr::select(Samples_Conv, Sample.Name, Component.Name, Amount_pg)
        Amount_convCoef$Amount_pg<-round(Amount_convCoef$Amount_pg, 3) #round to 3 digits
        Amount_convCoef_Table<- pivot_wider(Amount_convCoef, names_from = Component.Name, values_from = Amount_pg)
        Amount_convCoef_Table<-tibble::column_to_rownames(Amount_convCoef_Table, "Sample.Name") #make row names 
        
        #reorder samples names to be in order samples were run 
        Amount_convCoef_Table<-Amount_convCoef_Table[Sample_Order, ]
        
        output$AmountTable<-renderDataTable(Amount_convCoef_Table, extensions = c('FixedColumns',"FixedHeader"),
                                            options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, 
                                                           autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                            selection = "single", rownames = TRUE)
        #download table 
        output$DownloadAmountTable <- downloadHandler(filename = function(){"Amount_Table.csv"}, 
                                                      content = function(fname){
                                                        write.csv(Amount_convCoef_Table, fname, row.names = TRUE)})
        
        
        
        ###Step 5: Normalization coeeficent to account for sample loss during prep ----
        
        
        #select only internal standards
        Sample_Amount<-Amount_convCoef_Table
        Sample_Amount_IS <-dplyr::select(Sample_Amount, MS_CrossRef$Component.Name)
        
        #transpose the table and get the max value foe each column
        Sample_Amount_T<-as.data.frame(t(Sample_Amount_IS))
        
        #select only sample with the max value as coefficient 
        Sample_Amount_max<-dplyr::select(Sample_Amount_T, input$MaxCol_Name)
        colnames(Sample_Amount_max)[1]<-"Conversion_Coef"
        
        #merge max value with the IS for each sample
        IS_Recovery<-merge(Sample_Amount_max, Sample_Amount_T, by = "row.names" )
        row.names(IS_Recovery)<-IS_Recovery$Row.names
        IS_Recovery<-IS_Recovery[,-1]
        
        ISSample_Names<-colnames(IS_Recovery)
        
        for (i in 1:length(IS_Recovery$Conversion_Coef)){
          
          for (j in 2:length(ISSample_Names)){
            IS_Recovery[i,j]<-IS_Recovery[i,j]/IS_Recovery[i,1]
          }
        }
        
        #get table with the percentage recovery of the standards, rounding to 2 digits 
        Per_Recovery<-round(t(IS_Recovery[,-1]*100), 2)
        
        output$PerRecovery<-renderDataTable(Per_Recovery, extensions = c('FixedColumns',"FixedHeader"),
                                            options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, 
                                                           autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                            selection = "single", rownames = TRUE)
        
        #download recovery table 
        output$DownloadPerRecovery <- downloadHandler(filename = function(){"Percent_Recovery_Table.csv"}, 
                                                      content = function(fname){
                                                        write.csv(Per_Recovery, fname, row.names = TRUE)})
        
        
        
        
        
        #do a loop to process each sample  
        
        Sample_pg <- dplyr::select(Samples_Conv, IS.Name, Component.Name, Sample.Name, Amount_pg)
        
        
        #get list of internal standards 
        IS_Mediators <- rownames(Sample_Amount_T)
        IS_Mediators<-gsub(".", " ", IS_Mediators, fixed=TRUE)
        
        #remove Internal standard mediators from sample_pg df
        for (x in IS_Mediators){
          Sample_pg<-Sample_pg[Sample_pg["Component.Name"] != x,]
        }
        
        
        #get list of sample names 
        ISSample_Names<-ISSample_Names[-1] #remove cofficient name 
        
        #add metabolite names as a column
        IS_Recovery_Names<-rownames(IS_Recovery)
        IS_Recovery_Names<-gsub(".", " ", IS_Recovery_Names, fixed=TRUE)
        IS_Recovery["IS_Names"]<-as.character(IS_Recovery_Names)
        
        #emplty table
        Final_pg <- data.frame(IS.Name ="del",
                               Component.Name = 'del',
                               Sample.Name = 'del',
                               Amount_pg = 'del', 
                               Norm_Coef = 'del')
        
        
        #for each sample name merge the normalization coefficient 
        for (y in ISSample_Names){
          Table <-as.data.frame(Sample_pg[Sample_pg["Sample.Name"] == y,]) #get a table with one sample
          IS_Coef <- as.data.frame(dplyr::select(IS_Recovery, y, IS_Names)) #get coefficient with that sample
          colnames(IS_Coef)<-c("Norm_Coef", "IS_Names") 
          Table1 <- merge(Table, IS_Coef, by.x = "IS.Name", by.y = "IS_Names", all.x = TRUE) #merge two tables
          Table["Norm_Coef"]<-Table1$Norm_Coef #get the norm coefiiceint for the sample
          Table2<-dplyr::select(Table1, IS.Name, Component.Name, Sample.Name,  Amount_pg, Norm_Coef)
          Final_pg <-rbind(Final_pg, Table2) #add to final table 
        }
        
        #remove first row
        Final_pg<-Final_pg[-1, ]
        
        #divide conv_area by normaliation coefficient
        Final_pg$Amount_pg<-as.numeric(Final_pg$Amount_pg)
        Final_pg$Norm_Coef<-as.numeric(Final_pg$Norm_Coef)
        
        Final_pg["Metabolite_pg"]<-Final_pg["Amount_pg"]/Final_pg["Norm_Coef"]
        
        
        #obtain amount after normalization coefficent applied applied table 
        Amount_Norm<-dplyr::select(Final_pg, Sample.Name, Component.Name, Norm_Coef)
        Amount_Norm$Norm_Coef<-round(Amount_Norm$Norm_Coef, 3) #round to 3 digits
        Amount_Norm_Table<- pivot_wider(Amount_Norm, names_from = Component.Name, values_from = Norm_Coef)
        Amount_Norm_Table<-tibble::column_to_rownames(Amount_Norm_Table, "Sample.Name") #make row names 
        
        #reorder samples names to be in order samples were run 
        Amount_Norm_Table<-Amount_Norm_Table[Sample_Order, ]
        
        #output of datatable that user can scroll through and have fixed columns/header 
        output$NormTable<-renderDataTable(Amount_Norm_Table, extensions = c('FixedColumns',"FixedHeader"),
                                          options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, 
                                                         autoWidth = TRUE, fixedColumns = row.names, fixedHeader = TRUE), 
                                          selection = "single", rownames = TRUE)
        
        #download table 
        output$DownloadNormTable <- downloadHandler(filename = function(){"Normalized_Table.csv"}, 
                                                    content = function(fname){
                                                      write.csv(IS_Recovery, fname, row.names = TRUE)})
        
        ##Step 5: convert amount (pg) to concentration (pg/ul)------
        
        #if user selects volume
        if (input$SameVol_Data == TRUE){
          Final_pg["Concentration"] = (Final_pg["Metabolite_pg"]/(input$SameData_Amount))*(input$SameData_Standard)
        }
        
        #if select weight then work with uploaded file 
        
        
        if (input$DiffVol_Data == TRUE){
          
          #read sample mass file
          Sample_Mass<- read.csv(input$DiffVol_File$datapath,
                                 header = input$Header_DiffVol,
                                 sep = input$sep_DiffVol)
          
          #merge with the final time 
          Final_pg<-merge(Final_pg, Sample_Mass, by.x = "Sample.Name", by.y = "Sample_Name", all.x = TRUE)
          
          ##get the mean mass and round it to remove digits --- DO NOT CALC MEAN MASS
          #Mean_Mass<-round(mean(Final_pg$Weight), digits = 0)
          
          #get final amount 
          Final_pg["Concentration"] = (Final_pg["Metabolite_pg"]/Final_pg["Weight"])*(input$DiffData_Standard)
          
        }
        
        
        
        
        ###Step 6: Structuring output 
        
        #obtain columns of interest
        Final <- dplyr::select(Final_pg, Component.Name, Sample.Name, Concentration)
        
        #final pivot table
        Final_1<- dplyr::select(Final, Sample.Name, Component.Name, Concentration)
        Final_1$Concentration<-round(Final_1$Concentration, 3) #round to 3 digits 
        Final_1_Pivot<-pivot_wider(Final_1, names_from = Component.Name, values_from = Concentration) #pivot table to wider 
        Final_1_Pivot<-tibble::column_to_rownames(Final_1_Pivot, "Sample.Name") #make row names 
        
        
        
        Final_1_Pivot<-Final_1_Pivot[Sample_Order, ]
        #remove columns that are all NA
        Final_1_Pivot<- Final_1_Pivot[, colSums(is.na( Final_1_Pivot)) != nrow(Final_1_Pivot)]
        
        output$ConcTable<-renderDataTable(Final_1_Pivot, extensions = c('FixedColumns',"FixedHeader"),
                                          options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, fixedColumns = row.names,  fixedHeader = TRUE), 
                                          selection = "single", rownames = TRUE)
        
        output$DownloadConcTable <- downloadHandler(filename = function(){"ConcentrationTable.csv"}, 
                                                    content = function(fname){
                                                      write.csv(Final_1_Pivot, fname, row.names = TRUE)})
        #make error message outputs blank
        output$ConcTable_Error<- output$Raw_Table_Error<- output$CrossRefTable_Error<-output$AmountTable_Error<-
          output$PerRecovery_Error<-output$NormTable_Error<-renderUI(return())
      
      },
      #what happens when there is an error in the code 
      error = function(e) {
        #notification box with specific error
        showNotification(paste0(e), type = 'err')
        #error messages for each tab explaining what could be wrong
        output$ConcTable_Error<- output$Raw_Table_Error<- output$CrossRefTable_Error<-output$AmountTable_Error<-
          output$PerRecovery_Error<-output$NormTable_Error<-
          renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                               " - File Formats and seperator is correct and All inputs filled.", 
                               " - File contains the follwoing columns: Sample.Name,  Sample.Index, Component.Type, Component.Name,  IS.Name, Area",
                               " - If using set coefficients, confirm correct mass spectometry machine is selected",
                               " - If using cross reference Coefficient, ensure that there are 3 cross reference samples with the sample name",
                               sep = "<br/>"))})
        #make table outputs empty
        output$ConcTable<- output$Raw_Table<- output$CrossRefTable<-output$AmountTable<-
          output$PerRecovery<-output$NormTable<-renderDataTable(return())
      }, silent=FALSE)
    
    
   
   
  })
  
  
  ####Code to autopopulate mediator names once Data Format file is uploaded ----
  observe({
    req(input$Data_File2)
    
    options(digits = 3)
    
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {Conc_Table<-read.csv(input$Data_File2$datapath,
                           header = input$Header_Data2,
                           sep = input$sep_Data2,
                           row.names = 1)
      #get list of metabolites 
      Met_Names<-sort(unique(as.character(unname(unlist(colnames(Conc_Table))))))
      if (Met_Names == ''){
        updateSelectInput(session,"Met_DataList2", choices = "Error!!Please check file and Separator")
      }
      updateSelectInput(session,"Met_DataList2", choices = Met_Names)
      
    }, 
    
    #if there is an error in the above code, the error message will display 
    error = function(e) {
      showNotification('Error! Please check inputs are correct','',type = "error")
      updateSelectInput(session,"Met_DataList2", choices = "Error!!Please check file and Separator")
      
    }, silent= FALSE)
    
  })
  
  ###Code to process final concentration table to obtain the statistics and Machine Learning inpu tables----
  observeEvent(input$Go_Data_Format, {
    tryCatch(
      
      ##try to run the code to process the table 
      {Conc_Table<-read.csv(input$Data_File2$datapath,
                           header = input$Header_Data2,
                           sep = input$sep_Data2,
                           row.names = 1)
      
      if (input$Filtering_Data == TRUE){
        Met_List<-as.character(input$Met_DataList2)
        Conc_Table<-dplyr::select(Conc_Table, Met_List)}

      Stats_Final<-Conc_Table
      
      #remove the transition codes
      colnames(Stats_Final)<-gsub("\\.+[0-9]+[0-9]+[0-9]$","", colnames(Stats_Final))


      Met_Rename<-c()

      #remove spaces from column names
      colnames(Stats_Final)<-gsub("\\.","", colnames(Stats_Final))
      #loop to put correct final name
      for (x in 1:ncol(Stats_Final)){
        if (colnames(Stats_Final)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(Stats_Final)[x]]]
        }else{
          Met_Rename[x] <-colnames(Stats_Final)[x]
        }
      }
      #change column names to be in correct format
      colnames(Stats_Final)<-Met_Rename

      #outputs and download for stats table
      output$Stats_Table<-renderDataTable(Stats_Final, extensions = c('FixedColumns',"FixedHeader"),
                                          options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, fixedColumns = row.names,  fixedHeader = TRUE),
                                          selection = "single", rownames = TRUE)

      output$DownloadStats_Table <- downloadHandler(filename = function(){"Stats_Table.csv"},
                                                    content = function(fname){
                                                      write.csv(Stats_Final, fname, row.names = TRUE)})



      #Loop to get the fatty acids for the metabolite and add to first row ----
      Fatty_Acid<-c()

      #fatty acid points
      for (i in 1:ncol(Stats_Final)){
        Fatty_Acid[i]<-SPMs_FA[colnames(Stats_Final)[i]]
      }


      MLTable_Final<-rbind(Fatty_Acid, Stats_Final)
      row.names(MLTable_Final)[1]<-"FattyAcid"

      output$ML_Input_Table<-renderDataTable(MLTable_Final, extensions = c('FixedColumns',"FixedHeader"),
                                             options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE, fixedColumns = row.names,  fixedHeader = TRUE),
                                             selection = "single", rownames = TRUE)

      output$DownloadML_Input_Table <- downloadHandler(filename = function(){"ML_Table.csv"},
                                                       content = function(fname){
                                                         write.csv(MLTable_Final, fname, row.names = TRUE)})
      output$StatsInputTable_Error<-renderText(return())
      output$MLInputTable_Error<-renderText(return())
      
    },
    error = function(e) {
      showNotification(paste0(e), type = 'err')
      #showNotification('There was an error: Check if File and seperator are the correct format','',type = "error")
      output$StatsInputTable_Error<-renderText("Error: Check File Formats and seperator is correct.")
      output$MLInputTable_Error<-renderText("Error: Check File Formats and seperator is correct.")
      output$Stats_Table<-renderDataTable(return())
      output$ML_Input_Table<-renderDataTable(return())
    }, silent=FALSE)
    
    
  })
  
  
  observe({
    req(input$PCAPLSDA_file)
    
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {
        #form dataframe with file 
        df<- read.csv(input$PCAPLSDA_file$datapath,
                      header = input$header_PCA,
                      sep = input$sep_PCA)
        #transpose table if samples are at columns
        if(input$samp_PCA == "col_PCA"){
          df<-data.frame(t(df))
          df_samp<-rownames(df)
          df<-data.frame(cbind(df_samp, df))
          colnames(df)<- df[1,]
          df<- df[-1, ]}
        
        #get list of column names as options 
        Group_Names<-unique(as.character(unname(unlist(colnames(df)))))
        
        #if statement if the list of column names is empty
        if (Group_Names == ''){
          updateSelectInput(session,"Group_PCA", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"SampleCol_PCA", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
        }
        
        #update inputs if correct 
        updateSelectInput(session, "Group_PCA", choices = Group_Names, selected = Group_Names[[1]])
        updateSelectInput(session, "SampleCol_PCA", choices = Group_Names, selected = Group_Names[[1]])
        
      }, 
      
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Group_PCA", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"SampleCol_PCA", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
  })

  
  
  observeEvent(input$Go_PCAPLSDA|input$Go_PCAPLSDA2, {
    req(input$PCAPLSDA_file)
    
    tryCatch({
      #form dataframe with file 
      df<- read.csv(input$PCAPLSDA_file$datapath,
                    header = input$header_PCA,
                    sep = input$sep_PCA)
      
      #transpose table if samples are at columns
      if(input$samp_PCA == "col_PCA"){
        df<-as.data.frame(t(df))
        df_samp<-rownames(df)
        df<-data.frame(cbind(df_samp, df))
        colnames(df)<- df[1,]
        df<- df[-1, ]
        
        #remove spaces from column names 
        colnames(df)<-gsub("[^a-zA-Z0-9]","", colnames(df)) 
        
        for (x in 1:ncol(df)){
          if (substr(colnames(df)[x],1,1) %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
            colnames(df)[x]= paste("X", colnames(df)[x], sep = '')
          }
          else{colnames(df)[x] = colnames(df)[x]}
        }
      }
      
      
      #make sure metabilite names are standardized
      Met_Rename<-c()
      
      #remove spaces from column names 
      colnames(df)<-gsub("[^a-zA-Z0-9]","", colnames(df)) 
      
      
      for (x in 1:ncol(df)){
        if (colnames(df)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(df)[x]]]
        }else{
          Met_Rename[x] <-colnames(df)[x]
        }
      }
      colnames(df)<-Met_Rename
      
      #script for PCA analysis
      if(input$test_PCAPLSDA == "PCA"){
        
        GroupName<- as.character(input$Group_PCA)
        SampleName<-as.character(input$SampleCol_PCA)
        
        #save group columns
        Groups_PCA<-df[,GroupName]
        Sample_PCA<-df[,SampleName]
        
        #remove group column from dataframe
        #df_PCA<-df[,-c(input$SampCol_PCA)]
        df_PCA<-df[ , !(names(df) %in% c(GroupName, SampleName))]
        
        #ensure dataframe is numeric
        df_PCA<-data.frame(sapply(df_PCA, function(x) as.numeric(as.character(x))))
        
        #zero handling data
        cols<- 1:(ncol(df_PCA))
        #replace zeros for each column to 1/5 the smallest value for each column
        df_PCA[cols] <- lapply(df_PCA[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
        
        #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
        df_PCA<-as.matrix(df_PCA)
        df_PCA[is.infinite(df_PCA)] <- NA
        min_index <- which.min(df_PCA)
        zero_replace <- (df_PCA[min_index]/5)
        df_PCA <- as.data.frame(df_PCA)
        df_PCA[is.na(df_PCA)] <- zero_replace
        
        #exclude columns with the same value for all rows
        df_PCA<-Filter(var, df_PCA)
        
        
        #PCA test with normalization 
        PCA<-prcomp(df_PCA, scale = TRUE)
        
        #----------PCA Score (PC1vsPC2)
        
        #get PCA Scores into dataframe and change it to numeric
        PCA_Scores<-data.frame(PCA$x)
        PCA_Scores<-data.frame(sapply(PCA_Scores, function(x) as.numeric(as.character(x))))
        PCA_Scores<-(cbind(Sample_PCA, Groups_PCA, PCA_Scores))
        
        #PC1 vs PC2 scores plot
        if(input$ViewSamples_PCA == TRUE){
          output$Xscoresplot<-renderPlotly({ggplotly(ggplot(PCA_Scores, aes(x=PC1, y=PC2, color = Groups_PCA)) + geom_point() + 
                                                       stat_ellipse(level = as.numeric(input$PCA_Confidence), geom = "polygon", aes(fill = Groups_PCA), alpha = 0.25) + theme_bw() + 
                                                       geom_text(label=Sample_PCA, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T))})
        }
        else{
          output$Xscoresplot<-renderPlotly({ggplotly(ggplot(PCA_Scores, aes(x=PC1, y=PC2, color = Groups_PCA)) + geom_point() + 
                                                       stat_ellipse(level = as.numeric(input$PCA_Confidence), geom = "polygon", aes(fill = Groups_PCA), alpha = 0.25) + theme_bw())})
        }
        
        #download PC scores 
        output$DownloadScore <- downloadHandler(filename = function(){"PCA_Score.csv"}, 
                                                content = function(fname){
                                                  write.csv(PCA_Scores, fname, row.names = FALSE)})
        
        #----------PCA Loadings 
        
        #get loading scores, convert to numeric, get row names as column 
        loading<-data.frame(PCA$rotation)
        Variable<-rownames(loading)
        
        loading<-data.frame(sapply(loading, function(x) as.numeric(as.character(x))))
        loading<-cbind(Variable, loading)
        
        #change the names to be in the proper format 
        loading$Variable[grep("^X", loading$Variable)] <- paste(gsub("^X", "'", loading$Variable[grep("^X", loading$Variable)]), "'", sep = "")
        loading$Variable <- gsub("4$", "[4]", loading$Variable) #converts names ending with 4 to subscript 5
        loading$Variable <- gsub("B4'", "B'[4]", loading$Variable) #B4 converted to subscript 4
        loading$Variable <- gsub("\\.", "-", loading$Variable) # Replace "." with "-" when required 
        loading$Variable <- gsub("n.3.DPA", "[n-3~DPA]", loading$Variable) # Transform n-3 DPA as a subscript: 
        lipids <- as.character(loading$Variable)
        
        loading<-loading[order(-loading$PC1),] #order dataframe based on PC1 score
        loading$Variable<-factor(loading$Variable, levels = loading$Variable)
        
        Loading_plot<-ggplot(loading, aes(x=reorder(Variable, PC1), y=PC1)) + geom_bar(stat = "identity") + 
          scale_x_discrete(labels = parse(text = lipids)) +
          coord_flip() + labs(x = "Lipid Mediator", y= "Loading Value") +
          theme(axis.title = element_text(size = 20),
                axis.text.x  = element_text(size = 15, hjust = 0.5),
                axis.text.y  = element_text(size = 15, hjust = 1),
                panel.background = element_rect(fill = "white")) 
        
        
        #plot the loading plot
        output$Loadingplot<-renderPlot(Loading_plot)
        
        
        #download loading plot
        output$downloadLoading_Plot <- downloadHandler(filename = function(){paste("Loading_Plot",'.png',sep='')},
                                                       content = function(file){
                                                         ggsave(file,plot=Loading_plot, width = 12, height = 12)})
        
        #download Loading Data
        output$DownloadLoading <- downloadHandler(filename = function(){"PCA_LoadingScore.csv"}, 
                                                  content = function(fname){
                                                    write.csv(loading, fname, row.names = FALSE)})
        
        #----------Percentage of Variance 
        
        #percentage of variance
        Importance<-summary(PCA) #summary of PCA without normalization
        Importance<-t(Importance$importance)#save only the importance matrix
        Component<-rownames(Importance)
        Importance<-data.frame(cbind(Component, Importance[,2]))
        colnames(Importance)[2]<-"Percentage_of_Variance"
        Importance$Percentage_of_Variance<-(as.numeric(Importance$Percentage_of_Variance)*100)
        Importance$Component<-factor(Importance$Component, levels = Importance$Component)
        
        #plot the percentage variance 
        output$PerVarplot<- renderPlotly(ggplotly(ggplot(Importance, aes(Component, Percentage_of_Variance)) + 
                                                    geom_bar(stat = "identity", fill= "skyblue") + labs(y="Percentage of Variance") +
                                                    theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        
        #download Percent Variance Data
        output$DownloadVar <- downloadHandler(filename = function(){"PCA_PercentVariance.csv"}, 
                                              content = function(fname){
                                                write.csv(Importance, fname, row.names = FALSE)})
        
        #ensures that the VIP plot is 
        output$VIPplot<-NULL
        output$VIPtable<-NULL
        output$DownloadVIP <-NULL
        output$VIPtext <- renderText({ 
          "Please select PLS-DA Analysis to obtain VIP Scores "
        })
        
      }
      #script for PLS-DA analysis
      else{
        
        GroupName<- as.character(input$Group_PCA)
        SampleName<-as.character(input$SampleCol_PCA)
        
        #save group columns
        Groups_PLSDA = df[,GroupName]
        Sample_PLSDA = df[,SampleName]
        
        #remove group column from dataframe
        df_PLSDA<-df[ , !(names(df) %in% c(GroupName, SampleName))]
        Groups_PLSDA1<-as.factor(Groups_PLSDA)
        
        #ensure dataframe is numeric
        df_PLSDA<-data.frame(sapply(df_PLSDA, function(x) as.numeric(as.character(x))))
        
        #zero handling data
        cols<- 1:(ncol(df_PLSDA))
        #replace zeros for each column to 1/5 the smallest value for each column
        df_PLSDA[cols] <- lapply(df_PLSDA[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
        
        #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
        df_PLSDA<-as.matrix(df_PLSDA)
        df_PLSDA[is.infinite(df_PLSDA)] <- NA
        min_index <- which.min(df_PLSDA)
        zero_replace <- (df_PLSDA[min_index]/5)
        df_PLSDA <- as.data.frame(df_PLSDA)
        df_PLSDA[is.na(df_PLSDA)] <- zero_replace
        
        #exclude columns with the same value for all rows
        df_PLSDA<-Filter(var, df_PLSDA)
        
        #get the number of components using PCA test to get the length of the scores table = number of components 
        PCA<-prcomp(df_PLSDA, scale = TRUE)
        PCA_Scores<-data.frame(PCA$x)
        
        #run plsda() with data values and class
        PLSDA<- mixOmics::plsda(df_PLSDA, Groups_PLSDA, scale = TRUE, ncomp = (ncol(PCA_Scores)-1))
        
        
        #-------------PLS-DA X-scores
        
        #getting the scores and convert to numeric
        XScores<-data.frame(PLSDA$variates$X)
        XScores<-data.frame(sapply(XScores, function(x) as.numeric(as.character(x))))
        XScores<- cbind(Sample_PLSDA, Groups_PLSDA, XScores)
        
        #ploting comp1 vs comp 2 to show sample names 
        if(input$ViewSamples_PCA == TRUE){
          output$Xscoresplot<-renderPlotly({ggplotly(ggplot(XScores, aes(x=comp1, y=comp2, color = Groups_PLSDA)) + geom_point() + 
                                                       stat_ellipse(level = as.numeric(input$PCA_Confidence), geom = "polygon", aes(fill = Groups_PLSDA), alpha = 0.25) + 
                                                       theme_bw() + labs(x = "Component 1", y= "Component 2") + 
                                                       geom_text(label=Sample_PLSDA, nudge_x = 0.25, nudge_y = 0.25, check_overlap = T))})
        }
        #plotting comp1 vs comp 2 wihtout sample names 
        else{
          output$Xscoresplot<-renderPlotly({ggplotly(ggplot(XScores, aes(x=comp1, y=comp2, color = Groups_PLSDA)) + geom_point() + 
                                                       stat_ellipse(level = as.numeric(input$PCA_Confidence), geom = "polygon", aes(fill = Groups_PLSDA), alpha = 0.25) + 
                                                       labs(x = "Component 1", y= "Component 2") + theme_bw())})
        }
        
        #download PC scores 
        output$DownloadScore <- downloadHandler(filename = function(){"PLSDA_Score.csv"}, 
                                                content = function(fname){
                                                  write.csv(XScores, fname, row.names = FALSE)})
        
        ##get the p-value and the Q2 R2 plot 
        
        #get the R2 and Q2 values for the first component 
        PLSDA_R2Q2<-ropls::opls(df_PLSDA, Groups_PLSDA, predI = 2, fig.pdfC = "none", info.txtC = "none")
        Summ_R2Q2<-ropls::getSummaryDF(PLSDA_R2Q2)
        
        #output datatable
        R2Q2_Table<- data.frame(R2Y = Summ_R2Q2$`R2Y(cum)`, Q2Y = Summ_R2Q2$`Q2(cum)`, P.Value_R2Y = Summ_R2Q2$pR2Y , P.Value_Q2Y = Summ_R2Q2$pQ2)
        output$R2Q2_Table<-renderDataTable(R2Q2_Table, rownames = FALSE)
        output$R2Q2_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("R2, Q2, and P Values for PLSDA", align = "center") })
        
        #plot for R2Q2 values 
        R2Q2<-as.data.frame(t(R2Q2_Table[,1:2])) #get the R2 and Q2 values
        R2Q2<-cbind(rownames(R2Q2), data.frame(R2Q2, row.names=NULL)) #change rownames to the first column
        colnames(R2Q2)<-c("PLSDA_R2Q2", "Values") #change the column names 
        
        output$R2Q2_Plot<- renderPlotly(ggplotly(ggplot(R2Q2, aes(PLSDA_R2Q2, Values)) + 
                                                   geom_bar(stat = "identity", fill = c("gray", "black")) +
                                                   theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        
        #-------------PLS-DA Loading Scores
        
        #obtan loading dataframes and convert to numeric 
        Xloading<-data.frame(PLSDA$loadings$X) #get x loadings as dataframe 
        Variable<-rownames(Xloading)
        Xloading<-as.data.frame(sapply(Xloading, function (x) as.numeric(as.character(x))))
        Xloading<-cbind(Variable, Xloading)
        
        #change the names to be in the proper format 
        Xloading$Variable[grep("^X", Xloading$Variable)] <- paste(gsub("^X", "'", Xloading$Variable[grep("^X", Xloading$Variable)]), "'", sep = "")
        Xloading$Variable <- gsub("4$", "[4]", Xloading$Variable) #converts names ending with 4 to subscript 5
        Xloading$Variable <- gsub("B4'", "B'[4]", Xloading$Variable) #B4 converted to subscript 4
        Xloading$Variable <- gsub("\\.", "-", Xloading$Variable) # Replace "." with "-" when required 
        Xloading$Variable <- gsub("n.3.DPA", "[n-3~DPA]", Xloading$Variable) # Transform n-3 DPA as a subscript: 
        lipids <- as.character(Xloading$Variable)
        
        Xloading<-Xloading[order(-Xloading$comp1),] #order dataframe based on PC1 score
        Xloading$Variable<-factor(Xloading$Variable, levels = Xloading$Variable)
        
        XLoading_plot<-ggplot(Xloading, aes(x=reorder(Variable, comp1), y=comp1)) + geom_bar(stat = "identity") + 
          scale_x_discrete(labels = parse(text = lipids)) +
          coord_flip() + labs(x = "Lipid Mediator", y= "Loading Value") +
          theme(axis.title = element_text(size = 20),
                axis.text.x  = element_text(size = 15, hjust = 0.5),
                axis.text.y  = element_text(size = 15, hjust = 1),
                panel.background = element_rect(fill = "white")) 
        
        
        #plot the loading plot
        output$Loadingplot<-renderPlot(XLoading_plot)
        
        
        #Download loading plot
        output$downloadLoading_Plot <- downloadHandler(filename = function(){paste("Loading_Plot",'.png',sep='')},
                                                       content = function(file){
                                                         ggsave(file,plot=XLoading_plot, width = 12, height = 12)})
        #download Loading Data
        output$DownloadLoading <- downloadHandler(filename = function(){"PLSDA_LoadingScore.csv"}, 
                                                  content = function(fname){
                                                    write.csv(Xloading, fname, row.names = FALSE)})
        
        #-------------PLS-DA Percentage Variance
        
        #Variance dataframes 
        XVariance<-data.frame(PLSDA$prop_expl_var$X)
        colnames(XVariance)[1]<-"Percentage_of_Var"
        Component<-rownames(XVariance)
        XVariance<-cbind(Component, XVariance)
        XVariance$Percentage_of_Var<-as.numeric(XVariance$Percentage_of_Var)*100
        XVariance$Component<- factor(XVariance$Component, levels = XVariance$Component)
        
        #percentage variance plots
        output$PerVarplot<- renderPlotly(ggplotly(ggplot(XVariance, aes(Component, Percentage_of_Var)) + 
                                                    geom_bar(stat = "identity", fill = "skyblue") + labs(y = "Percentage of Variance") +
                                                    theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        #download Percent Variance Data
        output$DownloadVar <- downloadHandler(filename = function(){"PLSDA_PercentVariance.csv"}, 
                                              content = function(fname){
                                                write.csv(XVariance, fname, row.names = FALSE)})
        
        #-------------PLS-DA VIP Scores 
        
        #VIP scores plot
        VIP<-data.frame(vip(PLSDA)) #get VIP list
        Variable<-rownames(VIP) #get row names 
        VIP<-data.frame(sapply(VIP, function(x) as.numeric(as.character(x))))
        VIP<-data.frame(cbind(Variable, VIP)) #add variables as column
        
        
        Fatty_Acid <- 1:length(VIP$Variable)
        cbind(VIP, Fatty_Acid)
        
        #fatty acid points 
        for (i in 1:length(VIP$Variable)){
          VIP$Fatty_Acid[i]<-SPMs_FA[VIP$Variable[i]]
        }
        
        VIP$Fatty_Acid<-as.character(VIP$Fatty_Acid)
        
        VIP_Table <- VIP
        VIP_Table$Variable<- gsub("\\.", "-", VIP_Table$Variable) # Replace "." with "-" when required
        VIP_Table$Variable[grep("^X", VIP_Table$Variable)] <- paste(gsub("^X", "", VIP_Table$Variable[grep("^X", VIP_Table$Variable)]), "", sep = "") #remove x at the beginning
        
        VIP<-VIP[VIP["comp1"]>1,]
        
        #Change the naming format to be for a scientific paper
        VIP$Variable[grep("^X", VIP$Variable)] <- paste(gsub("^X", "'", VIP$Variable[grep("^X", VIP$Variable)]), "'", sep = "")
        VIP$Variable <- gsub("4$", "[4]", VIP$Variable) #everything ending with 4 will be subscript 4
        VIP$Variable <- gsub("B4'", "B'[4]", VIP$Variable) #4 will be changed with subscript
        VIP$Variable <- gsub("\\.", "-", VIP$Variable) # Replace "." with "-" when required 
        VIP$Variable <- gsub("n.3.DPA", "[n-3~DPA]", VIP$Variable) # Transform n-3 DPA as a subscript
        
        VIP<-VIP[order(VIP$comp1),] #order dataframe based on -comp1 score
        VIP$Variable <- factor(VIP$Variable, levels = VIP$Variable[order(VIP$comp1)]) #keeps variable column as factor
        
        
        lipids <- as.character(VIP$Variable)
        
        #add colors to SPM names
        lm_classes <- as.factor(VIP$Fatty_Acid)
        names(lm_classes) <- VIP$Variable
        
        dha_index<-which((lm_classes=="DHA")== TRUE)
        n_three_index<-which((lm_classes=="n3DPA")== TRUE)
        epa_index<-which((lm_classes=="EPA")== TRUE)
        aa_index<-which((lm_classes=="AA")== TRUE)
        
        lm_colors <- NULL
        lm_colors[dha_index] <- "blue"
        lm_colors[n_three_index] <- "brown"
        lm_colors[epa_index] <- "darkgoldenrod1"
        lm_colors[aa_index] <- "darkslategray"
        
        
        
        VIP_plot <- ggplot(data = VIP, aes(x = Variable, y = comp1, color = Fatty_Acid)) + geom_point(size =3) +
          scale_y_continuous(name = "VIP Score") +
          labs(x = "Lipid Mediators") +
          scale_x_discrete(labels = parse(text = lipids)) + 
          scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue", "n3DPA" = "brown", 
                                                                               "EPA"="darkgoldenrod1",  "AA" ="darkslategray")) +
          coord_flip() + 
          theme(axis.title = element_text(size = 25),
                axis.text.x  = element_text(size = 15, hjust = 0.5),
                axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                legend.position = "top",
                aspect.ratio = 2/1,
                legend.title = element_text(size = 20),
                legend.text  = element_text(size = 15),
                panel.background = element_rect(fill = "white"),
                panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
        
        
        
        #outputs 
        output$VIPtext<-NULL
        output$VIPplot<-renderPlot(VIP_plot)
        
        output$VIPtable<-renderDataTable(VIP, options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE), 
                                         selection = "single", rownames = FALSE)
        
        output$downloadPLSDA_VIP_Plot <- downloadHandler(filename = function(){paste("PLSDA_VIP_Plot",'.png',sep='')},
                                                         content = function(file){
                                                           ggsave(file,plot=VIP_plot, width = 12, height = 8)})
        output$DownloadVIP <- downloadHandler(filename = function(){"PLSDA_VIPScore.csv"}, 
                                              content = function(fname){
                                                write.csv(VIP_Table, fname, row.names = FALSE)})
        
        #Pathway Analysis PLSDA VIP score
        VIP_Table_Comp1<-VIP_Table[,1:2]
        #open csv with the pathways for each metabolite
        Pathway <-read.csv(file = "coefficients/Pathways2.csv", header = TRUE, sep = ",")
        colnames(Pathway)<- c("Metabolite", "Pathway")
        
        #set it up like a dictionary 
        MetPath<-setNames(Pathway$Pathway, Pathway$Metabolite)
        VIP_Table_Comp1["Pathways"]<-"N/A"#column for pathways with N/A as default
        
        #loop to put the pathway for the metabolite 
        for (x in 1:length(VIP_Table_Comp1$Variable)){
          if (VIP_Table_Comp1$Variable[x] %in% Pathway$Metabolite){
            VIP_Table_Comp1$Pathways[x] <- MetPath[[VIP_Table_Comp1$Variable[x]]]
          }else{
            VIP_Table_Comp1$Pathways[x]<-VIP_Table_Comp1$Pathways[x]
          }
        }
        
        #Keep rows that have a VIP Score greater than 1 and pathway != N/A
        VIP_Table_Comp1<-VIP_Table_Comp1[VIP_Table_Comp1["comp1"]>1,]
        VIP_Table_Comp1<-VIP_Table_Comp1[VIP_Table_Comp1["Pathways"] != "N/A",]
        
        #get list of unique pathways
        VIP_PathList<-unique(VIP_Table_Comp1$Pathways)
        #empty table
        VIPPath_Table<-data.frame(Pathway = "del", Number_SPM = 1)
        #loop to calculate the number of SPMs for each pathway
        for (x in VIP_PathList){
          x_VIP <- VIP_Table_Comp1[VIP_Table_Comp1["Pathways"] == x, ]
          Table<-data.frame(Pathway = x, Number_SPM = length(x_VIP$Pathways))
          VIPPath_Table<-rbind(VIPPath_Table, Table)
        }
        #remove plot 
        VIPPath_Table<-VIPPath_Table[-1,]
        VIPPath_Table<-VIPPath_Table[order(-VIPPath_Table$Number_SPM),] #reorder based on Number SPM
        
        #get columns in correct data format
        VIPPath_Table$Number_SPM<-as.numeric(VIPPath_Table$Number_SPM)
        VIPPath_Table$Pathway<-as.character(VIPPath_Table$Pathway)
        #save as factor to keep the order
        VIPPath_Table$Pathway<-factor(VIPPath_Table$Pathway, levels = VIPPath_Table$Pathway)
        
        #plot title and plot output
        output$VIP_Pathway_Title <- renderUI({req(input$Go_PCAPLSDA|input$Go_PCAPLSDA2); h2("VIP Pathways Plot", align = "center") })
        output$VIP_PathwayPlot<- renderPlotly(ggplotly(ggplot(VIPPath_Table, aes(x=reorder(Pathway, -Number_SPM), Number_SPM)) + 
                                                         geom_bar(stat = "identity", fill = c("#00C19A")) +
                                                         labs(x = "Pathway", y = "Number of SPMs")+
                                                         theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        output$Xscoresplot_Error<- output$Loadingplot_Error<- output$PerVarplot_Error<-output$VIPplot_Error<-
          renderUI(return())
        
      }
    },
    #what happens when there is an error in the code 
    error = function(e) {
      #notification box with specific error
      showNotification(paste0(e), type = 'err')
      #error messages for each tab explaining what could be wrong
      output$Xscoresplot_Error<- output$Loadingplot_Error<- output$PerVarplot_Error<-output$VIPplot_Error<-
        renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                             " - File Formats and seperator is correct and All inputs filled.", 
                             " - The correct group information and sample information column is selected",
                             " - The correct sample location is selected",
                             sep = "<br/>"))})
      #make table/plots outputs empty
      output$Xscoresplot<-output$R2Q2_Plot<-output$PerVarplot<-output$VIP_PathwayPlot<-renderPlotly(return())
      output$Loadingplot<-output$VIPplot<-plotOutput(return())
      output$R2Q2_Table<-renderDataTable(return())
      output$R2Q2_Title<-output$VIPtext<-output$VIP_Pathway_Title<-renderUI(return())
    }, silent=FALSE)
               
             
    
  })
  
  observe({
    inFile <- input$lm_profiles
    
    print(inFile)
    
    if(is.null(inFile))
      
      return(NULL)
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {
        # Open the file if samples are in column
        if(input$samp_Diff == "col_Diff"){
          lm_profiles<- read.csv(inFile$datapath,
                                 header = input$header_Diff,
                                 sep = input$sep_Diff,
                                 row.names = 1)
        }
        
        #transpose table if samples are in the rows 
        if(input$samp_Diff == "row_Diff"){
          lm_profiles = read.csv(inFile$datapath,
                                 header=input$header_Diff,
                                 sep=input$sep_Diff,
                                 row.names = 1)
          lm_profiles<-data.frame(t(lm_profiles))}
        
        # This option allows to chose the groups to make the comparison:
        
        options <- unique(as.character(unname(unlist(lm_profiles[1, ]))))
        
        #if statement if the list of column names is empty
        if (options == ''){
          updateSelectInput(session,"group_A", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"group_B", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
        }
        
        updateSelectInput(session, "group_A", choices = options, selected = options[[2]])
        updateSelectInput(session, "group_B", choices = options)
      }, 
      
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"group_A", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"group_B", choices = "Error!!Please chekc file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
  })
  
  observeEvent(input$Go_Diff|input$Go_DiffPlot, {
    
    inFile <- input$lm_profiles
    
    print(inFile)
    
    if(is.null(inFile))
      
      return(NULL)
    
    tryCatch(
      {
        
        # Open the file if samples are in column
        if(input$samp_Diff == "col_Diff"){
          lm_profiles<- read.csv(inFile$datapath,
                                 header = input$header_Diff,
                                 sep = input$sep_Diff,
                                 row.names = 1)
        }
        
        #transpose table if samples are in the rows 
        if(input$samp_Diff == "row_Diff"){
          lm_profiles = read.csv(inFile$datapath,
                                 header=input$header_Diff,
                                 sep=input$sep_Diff,
                                 row.names = 1)
          lm_profiles<-data.frame(t(lm_profiles))}
        
        # Get the classification from the main file:
        
        classification <- data.frame(samples = colnames(lm_profiles),
                                     group = as.character(unname(unlist(lm_profiles[1, ]))))
        
        # Transform the lipid mediator concentrations to numeric values so we can do the calculations: 
        
        lipid_m_profile <- as.data.frame(sapply(lm_profiles[-1, ], function (x) as.numeric(as.character(x))))
        rownames(lipid_m_profile) <- rownames(lm_profiles[-1, ])
        
        #transpose dataframe
        lm_profiles_t <- data.frame(t(lipid_m_profile))
        
        # Replace all the zero values by NA (This is for the sake of replacing zero/missing values)
        lm_profiles_t[lm_profiles_t == 0] <- NA
        
        
        #zero handling data
        cols<- 1:(ncol(lm_profiles_t))
        #replace zeros for each column to 1/5 the smallest value for each column
        lm_profiles_t[cols] <- lapply(lm_profiles_t[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
        
        #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
        lm_profiles_t<-as.matrix(lm_profiles_t)
        lm_profiles_t[is.infinite(lm_profiles_t)] <- NA
        min_index <- which.min(lm_profiles_t)
        zero_replace <- (lm_profiles_t[min_index]/5)
        lm_profiles_t <- as.data.frame(lm_profiles_t)
        lm_profiles_t[is.na(lm_profiles_t)] <- zero_replace
        
        # MULTINORMALITY CHECK:
        
        # For normality it takes a transpose dataframe as the one we made before.
        lm_profiles_mvn <- as.data.frame(t(lipid_m_profile))
        
        # Since normality does not work with variables with zero values in all its samples, we delete this variables. 
        lm_profiles_mvn <- lm_profiles_mvn[, colSums(lm_profiles_mvn != 0) > 0]
        
        # Also, the function doesn't work well when variables has the same values (duplicates), so we delete them. 
        lm_profiles_mvn <- lm_profiles_mvn[!duplicated(as.list(lm_profiles_mvn))]
        
        # This function calculates multivariate normality:
        testing_mult_norm <- mvn(data = lm_profiles_mvn, mvnTest = input$mvn, covariance = TRUE, tol = 1e-100)
        
        if (input$mvn == "mardia") {mnv_result <- as.character(testing_mult_norm$multivariateNormality$Result[3])}
        else {mnv_result <- as.character(testing_mult_norm$multivariateNormality$MVN)}
        
        
        # Comparison:
        
        # Please especify the names of the groups that you want to compare as it is in the classification file. Here,
        # the comparison is going to be made as "Group A vs Group B". Meaning, that group B, generally, is going to be
        # your control group. 
        
        group_a <- input$group_A
        group_b <- input$group_B    
        
        # Specify classification table:
        
        classification <- classification[(classification$group == group_a) | (classification$group == group_b), ]
        
        lm_profiles_t <- lm_profiles_t[rownames(lm_profiles_t) %in% as.character(classification$samples), ]
        
        # Check if both of your group follows a normal distribution:
        
        lm_profiles_group_a <- lipid_m_profile[ ,colnames(lipid_m_profile) %in% 
                                                  as.character(classification$samples[classification$group == group_a])]
        
        lm_profiles_group_b <- lipid_m_profile[ ,colnames(lipid_m_profile) %in% 
                                                  as.character(classification$samples[classification$group == group_b])]
        
        toptable <- data.frame(name = "delete",
                               FC = 1,
                               logFC = 1,
                               P.Value = 1,
                               adj.P.Val = 1,
                               test = "delete")
        
        
        
        for(i in 1:nrow(lipid_m_profile)) {
          
          if (mnv_result == "NO")  {
            
            # Calculates Mann_Whitney test:
            
            mw_test <- wilcox.test(lm_profiles_t[, i]~classification$group)
            
            # Calculates the Fold Change (A vs B):
            
            a <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_a]),i])
            b <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_b]),i])
            
            # Log FC:
            
            logFC <- log2(a) - log2(b)
            FC <- abs(a/b)
            
            # Put results together in a data frame:
            
            results_table <- data.frame(name = colnames(lm_profiles_t[i]),
                                        FC = format(FC, digits = 4),
                                        logFC = format(logFC, digits = 4),
                                        P.Value = format(mw_test$p.value, digits = 3),
                                        adj.P.Val = format(1, digits = 3),
                                        test = "Mann Whitney test")
            
            toptable <- rbind(toptable, results_table)
            
          }
          
          if (mnv_result == "YES") {
            
            t_test <- t.test(lm_profiles_t[, i]~classification$group)
            
            # Calculates the Fold Change (A vs B):
            
            a <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_a]),i])
            b <- mean(lm_profiles_t[which(rownames(lm_profiles_t) %in% classification$samples[classification$group == group_b]),i])
            
            # Log FC:
            
            logFC <- log2(a) - log2(b)
            FC <- abs(a/b)
            
            # Put results together in a data frame:
            
            results_table <- data.frame(name = colnames(lm_profiles_t[i]),
                                        FC = format(FC, digits = 4),
                                        logFC = format(logFC, digits = 4),
                                        P.Value = format(t_test$p.value, digits = 3),
                                        adj.P.Val = format(1, digits = 3),
                                        test = "t test")
            
            toptable <- rbind(toptable, results_table)
            
          }
          
        }
        
        toptable <- toptable[-1, ]
        toptable$adj.P.Val <- format(p.adjust(toptable$P.Value, "BH"), digits = 3)
        toptable$P.Value[toptable$P.Value == "NaN"] <- "1.000"
        toptable$adj.P.Val[toptable$adj.P.Val == "   NaN"] <- "1.000"
        toptable <- toptable[order(toptable$adj.P.Val), ]
        
        # #change the name to be in correct format 
        toptable$name[grep("^X", toptable$name)] <- paste(gsub("^X", "", toptable$name[grep("^X", toptable$name)]), "", sep = "")
       
        toptable$name <- gsub("\\.", "-", toptable$name) # Replace "." with "-" when required 
        
        
        #Volcano Plot 
        toptable$logFC<-as.numeric(toptable$logFC)
        toptable$adj.P.Val<-as.numeric(toptable$adj.P.Val)
        p_cutoff <- input$Sig_DiffPlot
        fc_cutoff <- input$FC_DiffPlot
        topN <- input$TopN_DiffPlot
        
        VolPlot<- toptable %>% 
          mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
          mutate(Rank = 1:n(), Label = ifelse(Rank < topN, name,"")) %>% 
          ggplot(aes(x = logFC, y = -log(adj.P.Val), col=Significant, label=Label)) + geom_point(size = 3) 
        
        VolPlot<-VolPlot + geom_text_repel(col="black", size = 5) +labs(x=" LogFC", y="Adjusted p-value" ) + 
          theme_minimal(base_size = 16) +
          scale_color_manual(values = c('black', "red"))
        
        # OUTPUTS:
        
        output$DiffTable <- renderDataTable(toptable,selection = "single", rownames = FALSE)
        output$DiffPlot<-renderPlot(VolPlot)
        
        output$DownloadDiffTable <- downloadHandler(filename = function(){"Differential_Analysis.csv"}, 
                                                    content = function(fname){
                                                      write.csv(toptable, fname)})
        
        output$downloadDiffPlot <- downloadHandler(filename = function(){paste("VolcanoPlot",'.png',sep='')},
                                                   content = function(file){
                                                     ggsave(file,plot=VolPlot, width = 12, height = 8)})
        
        
        ##calculating the pathway analysis information-----
        
        #open csv with the pathways for each metabolite
        Pathway <-read.csv(file = "coefficients/Pathways2.csv", header = TRUE, sep = ",")
        colnames(Pathway)<- c("Metabolite", "Pathway")
        
        #set it up like a dictionary 
        MetPath<-setNames(Pathway$Pathway, Pathway$Metabolite)
        toptable_P<-toptable
        toptable_P["Pathways"]<-"N/A"#column for pathways with N/A as default
        
        #loop to put the pathway for the metabolite 
        for (x in 1:length(toptable_P$name)){
          if (toptable_P$name[x] %in% Pathway$Metabolite){
            toptable_P$Pathways[x] <- MetPath[[toptable_P$name[x]]]
          }else{
            toptable_P$Pathways[x]<-toptable_P$Pathways[x]
          }
        }
        #remove all the rows missing the metabolite names
        toptable_Path<-toptable_P[toptable_P["Pathways"] != "N/A",]
        
        ###UPREGULATED PATHWAYS----
        toptable_UpReg<-toptable_Path[toptable_Path["logFC"] >0,]
        toptable_DownReg<-toptable_Path[toptable_Path["logFC"]<0,]
        
        UpReg_PathList<-unique(toptable_UpReg$Pathways)
        
        UpReg_Table<-data.frame(Pathway = "del", meanlogFC = 1)
        
        for (x in UpReg_PathList){
          x_UpReg <- toptable_UpReg[toptable_UpReg["Pathways"] == x, ]
          Avg_LogFC<-mean(x_UpReg$logFC)
          Table<-data.frame(Pathway = x, meanlogFC = Avg_LogFC)
          UpReg_Table<-rbind(UpReg_Table, Table)
        }
        
        UpReg_Table<-UpReg_Table[-1,]
        UpReg_Table<-UpReg_Table[order(-UpReg_Table$meanlogFC),]
        
        UpReg_Table$meanlogFC<-as.numeric(UpReg_Table$meanlogFC)
        UpReg_Table$Pathway<-as.character(UpReg_Table$Pathway)
        UpReg_Table$Pathway<-factor(UpReg_Table$Pathway, levels = UpReg_Table$Pathway)
        
        
        
        #output for title and plot 
        output$Diff_UpRegPath_Title <- renderUI({req(input$Go_Diff|input$Go_DiffPlot); h2("Upregulated Pathways Plot", align = "center") })
        
        output$Diff_UpRegPath_Plot<-renderPlotly(ggplotly(ggplot(UpReg_Table, aes(x=reorder(Pathway, -meanlogFC), meanlogFC)) + 
                                                            geom_bar(stat = "identity", fill = c("#00C19A")) +
                                                            labs(x = "Pathway", y = "Mean LogFC")+
                                                            theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        
        ##DOWNREGULATED PATHWAYS-----
        DownReg_PathList<-unique(toptable_DownReg$Pathways)
        DownReg_Table<-data.frame(Pathway = "del", meanlogFC = 1)
        
        for (x in DownReg_PathList){
          x_DownReg <- toptable_DownReg[toptable_DownReg["Pathways"] == x, ]
          Avg_LogFC<-mean(x_DownReg$logFC)
          Table<-data.frame(Pathway = x, meanlogFC = Avg_LogFC)
          DownReg_Table<-rbind(DownReg_Table, Table)
        }
        
        DownReg_Table<-DownReg_Table[-1,]
        DownReg_Table<-DownReg_Table[order(DownReg_Table$meanlogFC),]
        
        DownReg_Table$meanlogFC<-as.numeric(DownReg_Table$meanlogFC)
        DownReg_Table$Pathway<-as.character(DownReg_Table$Pathway)
        DownReg_Table$Pathway<-factor(DownReg_Table$Pathway, levels = DownReg_Table$Pathway)
        
        output$Diff_DownRegPath_Title <- renderUI({req(input$Go_Diff|input$Go_DiffPlot); h2("Downregulated Pathways Plot", align = "center") })
        output$Diff_DownRegPath_Plot<-renderPlotly(ggplotly(ggplot(DownReg_Table, aes(x=reorder(Pathway, meanlogFC), meanlogFC)) + 
                                                              geom_bar(stat = "identity", fill = c("#00C19A")) +
                                                              labs(x = "Pathway", y = "Mean LogFC")+
                                                              theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        
        
        
        ##P VALUES PATHWAYS
        PVal_toptable<-toptable_Path[toptable_Path["P.Value"] <0.05, ]
        PVal_PathList<-unique(PVal_toptable$Pathways)
        
        PVal_Table<-data.frame(Pathway = "del", Number_SPM = 1)
        
        for (x in PVal_PathList){
          x_PVal <- PVal_toptable[PVal_toptable["Pathways"] == x, ]
          Table<-data.frame(Pathway = x, Number_SPM = length(x_PVal$Pathways))
          PVal_Table<-rbind(PVal_Table, Table)
        }
        
        PVal_Table<-PVal_Table[-1,]
        PVal_Table<-PVal_Table[order(-PVal_Table$Number_SPM),]
        
        PVal_Table$Number_SPM<-as.numeric(PVal_Table$Number_SPM)
        PVal_Table$Pathway<-as.character(PVal_Table$Pathway)
        
        PVal_Table$Pathway<-factor(PVal_Table$Pathway, levels = PVal_Table$Pathway)
        
        output$Diff_PValPath_Title <- renderUI({req(input$Go_Diff|input$Go_DiffPlot); h2("Significant SPMs Pathways Plot", align = "center") })
        output$Diff_PValPath_Plot<-renderPlotly(ggplotly(ggplot(PVal_Table, aes(x=reorder(Pathway, -Number_SPM), Number_SPM)) + 
                                                           geom_bar(stat = "identity", fill = c("#00C19A")) +
                                                           labs(x = "Pathway", y = "Number of Signficant SPMs") +
                                                           theme(axis.text.x = element_text(angle = 45), panel.background = element_blank())))
        
        #make error messages blank if code works
        output$DiffTable_Error<- output$DiffPlot_Error<- output$Diff_UpRegPath_Plot_Error<-renderUI(return())
        
        
      },
      #what happens when there is an error in the code 
      error = function(e) {
        #notification box with specific error
        showNotification(paste0(e), type = 'err')
        #error messages for each tab explaining what could be wrong
        output$DiffTable_Error<- output$DiffPlot_Error<- output$Diff_UpRegPath_Plot_Error<-
          renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                               " - File Formats and seperator is correct and All inputs filled.", 
                               " - The correct Group B and Group B is selected",
                               " - The correct sample location is selected",
                               sep = "<br/>"))})
        #make table/plots outputs empty
        output$Diff_UpRegPath_Plot<-output$Diff_DownRegPath_Plot<-output$Diff_PValPath_Plot<-renderPlotly(return())
        output$DiffPlot<-plotOutput(return())
        output$DiffTable<-renderDataTable(return())
        output$Diff_UpRegPath_Title<-output$Diff_DownRegPath_Title<-output$Diff_PValPath_Title<-renderUI(return())
      }, silent=FALSE
    )
    
    
    
  })
  
  observe({
    req(input$Corr_file)
    
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {
        df_Corr<- read.csv(input$Corr_file$datapath,
                           header = input$header_Corr,
                           sep = input$sep_Corr)
        
        
        #transpose table if samples are at columns
        if(input$samp_Corr == "col_Corr"){
          df_Corr<-data.frame(t(df_Corr))
          df_samp<-rownames(df_Corr)
          df_Corr<-data.frame(cbind(df_samp, df_Corr))
          colnames(df_Corr)<- df_Corr[1,]
          df_Corr<- df_Corr[-1, ]}
        
        #get list of rows as options  
        Col_Names<-unique(as.character(unname(unlist(colnames(df_Corr)))))
        
        #if the colnames list is empty error messages shows
        if (Col_Names == ''){
          updateSelectInput(session,"Sample_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"Cond_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"GrpAstart_Cor", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"GrpAend_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"GrpBstart_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"GrpBend_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        }
        
        
        #updare input name 
        updateSelectInput(session, "Sample_Corr", choices = Col_Names, selected = Col_Names[[1]])
        updateSelectInput(session, "Cond_Corr", choices = Col_Names, selected = Col_Names[[1]])
        updateSelectInput(session, "GrpAstart_Corr", choices = Col_Names, selected = Col_Names[[1]])
        updateSelectInput(session, "GrpAend_Corr", choices = Col_Names, selected = Col_Names[[1]])
        updateSelectInput(session, "GrpBstart_Corr", choices = Col_Names, selected = Col_Names[[1]])
        updateSelectInput(session, "GrpBend_Corr", choices = Col_Names, selected = Col_Names[[1]])
        #make error messages blank
        output$CorrTable_Error<- output$CorrPlot_Error<-renderUI(return())
      }, 
      
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Sample_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"Cond_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"GrpAstart_Cor", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"GrpAend_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"GrpBstart_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"GrpBend_Corr", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
  })
  

  observe({
    
    #when Corr file for Group A is uploaded automatically get column names 
    req(input$Corr_fileA)
    
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {
        df_Corr<- read.csv(input$Corr_fileA$datapath,
                           header = input$header_CorrA,
                           sep = input$sep_CorrA)
        
        #transpose table if samples are at columns
        if(input$samp_CorrA == "col_CorrA"){
          df_Corr<-data.frame(t(df_Corr))
          df_samp<-rownames(df_Corr)
          df_Corr<-data.frame(cbind(df_samp, df_Corr))
          colnames(df_Corr)<- df_Corr[1,]
          df_Corr<- df_Corr[-1, ]}
        
        
        
        #get list of rows as options  
        Col_Names<-unique(as.character(unname(unlist(colnames(df_Corr)))))
        
        #if statement if the list of column names is empty
        if (Col_Names == ''){
          updateSelectInput(session,"Cond_CorrA", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          
        }
        
        #updare input name 
        updateSelectInput(session, "Cond_CorrA", choices = Col_Names, selected = Col_Names[[1]])
        
        #make error messages blank
        output$CorrTable_Error<- output$CorrPlot_Error<-renderUI(return())
        
      }, 
      
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Cond_CorrA", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
    
  })
  
  observe({
    
    #when Corr file for Group B is uploaded automatically get column names 
    req(input$Corr_fileB)
    
    #Error handling using trycatch
    tryCatch(
      
      #try to run the code to open file and get list of mediators 
      {
        
        df_Corr<- read.csv(input$Corr_fileB$datapath,
                           header = input$header_CorrB,
                           sep = input$sep_CorrB)
        
        #transpose table if samples are at columns
        if(input$samp_CorrB == "col_CorrB"){
          df_Corr<-data.frame(t(df_Corr))
          df_samp<-rownames(df_Corr)
          df_Corr<-data.frame(cbind(df_samp, df_Corr))
          colnames(df_Corr)<- df_Corr[1,]
          df_Corr<- df_Corr[-1, ]}
        
        #get list of rows as options  
        Col_Names<-unique(as.character(unname(unlist(colnames(df_Corr)))))
        
        #if statement if the list of column names is empty
        if (Col_Names == ''){
          updateSelectInput(session,"Cond_CorrB", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          
        }
        
        #updare input name 
        updateSelectInput(session, "Cond_CorrB", choices = Col_Names, selected = Col_Names[[1]])
        
        #make error messages blank
        output$CorrTable_Error<- output$CorrPlot_Error<-renderUI(return())
        
      }, 
      
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Cond_CorrB", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
  })
  
  #when Go button for correlation test is clicked 
  observeEvent(input$Go_Corr, {
    
    tryCatch(
      
      #try the code to see if it works 
      {
        if(input$Corr_OneTable == TRUE){
          req(input$Corr_file)
          
          #form dataframe with file
          df_Corr<- read.csv(input$Corr_file$datapath,
                             header = input$header_Corr,
                             sep = input$sep_Corr)
          #transpose table if samples are at columns
          if(input$samp_Corr == "col_Corr"){df_Corr<-data.frame(t(df_Corr))}
          inputCond<- as.character(input$Cond_Corr)
        }
        
        else{
          
          #dataframe for Group A-----
          req(input$Corr_fileA)
          
          df_Corr<- read.csv(input$Corr_fileA$datapath,
                             header = input$header_CorrA,
                             sep = input$sep_CorrA)
          #transpose table if samples are at columns
          if(input$samp_CorrA == "col_CorrA"){df_Corr<-data.frame(t(df_CorrA))}
          
          inputCond<- input$Cond_CorrA
        }
        
        
        
        #get list of conditions as options
        conditions<-unique(as.character(unname(unlist(df_Corr[,inputCond]))))
        
        #if the conditions list is empty
        if (conditions ==''){
          updateSelectInput(session,"Corr_Table", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"Corr_Plot", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        }
        
        #update table selection list based on options
        updateSelectInput(session, "Corr_Table", choices = conditions, selected = conditions[1])
        
        #update plot selection list based on options
        updateSelectInput(session, "Corr_Plot", choices = conditions, selected = conditions[1])
        
        #make error messages blank
        output$CorrTable_Error<- output$CorrPlot_Error<-renderUI(return())
        
      },
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Corr_Table", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"Corr_Plot", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)

    

  })
  

  
  observeEvent(input$Go_CorrTable|input$Go_CorrPlot, {
    
    if(input$Corr_OneTable == TRUE){
      req(input$Corr_file)}
    else{req(input$Corr_fileA)}
    
    #error handling trying to see if code would work
    tryCatch(
      
      #try to run the code 
      {
        
        if(input$Corr_OneTable == TRUE){
          req(input$Corr_file)
          
          #form dataframe with file 
          df_Corr<- read.csv(input$Corr_file$datapath,
                             header = input$header_Corr,
                             sep = input$sep_Corr)
          
          #transpose table if samples are at columns
          if(input$samp_Corr == "col_Corr"){df_Corr<-data.frame(t(df_Corr))}
          
          #### include option if there multiple condiditons or one conditions
          
          
          #get the column number for each of the selected column names 
          inputCond<-match(as.character(input$Cond_Corr), names(df_Corr))
          SampleCorr<-match(as.character(input$Sample_Corr), names(df_Corr))
          GrpAstartCorr<-match(as.character(input$GrpAstart_Corr), names(df_Corr))
          GrpAendCorr<-match(as.character(input$GrpAend_Corr), names(df_Corr))
          GrpBstartCorr<-match(as.character(input$GrpBstart_Corr), names(df_Corr))
          GrpBendCorr<-match(as.character(input$GrpBend_Corr), names(df_Corr))
          
          #zero handling 
          df_Corr_Zero<-df_Corr[,-c(SampleCorr, inputCond)]
          #ensure dataframe is numeric
          df_Corr_Zero<-data.frame(sapply(df_Corr_Zero, function(x) as.numeric(as.character(x))))
          
          #zero handling data
          cols<- 1:(ncol(df_Corr_Zero))
          #replace zeros for each column to 1/5 the smallest value for each column
          df_Corr_Zero[cols] <- lapply(df_Corr_Zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
          
          #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
          df_Corr_Zero<-as.matrix(df_Corr_Zero)
          df_Corr_Zero[is.infinite(df_Corr_Zero)] <- NA
          min_index <- which.min(df_Corr_Zero)
          zero_replace <- (df_Corr_Zero[min_index]/5)
          df_Corr_Zero <- as.data.frame(df_Corr_Zero)
          df_Corr_Zero[is.na(df_Corr_Zero)] <- zero_replace
          
          df_Corr<-cbind(df_Corr[,c(SampleCorr, inputCond)], df_Corr_Zero)
          
          #get datatables for each group 
          GrpA_Corr<-df_Corr[,c(SampleCorr, inputCond, (GrpAstartCorr:GrpAendCorr))]
          GrpB_Corr<-df_Corr[,c(SampleCorr, inputCond, (GrpBstartCorr:GrpBendCorr))]
          
          
        }
        
        else{
          
          #dataframe for Group A-----
          req(input$Corr_fileA)
          
          GrpA_Corr<- read.csv(input$Corr_fileA$datapath,
                               header = input$header_CorrA,
                               sep = input$sep_CorrA)
          
          # #transpose table if samples are at columns
          # if(input$samp_CorrA == "col_CorrA"){GrpA_Corr<-data.frame(t(GrpA_Corr))}
          
          #dataframe for Group B-----
          req(input$Corr_fileB)
          
          GrpB_Corr<- read.csv(input$Corr_fileB$datapath,
                               header = input$header_CorrB,
                               sep = input$sep_CorrB)
          inputCond<-match(as.character(input$Cond_CorrA), names(GrpA_Corr))
          
          # #transpose table if samples are at columns
          if(input$samp_CorrB == "col_CorrA"){GrpB_Corr<-data.frame(t(GrpB_Corr))}
        }
        
        #empty list for results in loop
        corrlist<-list()
        corrlist<-list()
        adjplist<-list()
        pos = 1
        
        #for loop to do the spearmans correlation for each condition
        
        for (cond in unique(as.character(unname(unlist(GrpA_Corr[,inputCond]))))) {
          
          
          #get the rows that match condition
          GrpA_Cond<-GrpA_Corr[GrpA_Corr[,2] == cond, ]
          GrpB_Cond<-GrpB_Corr[GrpB_Corr[,2] == cond, ]
          
          #make first row the rownames
          rownames(GrpA_Cond)<-GrpA_Cond[,1]
          rownames(GrpB_Cond)<-GrpB_Cond[,1]
          
          #remove excess columns
          GrpA_Cond[, c(1,2)]<-NULL
          GrpB_Cond[, c(1,2)]<-NULL
          
          #spearmans correlation
          cor_final<-cor(GrpA_Cond, GrpB_Cond, method = "spearman")
          
          #ad
          #p-values table
          pvalue_vector <- NULL # Creates empty vector
          pvalue_df <- data.frame(GrpA_Cond = colnames(GrpA_Cond)) # Creates p value data frame
          padj_value_df <- data.frame(GrpA_Cond =  colnames(GrpA_Cond)) # Creates adj p value data frame
          
          #populates p-value dataframe
          for (i in 1:ncol(GrpB_Cond)) {
            x <- GrpB_Cond[, i]
            for (j in 1:ncol(GrpA_Cond)) {
              y <- GrpA_Cond[, j]
              correlation <- cor.test(x, y, method = "spearman", exact = TRUE)
              pvalue_vector[j] <- correlation$p.value
            }
            padjust_vector <- p.adjust(pvalue_vector, method = "BH")
            pvalue_df[, colnames(GrpB_Cond[i])] <- pvalue_vector
            padj_value_df[, colnames(GrpB_Cond[i])] <- padjust_vector
          }
          
          # Update table with the p values associated to every correlation.
          rownames(pvalue_df) <- pvalue_df$GrpA_Cond
          pvalue_df$GrpA_Cond <- NULL
          pvalue_df <- as.matrix(pvalue_df)
          
          
          # Update table with the adj p values associated to every correlation.
          rownames(padj_value_df) <- padj_value_df$GrpA_Cond
          padj_value_df$GrpA_Cond <- NULL
          padj_value_df <- as.matrix(padj_value_df)
          
          #add the dataframe from loop to list
          corrlist[[pos]]<-cor_final
          adjplist[[pos]]<-padj_value_df
          corrlist[[pos]]<-cor_final
          
          
          #move to next position in list
          pos<-(pos + 1)
        }
        
        #list of conditions
        cond_list<-list(unique(GrpA_Corr[,inputCond]))
        #list of numbers as long as the conditions list 
        cond_num<-list(seq(1, (lengths(cond_list))))
        #create a dataframe with the conditions and associated number 
        cond_input<-cbind(data.frame(cond_list), data.frame(cond_num))
        colnames(cond_input)<-c("Conditions", "Key")
        
        
        #get the number once the condiition is inputed
        cond_val<-cond_input$Key[cond_input$Conditions==input$Corr_Table]
        cond_val2<-cond_input$Key[cond_input$Conditions==input$Corr_Plot]
        
        #test to see if it worked by outputing one table 
        output$CorrTable <- renderDataTable(corrlist[[cond_val]], caption = input$Corr_Table,
                                            options = list(paging = TRUE, pagelength = 20, 
                                                           scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE), 
                                            selection = "single", rownames = TRUE)
        output$DownloadCorrTable <- downloadHandler(filename = function(){paste(input$Corr_Table,"_Correlation_Analysis.csv","")}, 
                                                    content = function(fname){
                                                      write.csv(data.frame(corrlist[[cond_val]]), fname)})
        
        output$CorrPlot<-renderPlot(corrplot::corrplot(corrlist[[cond_val2]], method = "circle", p.mat = adjplist[[cond_val2]], 
                                                       sig.level = input$Sig_CorrPlot, insig = "blank", tl.col = "black", 
                                                       tl.srt = 45, cl.pos = "r"))
        
        #download function as pdf 
        output$DownloadCorrPlot <- downloadHandler(filename = function(){paste(input$Corr_Plot,"_Correlation_Plot.pdf","")}, 
                                                   content = function(fname){
                                                     pdf(fname, width = 5, height = 5)
                                                     corrplot::corrplot(corrlist[[cond_val2]], method = "circle", p.mat = adjplist[[cond_val2]], 
                                                                        sig.level = input$Sig_CorrPlot, insig = "blank", tl.col = "black", 
                                                                        tl.srt = 45, cl.pos = "r", tl.cex = 0.5)
                                                     dev.off()
                                                   })
        #make error messages blank
        output$CorrTable_Error<- output$CorrPlot_Error<-renderUI(return())
        
      },
      #what happens when there is an error in the code 
      error = function(e) {
        #notification box with specific error
        showNotification(paste0(e), type = 'err')
        #error messages for each tab explaining what could be wrong
        output$CorrTable_Error<- output$CorrPlot_Error<-
          renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                               " - File Formats and seperator is correct and All inputs filled.", 
                               " - The correct sample location is selected",
                               " - The first check box is only deselected if the groups are in 2 tables",
                               sep = "<br/>"))})
        #make table/plots outputs empty
        
        output$CorrPlot<-renderPlot(return())
        output$CorrTable<-renderDataTable(return())
        
      }, silent=FALSE
    )
    
    
   
  })


  



  observe({
    
    req(input$LM_ML_File)
    
    
    tryCatch(
      
      #try the code to see if it works 
      {
        
        lm_profile = read.table(input$LM_ML_File$datapath, 
                                sep=input$sep_LM_ML,
                                header = TRUE,
                                row.names = 1,
                                stringsAsFactors = FALSE)
        
        if (input$samp_ML == "col_ML"){
          
          lm_profile <- data.frame(t(lm_profile))}
        
        
        #get the column names and update 
        Col_NamesML<-unique(as.character(unname(unlist(colnames(lm_profile)))))
        
        #if error in file and column names are empty
        if (Col_NamesML == ""){
          updateSelectInput(session,"Group_ML", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        }
        
        updateSelectInput(session, "Group_ML", choices = Col_NamesML, selected = Col_NamesML[[1]])
        
        
        
        
        # Separates the profiles data by lipid mediators types:
        
        
        name_List <- unique(unname(unlist(lm_profile[1, -1])))
        
        name_List<-c("All_LM", name_List)
        
        #update the names of the models available in the datafile 
        updateSelectInput(session, "RF_Mod_Num", choices = name_List)
        updateSelectInput(session, "RF_Mod_Num1", choices = name_List)
        updateSelectInput(session, "XGB_Mod_Num", choices = name_List)
        updateSelectInput(session, "XGB_Mod_Num1", choices = name_List)
        updateSelectInput(session, "SVM_Mod_Num", choices = name_List)
        updateSelectInput(session, "LA_Mod_Num", choices = name_List)
        updateSelectInput(session, "BC_Mod_Num", choices = name_List)
        
      },
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Group_ML", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
    
  })
  
  observeEvent(input$Go_ML, {
    
    req(input$LM_ML_File)
    
    #Error handoling 
    tryCatch(
      {
      
      # Call the lipid mediator profiling file variant:
      inFile <- input$LM_ML_File
      
      print(inFile)
      
      if(is.null(inFile))
        
        return(NULL)
      
      
      set.seed(415) # To get same results even with the random part.
      options(digits = 3) # To get only part of the decimals. 
      
      # Open the file:
      
      
      lm_profile = read.table(inFile$datapath, 
                              sep=input$sep_LM_ML,
                              header = TRUE,
                              row.names = 1,
                              stringsAsFactors = FALSE)
      
      if (input$samp_ML == "col_ML"){
        
        lm_profile <- data.frame(t(lm_profile))
        
        #remove spaces from column names 
        colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
        
        for (x in 1:ncol(lm_profile)){
          if (substr(colnames(lm_profile)[x],1,1) %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
            colnames(lm_profile)[x]= paste("X", colnames(lm_profile)[x], sep = '')
          }
          else{colnames(lm_profile)[x] = colnames(lm_profile)[x]}
        }
        
      }
      
      
      Met_Rename<-c()
      #remove uncessary stull from names 
      colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
      
      
      for (x in 1:ncol(lm_profile)){
        if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
          Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
        }else{
          Met_Rename[x] <-colnames(lm_profile)[x]
        }
      }
      colnames(lm_profile)<-Met_Rename
      ##put coefficents into file!! 
      #lis of the unique substrates for the plots 
      name_List <- unique(unname(unlist(lm_profile[1, -1])))
      name_List<-c("All_LM", name_List)
      
      # Separates the profiles data by lipid mediators types:
      
      substrates <- unique(unname(unlist(lm_profile[1, -1])))
      
      
      
      # Creates data frames for each subset:
      
      dataframes_list <- list("ALL LM" = lm_profile[, -1])
      
      for (i in 1:length(substrates)) {
        substrate <- lm_profile[ ,lm_profile[1, ]== substrates[[i]]]
        assign(substrates[[i]], substrate)
        
      }
      
      if (exists('DHA') == TRUE) {dataframes_list[["DHA"]] <- DHA}
      if (exists('n3DPA') == TRUE) {dataframes_list[["n3DPA"]] <- n3DPA}
      if (exists('EPA') == TRUE) {dataframes_list[["EPA"]] <- EPA}
      if (exists('AA') == TRUE) {dataframes_list[["AA"]] <- AA}
      
      
      # Creates the data frame that is going to collect all the info regarding the models:
      
      final_table <- data.frame(machine_learning = "del",
                                groups = "del",
                                percentage_accuracy = 1,
                                sensitivity = 1,
                                specificity = 1,
                                TP_per = 1,
                                FP_per = 1,
                                TN_per = 1,
                                FN_per = 1,
                                stringsAsFactors = FALSE)
      
      #dataframe with top 5 models for XGBoost
      XGB_Model_Table <- data.frame(model = "del",
                                    Accuracy = 1,
                                    specificity = 1,
                                    sensitivity = 1,
                                    TP_per = 1, 
                                    FP_per = 1,
                                    TN_per = 1,
                                    FN_per = 1, 
                                    nrounds = 1,
                                    eta= 1,
                                    max_depth= 1,
                                    gamma= 1,
                                    min_child_weight= 1
      ) 
      
      #selecting tests
      RF_ML <- "RF_ML" %in% input$Model_ML
      XGB_ML <-"XGB_ML" %in% input$Model_ML
      SVM_ML <- "SVM_ML" %in% input$Model_ML
      BC_ML <- "BC_ML" %in% input$Model_ML
      LA_ML <- "LA_ML" %in% input$Model_ML
      
      #progress bar
      #withProgress(message = 'Building Models', style = style, value = 0.1, {
      # Sys.sleep(0.25)})
      
      
      withProgress(message = 'Building Models', value = 0, {
        
        for (j in 1:length(dataframes_list)) {
          
          lm_profiles <- dataframes_list[[j]]
          
          Fatty_Acid<-lm_profiles[1,]
          
          # Save all the values as numeric:
          lm_profile_number <- sapply(lm_profiles[-1, ], function(x) as.numeric(x))
          row.names(lm_profile_number) <- row.names(lm_profiles[-1, ])
          
          # Scale the data because is better when you are running Machine Learning models:
          lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))
          
          # If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors
          # replace the NA for zeros. 
          lm_profiles_scale[is.na(lm_profiles_scale)] <- 0
          
          #--------------------
          ### change 0 to 1/5 of lowest value 
          
          
          lm_profile_zero<-lm_profiles_scale[-1,]
          
          
          
          #get vector the same length as number of columns 
          cols<- 1:(ncol(lm_profile_zero))
          #replace zeros for each column to 1/5 the smallest value for each column 
          
          lm_profile_zero[cols] <- lapply(lm_profile_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
          
          #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
          lm_profile_zero<-as.matrix(lm_profile_zero)
          lm_profile_zero[is.infinite(lm_profile_zero)] <- NA
          min_index <- which.min(lm_profile_zero)
          zero_replace <- (lm_profile_zero[min_index]/5)
          lm_profile_zero <- as.data.frame(lm_profile_zero)
          lm_profile_zero[is.na(lm_profile_zero)] <- zero_replace  
          
          
          
          #add first row to dataframe 
          lm_profiles_scale<-rbind(lm_profiles_scale[1,] ,lm_profile_zero)
          
          
          
          
          
          #---------------------
          # Add the classification variable to the data frame (Responder and non responder):
          ResponseName<-as.character(input$Group_ML)
          lm_profiles_scale$responses <- factor(lm_profile[-1, ResponseName])
          
          # Make sure that column names do not represent a problem to randomForest making them a valid name to R.
          names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))
          
          #---> MACHINE LEARNING (randomForest R):
          if (RF_ML){
            # In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model. 
            oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.
            for (mtry in 1:(ncol(lm_profiles_scale) - 1)) { 
              rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                                    importance = TRUE, ntree = 10000)
              oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
            }
            # Define the best mtry according to the best prediction value. 
            final_mtry <- which.max(oob_error)
            
            # Run the model again with the right mtry value. 
            rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                                 importance = TRUE, ntree = 10000)
            # Get the confusion matrix of the model, sensitivity and specificity: 
            confusion_lm_profiles <- as.data.frame(rf_lm_profiles_final$confusion)
            confusion_lm_profiles[is.na(confusion_lm_profiles)] <- 0
            
            # Calculates sensitivity, specificity and AUC.
            sensitivity_lm_profiles <- confusion_lm_profiles[2, 2]/(confusion_lm_profiles[2, 2] + confusion_lm_profiles[2, 1])
            specificity_lm_profiles <- confusion_lm_profiles[1, 1]/(confusion_lm_profiles[1, 1] + confusion_lm_profiles[1, 2])
            
            # Final table for random forest:
            no_oob_error_table <- data.frame(machine_learning = "RF",
                                             groups = names(dataframes_list)[j],
                                             percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[10000])*100),
                                             sensitivity = sensitivity_lm_profiles,
                                             specificity = specificity_lm_profiles,
                                             TP_per = (1 - confusion_lm_profiles[2, 3])*100,
                                             FP_per = confusion_lm_profiles[1, 3]*100,
                                             TN_per = (1 -confusion_lm_profiles[1, 3])*100,
                                             FN_per = confusion_lm_profiles[2, 3]*100,
                                             stringsAsFactors = FALSE)
            final_table <- rbind(final_table, no_oob_error_table)
            # Number of trees plot:
            tree_table <- data.frame(Ensemble = c(1:10000),
                                     OBB = rf_lm_profiles_final$err.rate[, 1],
                                     AvgAcc = 100 - ((rf_lm_profiles_final$err.rate[, 1])*100))
            ntree_plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) +
              geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
              ggtitle(names(dataframes_list)[j]) +
              scale_x_continuous(name = "Trees") +
              scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
              theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 25, colour = "black"),
                    axis.title.x = element_text(size = 25, colour = "black"),
                    axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                    legend.position = ("none"),
                    panel.background = element_rect(fill = "white"),
                    axis.ticks.length = unit(0.4, "cm"))
            assign(paste("ntree_", names(dataframes_list)[j], sep = ""), ntree_plot)
            
            #get importance dataframe 
            rf_Importance<-as.data.frame(importance(rf_lm_profiles_final, type = 1))
            rf_Importance$lipid_mediators <- rownames(rf_Importance)
            
            Fatty_Acid <- 1:length(rf_Importance$lipid_mediators)
            cbind(rf_Importance, Fatty_Acid)
            
            
            # #fatty acid points 
            for (i in 1:length(rf_Importance$lipid_mediators)){
              rf_Importance$Fatty_Acid[i]<-SPMs_FA[rf_Importance$lipid_mediators[i]]
            }
            
            #change the names to be in the proper format 
            rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)] <- paste(gsub("^X", "'", rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)]), "'", sep = "")
            rf_Importance$lipid_mediators <- gsub("4$", "[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
            rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
            rf_Importance$lipid_mediators<- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required 
            rf_Importance$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", rf_Importance$lipid_mediators) # Transform n-3 DPA as a subscript: 
            
            # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
            # (this with the purpose of get the LM in decreasing order in the figure):
            
            rf_Importance <- rf_Importance[order(rf_Importance$MeanDecreaseAccuracy), ]
            
            rf_Importance$lipid_mediators <- factor( rf_Importance$lipid_mediators, 
                                                     levels = rf_Importance$lipid_mediators[
                                                       order(rf_Importance$MeanDecreaseAccuracy)])
            
            lipids <- as.character(rf_Importance$lipid_mediators)
            #add colors to SPM names
            lm_classes <- as.factor(rf_Importance$Fatty_Acid)
            names(lm_classes) <- rf_Importance$lipid_mediators
            
            dha_index<-which((lm_classes=="DHA")== TRUE)
            n_three_index<-which((lm_classes=="n3DPA")== TRUE)
            epa_index<-which((lm_classes=="EPA")== TRUE)
            aa_index<-which((lm_classes=="AA")== TRUE)
            
            lm_colors <- NULL
            lm_colors[dha_index] <- "blue"
            lm_colors[n_three_index] <- "brown"
            lm_colors[epa_index] <- "darkgoldenrod1"
            lm_colors[aa_index] <- "darkslategray"
            
            #rf_Importance1<- subset(rf_Importance, MeanDecreaseAccuracy >1)
            rf_VIP_plot <- ggplot(data = rf_Importance, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = Fatty_Acid)) + geom_point(size = 3) +
              scale_y_continuous(name = "Mean Decrease Accuracy") +
              labs(x = "Lipid Mediators", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
              scale_x_discrete(labels = parse(text = lipids)) + 
              scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                   "EPA" ="darkgoldenrod1", "AA" = "darkslategray"
              )) +
              coord_flip() + 
              theme(axis.title = element_text(size = 20),
                    axis.text.x  = element_text(size = 15, hjust = 0.5),
                    axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                    legend.position = "top",
                    aspect.ratio = 2/1,
                    legend.title = element_text(size = 20),
                    legend.text  = element_text(size = 15),
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
            
            assign(paste("rf_VIP_", names(dataframes_list)[j], sep = ""), rf_VIP_plot)
            assign(paste("rf_Mod_", names(dataframes_list)[j], sep = ""), rf_lm_profiles_final)
            
            
          }
          
          #---> MACHINE LEARNING (EXTREME GRADIENT BOOSTING)
          if (XGB_ML){
            df<-lm_profiles_scale
            #need to convert responses to integer of 0 and 1
            
            df$responses<-as.factor(df$responses)
            responses<-df$responses 
            label<-as.integer(df$responses) -1
            
            #remove labels from dataframe
            df$responses<-NULL
            df<-sapply(df, function(x) as.numeric(x))
            
            #
            #split it into training and testing dataset with 75/20
            n = nrow(df) #get number of rows for dataframe
            train.index = sample(n,floor(0.7*n)) #randomize rows to get dataset 
            train.data = as.matrix(df[train.index,])
            train.label = label[train.index]
            test.data = as.matrix(df[-train.index,])
            test.label = label[-train.index]
            
            #transform dataset to xgb matrix 
            xgb.train = xgb.DMatrix(data=train.data,label=train.label)
            xgb.test = xgb.DMatrix(data=test.data, label=test.label)
            
            #create grid with all the options for xgboost
            searchGridSubCol <- expand.grid(max_depth = c(1, 2, 3, 4, 6),
                                            gamma_val = c(0.4, 1, 1.5, 2, 2, 5), 
                                            eta = c(0.01, 0.02, 0.03, 0.4, 0.05),
                                            min_child = c(1, 2, 3, 4, 5)
            )
            
            #run through each combination of the grid 
            system.time(
              ErrorsHyperparameters <- apply(searchGridSubCol, 1, function(parameterList){
                
                #Extract Parameters to test
                currentDepth <- parameterList[["max_depth"]]
                currentEta <- parameterList[["eta"]]
                currentGamma <- parameterList[["gamma_val"]]
                currentMinChild <- parameterList[["min_child"]]
                
                #run selected parameter through the model 
                xgboostModelCV <- xgb.cv(data =  xgb.train, nrounds = 10000, nfold = 5, showsd = TRUE, 
                                         metrics = "error", verbose = TRUE, "eval_metric" = "error",
                                         "objective" = "binary:logistic", "max.depth" = currentDepth, "eta" = currentEta,       
                                         "subsample" = 1, "colsample_bytree" = 1,
                                         print_every_n = 10, booster = "gbtree",
                                         early_stopping_rounds = 10, "gamma" = currentGamma, "min_child_weight" = currentMinChild)
                
                #have error evaluation score as dataframe 
                xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
                test_error <- tail(xvalidationScores$test_error_mean, 1)
                train_error <- tail(xvalidationScores$train_error_mean,1)
                output <- return(c(test_error, train_error, currentDepth, currentEta, currentGamma, currentMinChild))})) #add to final table called output
            
            
            output_Error <- as.data.frame(t(ErrorsHyperparameters))
            
            varnames <- c("TestError", "TrainError","Depth", "eta", "gamma", "MinChild")
            names(output_Error) <- varnames #change the variable names 
            top_output<-output_Error[order(output_Error$TestError, output_Error$TrainError), ] #order it by the lowest error rate 
            
            #empty table to put the top 10 model information into 
            Final_Acc <- data.frame(model = "test",
                                    Accuracy = 0,
                                    specificity = 0, 
                                    sensitivity = 0,
                                    TP_per = 0, 
                                    FP_per = 0,
                                    TN_per = 0,
                                    FN_per = 0, 
                                    nrounds = 0,
                                    eta= 0,
                                    max_depth= 0,
                                    gamma= 0,
                                    min_child_weight=0
            )
            
            top<-1:50
            
            #loop to build models for the parameters with the lowest error rate 
            for (i in top){
              #parameter list with the lowest output 
              params <- list(booster = "gbtree", objective = "binary:logistic", eta= top_output[i,4], gamma= top_output[i,5],
                             max_depth= top_output[i,3], min_child_weight=top_output[i,6], subsample=1, 
                             colsample_bytree= 1)
              #Build model based on the training data 
              xgb<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
              
              #run model on test data = see how good it is at predicting
              pred<-predict(xgb, test.data)
              
              #obtain confusion matrix 
              pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1 
              pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
              pred<-as.numeric(pred) #makes it numeric
              test.label_num <-as.numeric(test.label) #makes it numeric
              response_List<-unique(responses)#get list of the 2 conditions 
              pred_y <- cut(pred, 2, label = c(response_List[1], response_List[2])) #converts to factor of responder or non responder 
              test_y <-cut(test.label_num, 2, label = c(response_List[1], response_List[2])) #convers to factor of responder or non responder 
              ConMatrix <- confusionMatrix(reference = test_y, data = pred_y) #build confusion matrix 
              TP = as.numeric(ConMatrix$table[2,2])
              FP = as.numeric(ConMatrix$table[2,1])
              TN = as.numeric(ConMatrix$table[1,1])
              FN = as.numeric(ConMatrix$table[1,2])
              
              #convert build table with accuracy for each model 
              Table <- data.frame(model = names(dataframes_list)[j],
                                  Accuracy = as.numeric(ConMatrix$overall[1])*100,
                                  specificity = as.numeric(ConMatrix$byClass[2]),
                                  sensitivity = as.numeric(ConMatrix$byClass[1]),
                                  TP_per = (TP/(TP+FN)*100), 
                                  FP_per = (FP/(FP+TN)*100),
                                  TN_per = (TN/(FP+TN)*100),
                                  FN_per = (FN/(TP+FN)*100), 
                                  nrounds = 10000,
                                  eta= top_output[i,4],
                                  max_depth= top_output[i,3],
                                  gamma= top_output[i,5],
                                  min_child_weight=top_output[i,6])
              Final_Acc<-rbind(Final_Acc, Table)
              
            }
            
            #remove the first row of random values 
            Final_Acc<-Final_Acc[-1,]
            
            #order the top 10 models based on highest accuracy, sensitivity, and specificity 
            Top_Acc<-Final_Acc[with(Final_Acc, order(-Accuracy, -sensitivity, -specificity)), ]
            
            # model<- c(paste(names(dataframes_list)[j],"_Model_1", sep = ""), 
            #            paste(names(dataframes_list)[j],"_Model_2", sep = ""),
            #            paste(names(dataframes_list)[j],"_Model_3", sep = ""),
            #            paste(names(dataframes_list)[j],"_Model_4", sep = ""),
            #            paste(names(dataframes_list)[j],"_Model_5", sep = ""))
            
            Top_Acc<-Top_Acc[1:5,]
            
            Top_Acc1 <- Top_Acc
            
            
            top5<-1:5
            for (i in top5){
              Top_Acc1[i,1] <-paste(names(dataframes_list)[j],"_Model", i,  sep = "")
            }
            
            
            # XGBoost Table:
            XGBoost_table <- data.frame(machine_learning = "XGB",
                                        groups = names(dataframes_list)[j],
                                        percentage_accuracy =Top_Acc[1,2],
                                        sensitivity = Top_Acc[1,3],
                                        specificity = Top_Acc[1,4],
                                        TP_per = Top_Acc[1,5], 
                                        FP_per = Top_Acc[1,6],
                                        TN_per = Top_Acc[1,7],
                                        FN_per = Top_Acc[1,8], 
                                        stringsAsFactors = FALSE)
            
            final_table <- rbind(final_table, XGBoost_table)
            
            XGB_Model_Table <- rbind(XGB_Model_Table, Top_Acc1)
            
            
            
            #preparing nrounds plots for top 5 models 
            iterations<-1:10000
            top_table<- data.frame(iterations)
            
            for (i in top5){
              
              #paramters for the i model 
              params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[i,10], gamma= Top_Acc[i,12], 
                             max_depth= Top_Acc[i,11], min_child_weight= Top_Acc[i,13], subsample= 1, colsample_bytree= 1)
              
              #run the cv to get the error rate for each nround 
              xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = 10000, nfold = 5, showsd = T, 
                               stratified = T, early_stop_round = 20, maximize = F, metrics = "error", verbose = FALSE)
              
              Table1<-data.frame( test = as.numeric(xgbcv$evaluation_log$test_error_mean))
              colnames(Table1)[1] <-paste(names(dataframes_list)[j],"_Model_", i, sep = "")
              #colnames(Table1)[2] <-paste(names(dataframes_list)[j],"Test_model_", i, sep = "")
              top_table<-cbind(top_table, Table1) #add to table
              
            }
            
            #plot using plotly 
            nround_plot<-plot_ly(data = top_table, x = ~iterations) 
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[2])), name = colnames(top_table)[2], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[3])), name = colnames(top_table)[3], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[4])), name = colnames(top_table)[4], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[5])), name = colnames(top_table)[5], mode = 'lines')
            nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[6])), name = colnames(top_table)[6], mode = 'lines')
            nround_plot <- nround_plot %>% layout(title = paste(names(dataframes_list)[j], "Plot", sep =" "), xaxis = list(title = 'Number of Rounds'), 
                                                  yaxis = list(title = 'Test Error Rate'))
            #makes download quality of plot better
            
            nround_plot <- nround_plot%>% config(toImageButtonOptions = list(format = "jpeg", width = 1500, height = 750))
            
            assign(paste("nround_", names(dataframes_list)[j], sep = ""), nround_plot)
            
            
            
            #importance plot 
            
            #build model for the top model 
            params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[1,10], gamma= Top_Acc[1,12], 
                           max_depth= Top_Acc[1,11], min_child_weight= Top_Acc[1,13], subsample= 1, colsample_bytree= 1)
            xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
            
            
            XGB_Importance = xgb.importance(model = xgb_Mod1) #important matrix
            
            Fatty_Acid <- 1:length(XGB_Importance$Feature)
            cbind(XGB_Importance, Fatty_Acid)
            
            #fatty acid points 
            for (i in 1:length(XGB_Importance$Feature)){
              XGB_Importance$Fatty_Acid[i]<-SPMs_FA[XGB_Importance$Feature[i]]
            }
            
            # #change the names to be in the proper format
            XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)] <- paste(gsub("^X", "'", XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)]), "'", sep = "")
            XGB_Importance$Feature <- gsub("4$", "[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 4
            XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature ) #B4 converted to subscript 4
            XGB_Importance$Feature  <- gsub("\\.", "-", XGB_Importance$Feature ) # Replace "." with "-" when required
            XGB_Importance$Feature <- gsub("n.3.DPA", "[n-3~DPA]", XGB_Importance$Feature ) # Transform n-3 DPA as a subscript:
            
            
            # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
            # (this with the purpose of get the LM in decreasing order in the figure):
            
            XGB_Importance <- XGB_Importance[order(XGB_Importance$Gain), ]
            
            XGB_Importance$Feature<- factor(XGB_Importance$Feature, 
                                            levels = XGB_Importance$Feature[
                                              order(XGB_Importance$Gain)])
            lipids <- as.character(XGB_Importance$Feature )
            
            #add colors to x axis labels
            lm_classes <- as.factor(XGB_Importance$Fatty_Acid)
            names(lm_classes) <- XGB_Importance$Feature
            
            dha_index<-which((lm_classes=="DHA")== TRUE)
            n_three_index<-which((lm_classes=="n3DPA")== TRUE)
            epa_index<-which((lm_classes=="EPA")== TRUE)
            aa_index<-which((lm_classes=="AA")== TRUE)
            
            lm_colors <- NULL
            lm_colors[dha_index] <- "blue"
            lm_colors[n_three_index] <- "brown"
            lm_colors[epa_index] <- "darkgoldenrod1"
            lm_colors[aa_index] <- "darkslategray"
            
            XGB_VIP_Plot <- ggplot(data = XGB_Importance, mapping = aes(x = Feature, y = Gain, color = Fatty_Acid)) + geom_point(size = 3) +
              scale_y_continuous(name = "Gain") +
              labs(x = "Features", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
              scale_x_discrete(labels = parse(text = lipids)) +
              scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                   "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
              coord_flip() +
              theme(axis.title = element_text(size = 20),
                    axis.text.x  = element_text(size = 15, hjust = 0.5),
                    axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                    legend.position = "top",
                    panel.background = element_rect(fill = "white"),
                    panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
            
            
            #save VIP plot
            assign(paste("XGB_VIP_", names(dataframes_list)[j], sep = ""), XGB_VIP_Plot)
            
            #Save Xgboost model 
            assign(paste("xgb_Mod_", names(dataframes_list)[j], sep = ""), xgb_Mod1)
          }
          #---> MACHINE LEARNING (Classyfire R): 
          
          # Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
          # and creates a novelty detection for the creation of the model. 
          
          # The idea is to create several models and see which one fits the best. The models will be based on the whole
          # lipid profiles and the different groups based on substrates. 
          
          # "cfBuild" to create the SVM:
          # Clasyfire requieres matrix: 
          if (SVM_ML){
            #convert dataframe to matrix
            lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])
            
            #building SVM model
            support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses, 
                                                bootNum = 70,ensNum = 70, cpus = 4) 
            #obtaining confusion matrix 
            conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale))
            
            # SVM table:
            Support_vector_table <- data.frame(machine_learning = "SVM",
                                               groups = names(dataframes_list)[j],
                                               percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                               sensitivity = conf_matrix[4, 3]/100,
                                               specificity = conf_matrix[1, 3]/100,
                                               TP_per = conf_matrix[4, 3],
                                               FP_per = conf_matrix[2, 3],
                                               TN_per = conf_matrix[1, 3],
                                               FN_per = conf_matrix[3, 3],
                                               stringsAsFactors = FALSE)
            final_table <- rbind(final_table, Support_vector_table)
            
            # Ensemble Plot:
            ensAcc   <- getAcc(support_lmprofiles_scale)$Test
            meanVal  <- ensAcc[1]
            for (i in 2:length(ensAcc)) {
              meanVal <- c(meanVal, mean(ensAcc[1:i]))
            } 
            ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc), 
                                        AvgAcc = meanVal)
            ensemble_plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) + 
              geom_point(aes(colour = AvgAcc), size = 5) + 
              geom_line(linetype = "dotted", size = 1) +
              ggtitle(names(dataframes_list)[j]) +
              scale_x_continuous(name = "Ensemble interaction") +
              scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
              theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 25, colour = "black"),
                    axis.title.x = element_text(size = 25, colour = "black"),
                    axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                    legend.position = ("none"),
                    panel.background = element_rect(fill = "white"),
                    axis.ticks.length = unit(0.4, "cm"))
            assign(paste("ensemble_", names(dataframes_list)[j], sep = ""), ensemble_plot)
            
            #Save SVM model 
            assign(paste("svm_Mod_", names(dataframes_list)[j], sep = ""), support_lmprofiles_scale)
          }
          
          #---> ELASTIC NET REGRESSION (caret R):
          if (LA_ML){
            # Get the explanatory variables as a matrix:
            explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
            
            # LASSO Analysis:
            model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet", 
                               trControl = trainControl("boot", number = 70))
            # Get the confusion matrix of the model:
            conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
            
            # Final Elastic net model table:
            net_table <- data.frame(machine_learning = "GLMNET",
                                    groups = names(dataframes_list)[j],
                                    percentage_accuracy = max(model_net$results$Accuracy)*100,
                                    sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                                    specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                                    TP_per = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                    FP_per = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                    TN_per = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                    FN_per = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                    stringsAsFactors = FALSE)
            final_table <- rbind(final_table, net_table)
            
            # Parameter tuning figure: 
            scaleFUN <- function(x) sprintf("%.2f", x)
            
            boot_net_plot <- ggplot(model_net, highlight = TRUE) +
              scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
              scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
              ggtitle(names(dataframes_list)[j]) +
              scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
              scale_shape_manual(values=c(16, 16, 16)) +
              labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
              guides(shape = FALSE) +
              theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                    axis.title.y = element_text(size = 25, colour = "black"),
                    axis.title.x = element_text(size = 25, colour = "black"),
                    axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                    axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
                    axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep.
                    legend.title = element_text(size = 25),
                    legend.text  = element_text(size = 20),
                    legend.key = element_rect(fill = NA),
                    legend.key.size = unit(1.3, "cm"), 
                    panel.background = element_rect(fill = "white"),
                    axis.ticks.length = unit(0.4, "cm"))
            assign(paste("net_", names(dataframes_list)[j], sep = ""), boot_net_plot)
            
            #rename and save elastic net model 
            assign(paste("la_Mod_", names(dataframes_list)[j], sep = ""), model_net)
          }
          
          #---> BAYESIAN MODEL (Caret R): 
          
          if (BC_ML){
            explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
            
            bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm", 
                              trControl = trainControl("boot", number = 70))
            conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)
            
            # Final Bayesian model table:
            bay_table <- data.frame(machine_learning = "BAYES",
                                    groups = names(dataframes_list)[j],
                                    percentage_accuracy = max(bayesian$results$Accuracy)*100,
                                    sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                                    specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                                    TP_per = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                    FP_per = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                    TN_per = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                    FN_per = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                    stringsAsFactors = FALSE)
            final_table <- rbind(final_table, bay_table)
            
            #rename and save bayesian model 
            assign(paste("bc_Mod_", names(dataframes_list)[j], sep = ""), bayesian)
          }
          
          incProgress(1/length(dataframes_list), detail = paste(names(dataframes_list)[j]))
        }
        
        final_table <- final_table[-1, ]
        table <- data.frame(Model = factor(final_table$groups, levels = unique(final_table$groups)),
                            methodology = final_table$machine_learning,
                            accuracy = round(final_table$percentage_accuracy, 0))
        accuracy <- ggplot(data = table, aes(x = Model, y = accuracy, fill = methodology)) +
          geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black", size = 2) +
          scale_fill_manual(values=c("dodgerblue2","firebrick2",'goldenrod1', 'lightslategray', "purple")) + 
          geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
                    aes(label = paste(accuracy, "%", sep = "")), size = 8) +
          scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20),
                             expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
          coord_cartesian(ylim = c(1, 120)) +
          theme(axis.title = element_text(size = 40),
                axis.title.x = element_blank(),
                axis.text.x  =  element_text(size = 40, hjust = 0.5, colour = "black"), # Put color to the labels
                axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
                axis.line = element_line(colour = 'black', size = 1.0), # Color and thickness of axis
                axis.ticks = element_line(colour = "black", size = 1.0), # Color and thickness of every axis sep. 
                panel.background = element_rect(fill = "white"),
                legend.title = element_blank(),
                legend.position = "top",
                legend.key.size = unit(1.3, "cm"), 
                legend.text  = element_text(size = 25),
                legend.spacing.x = unit(1, "cm"),
                axis.ticks.length = unit(0.4, "cm"))
        
        #Put figures and models together together: 
        if (RF_ML){
          rf_plot_list <- list("ALL LM" = `ntree_ALL LM`)
          rf_VIPplot_list <- list("ALL LM" = `rf_VIP_ALL LM`)
          rf_Mod_list<-list("ALL LM" = `rf_Mod_ALL LM`)
        }
        if (XGB_ML){
          XGB_VIPplot_list <- list("ALL LM" = `XGB_VIP_ALL LM`)
          xgb_Mod_list<-list("ALL LM" = `xgb_Mod_ALL LM`)}
        if (SVM_ML){
          svm_plot_list <- list("ALL LM" = `ensemble_ALL LM`)
          svm_Mod_list<-list("ALL LM" = `svm_Mod_ALL LM`)
        }
        if (LA_ML){
          net_plot_list <- list("ALL LM" = `net_ALL LM`)
          la_Mod_list<-list("ALL LM" = `la_Mod_ALL LM`)}
        if (BC_ML){bc_Mod_list<-list("ALL LM" = `bc_Mod_ALL LM`)}
        
        
        column <- 2
        if (exists('DHA') == TRUE) {
          if (RF_ML){
            rf_plot_list[["DHA"]] <- ntree_DHA
            rf_VIPplot_list[["DHA"]] <- rf_VIP_DHA
            rf_Mod_list[["DHA"]] <- rf_Mod_DHA
          }
          
          if (XGB_ML){
            XGB_VIPplot_list[["DHA"]] <- XGB_VIP_DHA
            xgb_Mod_list[["DHA"]] <- xgb_Mod_DHA
          }
          if (SVM_ML){
            svm_plot_list[["DHA"]] <- ensemble_DHA
            svm_Mod_list[["DHA"]] <- svm_Mod_DHA}
          if (LA_ML){
            net_plot_list[["DHA"]] <- net_DHA
            la_Mod_list[["DHA"]] <- la_Mod_DHA}
          if (BC_ML){bc_Mod_list[["DHA"]] <- bc_Mod_DHA}
          column <- column + 1}
        
        if (exists('n3DPA') == TRUE) {
          if (RF_ML){
            rf_plot_list[["n3DPA"]] <- ntree_n3DPA
            rf_VIPplot_list[["n3DPA"]] <- rf_VIP_n3DPA
            rf_Mod_list[["n3DPA"]] <- rf_Mod_n3DPA
          }
          if (XGB_ML){
            XGB_VIPplot_list[["n3DPA"]] <- XGB_VIP_n3DPA
            xgb_Mod_list[["n3DPA"]] <- xgb_Mod_n3DPA}
          if (SVM_ML){
            svm_plot_list[["n3DPA"]] <- ensemble_n3DPA
            svm_Mod_list[["n3DPA"]] <- svm_Mod_n3DPA}
          if (LA_ML){
            net_plot_list[["n3DPA"]] <- net_n3DPA
            la_Mod_list[["n3DPA"]] <- la_Mod_n3DPA}
          if (BC_ML){bc_Mod_list[["n3DPA"]] <- bc_Mod_n3DPA}
          
          column <- column + 1}
        
        if (exists('EPA') == TRUE) {
          if (RF_ML){
            rf_plot_list[["EPA"]] <- ntree_EPA
            rf_VIPplot_list[["EPA"]] <- rf_VIP_EPA
            rf_Mod_list[["EPA"]] <- rf_Mod_EPA
          }
          if (XGB_ML){
            XGB_VIPplot_list[["EPA"]] <- XGB_VIP_EPA
            xgb_Mod_list[["EPA"]] <- xgb_Mod_EPA}
          if (SVM_ML){
            svm_plot_list[["EPA"]] <- ensemble_EPA
            svm_Mod_list[["EPA"]] <- svm_Mod_EPA}
          if (LA_ML){
            net_plot_list[["EPA"]] <- net_EPA
            la_Mod_list[["EPA"]] <- la_Mod_EPA}
          if (BC_ML){bc_Mod_list[["EPA"]] <- bc_Mod_EPA}
          
          column <- column + 1}
        
        if (exists('AA') == TRUE) {
          if (RF_ML){
            rf_plot_list[["AA"]] <- ntree_AA
            rf_VIPplot_list[["AA"]] <- rf_VIP_AA
            rf_Mod_list[["AA"]] <- rf_Mod_AA}
          if (XGB_ML){
            XGB_VIPplot_list[["AA"]] <- XGB_VIP_AA
            xgb_Mod_list[["AA"]] <- xgb_Mod_AA}
          if (SVM_ML){
            svm_plot_list[["AA"]] <- ensemble_AA
            svm_Mod_list[["AA"]] <- svm_Mod_AA}
          if (LA_ML){
            net_plot_list[["AA"]] <- net_AA
            la_Mod_list[["AA"]] <- la_Mod_AA}
          if (LA_ML){bc_Mod_list[["AA"]] <- bc_Mod_AA}
          
          column <- column + 1}
        if (column >= 6) {column <- 4}
        
        #if (RF_ML){grid.arrange(grobs = rf_plot_list, ncol = column/2)}
        else if (SVM_ML){grid.arrange(grobs = svm_plot_list, ncol = column/2)}
        else if (LA_ML){grid.arrange(grobs = net_plot_list, ncol = column/2)}
        
        # Final table otputs :
        final_table$`% Accuracy Score` <- round(final_table$percentage_accuracy, 0)
        final_table$Sensitivity <- round(final_table$sensitivity, 2)
        final_table$Specificity <- round(final_table$specificity, 2)
        final_table$TP <- round(final_table$TP_per, 0)
        final_table$FP <- round(final_table$FP_per, 0)
        final_table$TN <- round(final_table$TN_per, 0)
        final_table$FN <- round(final_table$FN_per, 0)
        
        final_table <- final_table[, c(1, 2, 10:16)]
        colnames(final_table)[1] <- "Machine Learning Methodology"
        colnames(final_table)[2] <- "Model"
        final_table<-as.data.frame(final_table)
        
        
        ##outputs for accuracy tab 
        output$AccPlot_Title <- renderUI({req(input$Go_ML); h2("% Accuracy Score Figure for different ML models", align = "center") })
        output$AccTable_Title <- renderUI({req(input$Go_ML); h2("Model Table Summary", align = "center") })
        
        output$Accuracy_ML <- renderPlot({return(accuracy)})
        output$downloadAcc_ML_Plot <- downloadHandler(filename = function(){paste("Accuracy_Plot",'.png',sep='')},
                                                      content = function(file){
                                                        ggsave(file,plot= accuracy, width = 15, height = 10)})
        
        output$ML_Table <- renderDataTable({return(final_table)})
        output$downloadAcc_ML_Table <- downloadHandler(filename = function(){"Accuracy_Table.csv"}, 
                                                       content = function(fname){
                                                         write.csv(final_table, fname)})
        
        
        
        #outputs if Random forest model is run
        if (RF_ML){
          #ui outputs to get the page names 
          output$RF_Plot_Title <- renderUI({req(input$Go_ML); h2("Optimal Parameters RandomForest", align = "center") })
          output$RF_VIPPlot_Title <- renderUI({req(input$Go_ML); h2("Importance of Variance Plot", align = "center") })
          
          #random forest optimal parameters plots output and download button 
          output$RF_Plot <- renderPlot({return(grid.arrange(grobs = rf_plot_list, ncol = 1))})
          rf_plot_grid<-grid.arrange(grobs = rf_plot_list, ncol = 1)
          output$download_RF_Plot <- downloadHandler(filename = function(){paste("Random_Forest_Plot",'.png',sep='')},
                                                     content = function(file){
                                                       ggsave(file, plot= rf_plot_grid, width = 12, height = 20)})
          
          #random forest importance plot output and download button
          output$RF_VIPPlot <- renderPlot({return(grid.arrange(grobs = rf_VIPplot_list, ncol = 2))})
          output$download_RF_VIPPlot <- downloadHandler(filename = function(){paste("Random_Forest_VIP_Plot",input$RF_Mod_Num,'.png',sep='')},
                                                        content = function(file){
                                                          ggsave(file,plot= rf_VIPplot_list[[match(input$RF_Mod_Num, name_List)]], width = 12, height = 12)})
          
          #downloading model information for paramters tab and importance tab
          output$download_RF_Mod <- downloadHandler(filename = function(){paste("Random_Forest_Model_",input$RF_Mod_Num,'.rds',sep='')},
                                                    content = function(fname){
                                                      saveRDS(rf_Mod_list[[match(input$RF_Mod_Num, name_List)]], fname)})
          
          output$download_RF_Mod1 <- downloadHandler(filename = function(){paste("Random_Forest_Model_",input$RF_Mod_Num1,'.rds',sep='')},
                                                     content = function(fname){
                                                       saveRDS(rf_Mod_list[[match(input$RF_Mod_Num1, name_List)]], compress = TRUE, fname)})
          
        }
        
        #outputs for xgboost tabs 
        if (XGB_ML){
          #ui output page titles - only appears once action button clicked
          output$XGB_Plot_Title <- renderUI({req(input$Go_ML); h2("Top 5 models for XGBoost test error rate", align = "center") })
          output$XGB_VIPPlot_Title <- renderUI({req(input$Go_ML); h2("Importance of Variance Plot", align = "center") })
          
          
          #Parameters plots and tablefor xgboost parameters page 
          output$XGB_Plot1 <- renderPlotly({return(`nround_ALL LM`)})
          if (exists('DHA') == TRUE){output$XGB_Plot2 <- renderPlotly({return(nround_DHA)})}
          if (exists('n3DPA') == TRUE) {output$XGB_Plot3 <- renderPlotly({return(nround_n3DPA)})}
          if (exists('EPA') == TRUE) {output$XGB_Plot4 <- renderPlotly({return(nround_EPA)})}
          if (exists('AA') == TRUE) {output$XGB_Plot5 <- renderPlotly({return(nround_AA)})}
          XGB_Model_Table<-XGB_Model_Table[-1,]
          XGB_Model_Table[,-1] <-round(XGB_Model_Table[,-1],2)
          output$XGB_Models<-renderDataTable(XGB_Model_Table, rownames = FALSE,
                                             options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
          output$download_XGB_Table <- downloadHandler(filename = function(){"XGB_Model_Table.csv"}, 
                                                       content = function(fname){
                                                         write.csv(XGB_Model_Table, fname)})
          
          #xgboost importance plots output 
          output$XGB_VIPPlot <- renderPlot({return(grid.arrange(grobs = XGB_VIPplot_list, ncol = 2))})
          output$download_XGB_VIPPlot <- downloadHandler(filename = function(){paste("XGBoost_VIP_Plot",'.png',sep='')},
                                                         content = function(file){
                                                           ggsave(file,plot= XGB_VIPplot_list[[input$XGB_Mod_Num]], width = 12, height = 12)})
          
          #downloading model information for paramters tab and importance tab for cgboost
          output$download_XGB_Mod <- downloadHandler(filename = function(){paste("XGBoost_Model_",input$XGB_Mod_Num,'.rds',sep='')},
                                                     content = function(fname){
                                                       saveRDS(xgb_Mod_list[[match(input$XGB_Mod_Num, name_List)]], fname)})
          
          output$download_XGB_Mod1 <- downloadHandler(filename = function(){paste("XGBoost_Model_",input$XGB_Mod_Num1,'.rds',sep='')},
                                                      content = function(fname){
                                                        saveRDS(xgb_Mod_list[[match(input$XGB_Mod_Num1, name_List)]], fname)})
          
        }
        
        
        #all the outputs for svm model 
        if (SVM_ML){
          
          #ui output for SVM title page 
          output$SVM_Plot_Title <- renderUI({req(input$Go_ML); h2("Optimal Parameters SVM", align = "center") })
          
          #plot output and download for svm optima parameters 
          output$SVM_Plot <- renderPlot({return(grid.arrange(grobs = svm_plot_list, ncol = 1))})
          svm_plot_grid<-grid.arrange(grobs = svm_plot_list, ncol = 1)
          output$download_SVM_Plot <- downloadHandler(filename = function(){paste("SVM_Ensemble_Plot",'.png',sep='')},
                                                      content = function(file){
                                                        ggsave(file,plot= svm_plot_grid, width = 12, height = 20)})
          
          #download button for svm model
          output$download_SVM_Mod <- downloadHandler(filename = function(){paste("SVM_Model_",input$SVM_Mod_Num,'.rds',sep='')},
                                                     content = function(fname){
                                                       saveRDS(svm_Mod_list[[match(input$SVM_Mod_Num, name_List)]], fname)})
          
        }
        
        #outputs for elastic net regreession model 
        if (LA_ML){
          #Elastic Net title 
          output$LA_Plot_Title <- renderUI({req(input$Go_ML); h2("Optimal Parameters Elastic Net Regression", align = "center") })
          
          #plot output for elastic net regression 
          output$LA_Plot <- renderPlot({return(grid.arrange(grobs = net_plot_list, ncol = 1))})
          net_plot_grid<-grid.arrange(grobs = net_plot_list, ncol = 1)
          output$download_LA_Plot <- downloadHandler(filename = function(){paste("Elasti_Net_Plot",'.png',sep='')},
                                                     content = function(file){
                                                       ggsave(file,plot= net_plot_grid, width = 12, height = 20)})
          
          #download button for LA model
          output$download_LA_Mod <- downloadHandler(filename = function(){paste("ElasticNet_Model_",input$LA_Mod_Num,'.rds',sep='')},
                                                    content = function(fname){
                                                      saveRDS(la_Mod_list[[match(input$LA_Mod_Num, name_List)]], fname)})
        }
        
        #output for bayesian model 
        if (BC_ML){
          
          #download button for BC model
          output$download_BC_Mod <- downloadHandler(filename = function(){paste("BayesianGLM_Model_",input$BC_Mod_Num,'.rds',sep='')},
                                                    content = function(fname){
                                                      saveRDS(bc_Mod_list[[match(input$BC_Mod_Num, name_List)]], fname)})
        }
        
        output$Accuracy_ML_Error<- output$RF_Plot_Error<- output$RF_VIPPlot_Error<-output$XGB_Models_Error<-
          output$XGB_VIPPlot_Error<-output$SVM_Plot_Error<-output$LA_Plot_Error<-output$BC_Mod_Error<-renderUI(return())
        
        
        
      })
      
      
    },
    
    #what happens when there is an error in the code 
    error = function(e) {
      #notification box with specific error
      showNotification(paste0(e), type = 'err')
      #error messages for each tab explaining what could be wrong
      output$Accuracy_ML_Error<- output$RF_Plot_Error<- output$RF_VIPPlot_Error<-output$XGB_Models_Error<-
        output$XGB_VIPPlot_Error<-output$SVM_Plot_Error<-output$LA_Plot_Error<-output$BC_Mod_Error<-
        renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                             " - File Formats and seperator is correct and All inputs filled.", 
                             " - The correct Group column is selected and only 2 groups",
                             " - The correct sample location is selected",
                             " - The first row is the fatty acid family",
                             sep = "<br/>"))})
      #make table/plots outputs empty
      output$XGB_Plot1<-output$XGB_Plot2<-output$XGB_Plot3<-output$XGB_Plot4<-output$XGB_Plot5<-renderPlotly(return())
      output$Accuracy_ML<-output$RF_Plot<-output$RF_VIPPlot<-output$SVM_Plot<-output$LA_Plot<-
        renderPlot(return())
      output$ML_Table<-output$XGB_Models<-renderDataTable(return())
      output$AccPlot_Title<-output$AccTable_Title<-output$RF_Plot_Title<-output$RF_VIPPlot_Title<-output$XGB_Plot_Title<-
        output$XGB_VIPPlot_Title<-output$SVM_Plot_Title<-output$LA_Plot_Title<-
        renderUI(return())
    }, silent=FALSE)
    
    
    
 
  })
  
  observe({
    
    req(input$LM_ML_File2)
    
    tryCatch(
      
      #try the code to see if it works 
      {
        
        lm_profile = read.table(input$LM_ML_File2$datapath, 
                                sep=input$sep_LM_ML2,
                                header = TRUE,
                                row.names = 1,
                                stringsAsFactors = FALSE)
        
        if (input$samp_ML2 == "col_ML2"){
          
          lm_profile <- data.frame(t(lm_profile))
          
        }
        
        
        #get the column names and update 
        Col_NamesML<-unique(as.character(unname(unlist(colnames(lm_profile)))))
        
        #if error in file and column names are empty
        if (Col_NamesML == ""){
          updateSelectInput(session,"Group_ML2", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        }
        updateSelectInput(session, "Group_ML2", choices = Col_NamesML, selected = Col_NamesML[[1]])
        
        
        
        
        # Separates the profiles data by lipid mediators types:
        
        
        name_List <- unique(unname(unlist(lm_profile[1, -1])))
        
        name_List<-c("All_LM", name_List)
        
        #update the names of the models available in the datafile 
        updateSelectInput(session, "Build_RF_Mod_Num", choices = name_List)
        updateSelectInput(session, "Build_RF_Mod_Num1", choices = name_List)
        updateSelectInput(session, "Build_XGB_Mod_Num", choices = name_List)
        updateSelectInput(session, "Build_XGB_Mod_Num1", choices = name_List)
        updateSelectInput(session, "Build_SVM_Mod_Num", choices = name_List)
        updateSelectInput(session, "Build_LA_Mod_Num", choices = name_List)
        updateSelectInput(session, "Build_BC_Mod_Num", choices = name_List)
        
      },
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Group_ML2", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
  })
  
  observeEvent(input$Build_ML, {
    
    req(input$LM_ML_File2)
    
    tryCatch(
      {
        # Call the lipid mediator profiling file variant:
        inFile <- input$LM_ML_File2
        
        print(inFile)
        
        if(is.null(inFile))
          
          return(NULL)
        
        
        set.seed(415) # To get same results even with the random part.
        options(digits = 3) # To get only part of the decimals. 
        
        # Open the file:
        
        
        lm_profile = read.table(inFile$datapath, 
                                sep=input$sep_LM_ML2,
                                header = TRUE,
                                row.names = 1,
                                stringsAsFactors = FALSE)
        
        if (input$samp_ML2 == "col_ML2"){
          
          lm_profile <- data.frame(t(lm_profile))
          
          #remove spaces from column names 
          colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
          
          for (x in 1:ncol(lm_profile)){
            if (substr(colnames(lm_profile)[x],1,1) %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
              colnames(lm_profile)[x]= paste("X", colnames(lm_profile)[x], sep = '')
            }
            else{colnames(lm_profile)[x] = colnames(lm_profile)[x]}
          }
          
        }
        
        
        Met_Rename<-c()
        #remove spaces from column names 
        colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
        
        for (x in 1:ncol(lm_profile)){
          if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
            Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
          }else{
            Met_Rename[x] <-colnames(lm_profile)[x]
          }
        }
        colnames(lm_profile)<-Met_Rename
        
        #list of the unique substrates for the plots 
        name_List <- unique(unname(unlist(lm_profile[1, -1])))
        name_List<-c("All_LM", name_List)
        
        # Separates the profiles data by lipid mediators types:
        
        substrates <- unique(unname(unlist(lm_profile[1, -1])))
        
        # Creates data frames for each subset:
        
        dataframes_list <- list("ALL LM" = lm_profile[, -1])
        
        for (i in 1:length(substrates)) {
          substrate <- lm_profile[ ,lm_profile[1, ]== substrates[[i]]]
          assign(substrates[[i]], substrate)
          
        }
        
        if (exists('DHA') == TRUE) {dataframes_list[["DHA"]] <- DHA}
        if (exists('n3DPA') == TRUE) {dataframes_list[["n3DPA"]] <- n3DPA}
        if (exists('EPA') == TRUE) {dataframes_list[["EPA"]] <- EPA}
        if (exists('AA') == TRUE) {dataframes_list[["AA"]] <- AA}
        
        
        # Creates the data frame that is going to collect all the info regarding the models:
        
        final_table <- data.frame(machine_learning = "del",
                                  groups = "del",
                                  percentage_accuracy = 1,
                                  sensitivity = 1,
                                  specificity = 1,
                                  TP_per = 1,
                                  FP_per = 1,
                                  TN_per = 1,
                                  FN_per = 1,
                                  stringsAsFactors = FALSE)
        #dataframe with top 5 models for XGBoost
        XGB_Model_Table <- data.frame(model = "del",
                                      Accuracy = 1,
                                      specificity = 1,
                                      sensitivity = 1,
                                      TP_per = 1, 
                                      FP_per = 1,
                                      TN_per = 1,
                                      FN_per = 1, 
                                      nrounds = 1,
                                      eta= 1,
                                      max_depth= 1,
                                      gamma= 1,
                                      min_child_weight= 1)
        
        
        #progress bar
        #withProgress(message = 'Building Models', style = style, value = 0.1, {
        # Sys.sleep(0.25)})
        
        
        withProgress(message = 'Building Models', value = 0, {
          
          for (j in 1:length(dataframes_list)) {
            
            lm_profiles <- dataframes_list[[j]]
            
            Fatty_Acid<-lm_profiles[1,]
            
            # Save all the values as numeric:
            lm_profile_number <- sapply(lm_profiles[-1, ], function(x) as.numeric(x))
            row.names(lm_profile_number) <- row.names(lm_profiles[-1, ])
            
            # Scale the data because is better when you are running Machine Learning models:
            lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))
            
            # If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors
            # replace the NA for zeros. 
            lm_profiles_scale[is.na(lm_profiles_scale)] <- 0
            
            #--------------------
            ### change 0 to 1/5 of lowest value 
            
            
            lm_profile_zero<-lm_profiles_scale[-1,]
            
            #get vector the same length as number of columns 
            cols<- 1:(ncol(lm_profile_zero))
            #replace zeros for each column to 1/5 the smallest value for each column 
            
            lm_profile_zero[cols] <- lapply(lm_profile_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
            
            #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
            lm_profile_zero<-as.matrix(lm_profile_zero)
            lm_profile_zero[is.infinite(lm_profile_zero)] <- NA
            min_index <- which.min(lm_profile_zero)
            zero_replace <- (lm_profile_zero[min_index]/5)
            lm_profile_zero <- as.data.frame(lm_profile_zero)
            lm_profile_zero[is.na(lm_profile_zero)] <- zero_replace  
            
            #add first row to dataframe 
            lm_profiles_scale<-rbind(lm_profiles_scale[1,] ,lm_profile_zero)
            
            #exclude columns with the same value for all rows
            #lm_profile_scale<-Filter(var, lm_profile_scale[-1,])
            
            #---------------------
            # Add the classification variable to the data frame (Responder and non responder):
            lm_profiles_scale$responses <- factor(lm_profile[-1, ]$groups)
            
            # Add the classification variable to the data frame (Responder and non responder):
            ResponseName<-as.character(input$Group_ML2)
            lm_profiles_scale$responses <- factor(lm_profile[-1, ResponseName])
            
            # Make sure that column names do not represent a problem to randomForest making them a valid name to R.
            names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))
            
            #---> MACHINE LEARNING (randomForest R):
            if (input$RF_Build == TRUE){
              # In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model. 
              oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.
              for (mtry in 1:(ncol(lm_profiles_scale) - 1)) { 
                rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry, 
                                                      importance = TRUE, ntree = input$RF_ntrees)
                oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[input$RF_ntrees])*100)
              }
              # Define the best mtry according to the best prediction value. 
              final_mtry <- which.max(oob_error)
              
              # Run the model again with the right mtry value. 
              rf_lm_profiles_final <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry, 
                                                   importance = TRUE, ntree = input$RF_ntrees)
              # Get the confusion matrix of the model, sensitivity and specificity: 
              confusion_lm_profiles <- as.data.frame(rf_lm_profiles_final$confusion)
              confusion_lm_profiles[is.na(confusion_lm_profiles)] <- 0
              
              # Calculates sensitivity, specificity and AUC.
              sensitivity_lm_profiles <- confusion_lm_profiles[2, 2]/(confusion_lm_profiles[2, 2] + confusion_lm_profiles[2, 1])
              specificity_lm_profiles <- confusion_lm_profiles[1, 1]/(confusion_lm_profiles[1, 1] + confusion_lm_profiles[1, 2])
              
              # Final table for random forest:
              no_oob_error_table <- data.frame(machine_learning = "RF",
                                               groups = names(dataframes_list)[j],
                                               percentage_accuracy = 100 - ((rf_lm_profiles_final$err.rate[input$RF_ntrees])*100),
                                               sensitivity = sensitivity_lm_profiles,
                                               specificity = specificity_lm_profiles,
                                               TP_per = (1 - confusion_lm_profiles[2, 3])*100,
                                               FP_per = confusion_lm_profiles[1, 3]*100,
                                               TN_per = (1 -confusion_lm_profiles[1, 3])*100,
                                               FN_per = confusion_lm_profiles[2, 3]*100,
                                               stringsAsFactors = FALSE)
              final_table <- rbind(final_table, no_oob_error_table)
              # Number of trees plot:
              tree_table <- data.frame(Ensemble = c(1:input$RF_ntrees),
                                       OBB = rf_lm_profiles_final$err.rate[, 1],
                                       AvgAcc = 100 - ((rf_lm_profiles_final$err.rate[, 1])*100))
              ntree_plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) +
                geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
                ggtitle(names(dataframes_list)[j]) +
                scale_x_continuous(name = "Trees") +
                scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
                theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                      axis.title.y = element_text(size = 25, colour = "black"),
                      axis.title.x = element_text(size = 25, colour = "black"),
                      axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                      axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep. 
                      legend.position = ("none"),
                      panel.background = element_rect(fill = "white"),
                      axis.ticks.length = unit(0.4, "cm"))
              assign(paste("ntree_", names(dataframes_list)[j], sep = ""), ntree_plot)
              
              
              #get importance dataframe 
              rf_Importance<-as.data.frame(importance(rf_lm_profiles_final, type = 1))
              rf_Importance$lipid_mediators <- rownames(rf_Importance)
              
              Fatty_Acid <- 1:length(rf_Importance$lipid_mediators)
              cbind(rf_Importance, Fatty_Acid)
              
              
              # #fatty acid points 
              for (i in 1:length(rf_Importance$lipid_mediators)){
                rf_Importance$Fatty_Acid[i]<-SPMs_FA[rf_Importance$lipid_mediators[i]]
              }
              
              #change the names to be in the proper format 
              rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)] <- paste(gsub("^X", "'", rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)]), "'", sep = "")
              rf_Importance$lipid_mediators <- gsub("4$", "[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
              rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
              rf_Importance$lipid_mediators<- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required 
              rf_Importance$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", rf_Importance$lipid_mediators) # Transform n-3 DPA as a subscript: 
              
              # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
              # (this with the purpose of get the LM in decreasing order in the figure):
              
              rf_Importance <- rf_Importance[order(rf_Importance$MeanDecreaseAccuracy), ]
              
              rf_Importance$lipid_mediators <- factor( rf_Importance$lipid_mediators, 
                                                       levels = rf_Importance$lipid_mediators[
                                                         order(rf_Importance$MeanDecreaseAccuracy)])
              
              lipids <- as.character(rf_Importance$lipid_mediators)
              #add colors to SPM names
              lm_classes <- as.factor(rf_Importance$Fatty_Acid)
              names(lm_classes) <- rf_Importance$lipid_mediators
              
              dha_index<-which((lm_classes=="DHA")== TRUE)
              n_three_index<-which((lm_classes=="n3DPA")== TRUE)
              epa_index<-which((lm_classes=="EPA")== TRUE)
              aa_index<-which((lm_classes=="AA")== TRUE)
              
              lm_colors <- NULL
              lm_colors[dha_index] <- "blue"
              lm_colors[n_three_index] <- "brown"
              lm_colors[epa_index] <- "darkgoldenrod1"
              lm_colors[aa_index] <- "darkslategray"
              
              #rf_Importance1<- subset(rf_Importance, MeanDecreaseAccuracy >1)
              rf_VIP_plot <- ggplot(data = rf_Importance, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = Fatty_Acid)) + geom_point(size = 3) +
                scale_y_continuous(name = "Mean Decrease Accuracy") +
                labs(x = "Lipid Mediators", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
                scale_x_discrete(labels = parse(text = lipids)) + 
                scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                     "EPA" ="darkgoldenrod1", "AA" = "darkslategray"
                )) +
                coord_flip() + 
                theme(axis.title = element_text(size = 20),
                      axis.text.x  = element_text(size = 15, hjust = 0.5),
                      axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                      legend.position = "top",
                      aspect.ratio = 2/1,
                      legend.title = element_text(size = 20),
                      legend.text  = element_text(size = 15),
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
              
              
              assign(paste("rf_VIP_", names(dataframes_list)[j], sep = ""), rf_VIP_plot)
              assign(paste("rf_Mod_", names(dataframes_list)[j], sep = ""), rf_lm_profiles_final)
              
              
            }
            
            #---> MACHINE LEARNING (EXTREME GRADIENT BOOSTING)
            if (input$XGB_Build == TRUE){
              df<-lm_profiles_scale
              #need to convert responses to integer of 0 and 1
              
              df$responses<-as.factor(df$responses)
              responses<-df$responses
              label<-as.integer(df$responses) -1
              
              #remove labels from dataframe
              df$responses<-NULL
              df<-sapply(df, function(x) as.numeric(x))
              
              #
              #split it into training and testing dataset with 75/20
              n = nrow(df) #get number of rows for dataframe
              train.index = sample(n,floor(0.7*n)) #randomize rows to get dataset
              train.data = as.matrix(df[train.index,])
              train.label = label[train.index]
              test.data = as.matrix(df[-train.index,])
              test.label = label[-train.index]
              
              #transform dataset to xgb matrix
              xgb.train = xgb.DMatrix(data=train.data,label=train.label)
              xgb.test = xgb.DMatrix(data=test.data, label=test.label)
              
              #create grid with all the options for xgboost
              searchGridSubCol <- expand.grid(max_depth = c(1, 2, 3, 4, 6),
                                              gamma_val = c(0.5, 1, 1.5, 2, 2, 5), 
                                              eta = c(0.01, 0.02, 0.03, 0.04, 0.05),
                                              min_child = c(1, 2, 3, 4, 5)
              )
              
              #run through each combination of the grid
              system.time(
                ErrorsHyperparameters <- apply(searchGridSubCol, 1, function(parameterList){
                  
                  #Extract Parameters to test
                  currentDepth <- parameterList[["max_depth"]]
                  currentEta <- parameterList[["eta"]]
                  currentGamma <- parameterList[["gamma_val"]]
                  currentMinChild <- parameterList[["min_child"]]
                  
                  #run selected parameter through the model
                  xgboostModelCV <- xgb.cv(data =  xgb.train, nrounds = input$XGB_nrounds, nfold = 5, showsd = TRUE,
                                           metrics = "error", verbose = TRUE, "eval_metric" = "error",
                                           "objective" = "binary:logistic", "max.depth" = currentDepth, "eta" = currentEta,
                                           "subsample" = 1, "colsample_bytree" = 1, print_every_n = 10, booster = "gbtree",
                                           early_stopping_rounds = 10, "gamma" = currentGamma, "min_child_weight" = currentMinChild)
                  
                  #have error evaluation score as dataframe
                  xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
                  test_error <- tail(xvalidationScores$test_error_mean, 1)
                  train_error <- tail(xvalidationScores$train_error_mean,1)
                  output <- return(c(test_error, train_error, currentDepth, currentEta, currentGamma, currentMinChild))})) #add to final table called output
              
              
              output_Error <- as.data.frame(t(ErrorsHyperparameters))
              
              varnames <- c("TestError", "TrainError","Depth", "eta", "gamma", "MinChild")
              names(output_Error) <- varnames #change the variable names
              top_output<-output_Error[order(output_Error$TestError, output_Error$TrainError), ] #order it by the lowest error rate
              
              #empty table to put the top 10 model information into
              Final_Acc <- data.frame(model = "test",
                                      Accuracy = 0,
                                      specificity = 0,
                                      sensitivity = 0,
                                      TP_per = 0,
                                      FP_per = 0,
                                      TN_per = 0,
                                      FN_per = 0,
                                      nrounds = 0,
                                      eta= 0,
                                      max_depth= 0,
                                      gamma= 0,
                                      min_child_weight=0
              )
              
              top<-1:50
              
              #loop to build models for the parameters with the lowest error rate
              for (i in top){
                #parameter list with the lowest output
                params <- list(booster = "gbtree", objective = "binary:logistic", eta= top_output[i,4], gamma= top_output[i,5],
                               max_depth= top_output[i,3], min_child_weight=top_output[i,6], subsample=1,
                               colsample_bytree= 1)
                #Build model based on the training data
                xgb<-xgb.train(data = xgb.train, params = params, nrounds = input$XGB_nrounds, eval.metric = "error", early.stop.rounds=10,  silent = 0)
                
                #run model on test data = see how good it is at predicting
                pred<-predict(xgb, test.data)
                
                #obtain confusion matrix
                pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1
                pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
                pred<-as.numeric(pred) #makes it numeric
                test.label_num <-as.numeric(test.label) #makes it numeric
                pred_y <- cut(pred, 2, label = c("Non-responder", "Responder")) #converts to factor of responder or non responder
                test_y <-cut(test.label_num, 2, label = c("Non-responder", "Responder")) #convers to factor of responder or non responder
                ConMatrix <- confusionMatrix(reference = test_y, data = pred_y) #build confusion matrix
                TP = as.numeric(ConMatrix$table[2,2])
                FP = as.numeric(ConMatrix$table[2,1])
                TN = as.numeric(ConMatrix$table[1,1])
                FN = as.numeric(ConMatrix$table[1,2])
                
                #convert build table with accuracy for each model
                Table <- data.frame(model = names(dataframes_list)[j],
                                    Accuracy = as.numeric(ConMatrix$overall[1])*100,
                                    specificity = as.numeric(ConMatrix$byClass[2]),
                                    sensitivity = as.numeric(ConMatrix$byClass[1]),
                                    TP_per = (TP/(TP+FN)*100),
                                    FP_per = (FP/(FP+TN)*100),
                                    TN_per = (TN/(FP+TN)*100),
                                    FN_per = (FN/(TP+FN)*100),
                                    nrounds = input$XGB_nrounds,
                                    eta= top_output[i,4],
                                    max_depth= top_output[i,3],
                                    gamma= top_output[i,5],
                                    min_child_weight=top_output[i,6]
                )
                Final_Acc<-rbind(Final_Acc, Table)
                
              }
              
              #remove the first row of random values
              Final_Acc<-Final_Acc[-1,]
              
              #order the top 10 models based on highest accuracy, sensitivity, and specificity
              Top_Acc<-Final_Acc[with(Final_Acc, order(-Accuracy, -sensitivity, -specificity)), ]
              
              Top_Acc<-Top_Acc[1:5,]
              
              Top_Acc1 <- Top_Acc
              top5<-1:5
              for (i in top5){
                Top_Acc1[i,1] <-paste(names(dataframes_list)[j],"_Model", i,  sep = "")
              }
              
              #selecting the model number based on the users input 
              if ( names(dataframes_list)[j] == "ALL LM"){
                Model_Num<-input$XGB_ALL_LM_Model
              }
              if ( names(dataframes_list)[j] == "DHA"){
                Model_Num<-input$XGB_DHA_Model
              }
              if ( names(dataframes_list)[j] == "n3DPA"){
                Model_Num<-input$XGB_n3DPA_Model
              }
              if ( names(dataframes_list)[j] == "EPA"){
                Model_Num<-input$XGB_EPA_Model
              }
              if ( names(dataframes_list)[j] == "AA"){
                Model_Num<-input$XGB_AA_Model
              }
              
              # XGBoost Table:
              XGBoost_table <- data.frame(machine_learning = "XGB",
                                          groups = names(dataframes_list)[j],
                                          percentage_accuracy =Top_Acc[Model_Num,2],
                                          sensitivity = Top_Acc[Model_Num,3],
                                          specificity = Top_Acc[Model_Num,4],
                                          TP_per = Top_Acc[Model_Num,5],
                                          FP_per = Top_Acc[Model_Num,6],
                                          TN_per = Top_Acc[Model_Num,7],
                                          FN_per = Top_Acc[Model_Num,8],
                                          stringsAsFactors = FALSE)
              
              
              
              final_table <- rbind(final_table, XGBoost_table)
              XGB_Model_Table <- rbind(XGB_Model_Table, Top_Acc1)
              
              #preparing nrounds plots for top 5 models
              iterations<-1:input$XGB_nrounds
              top_table<- data.frame(iterations)
              
              for (i in top5){
                
                #paramters for the i model
                params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[i,10], gamma= Top_Acc[i,12],
                               max_depth= Top_Acc[i,11], min_child_weight= Top_Acc[i,13], subsample= 1)
                
                #run the cv to get the error rate for each nround
                xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = input$XGB_nrounds, nfold = 5, showsd = T,
                                 stratified = T, early_stop_round = 20, maximize = F, metrics = "error", verbose = FALSE)
                
                Table1<-data.frame( test = as.numeric(xgbcv$evaluation_log$test_error_mean))
                colnames(Table1)[1] <-paste(names(dataframes_list)[j],"_Model_", i, sep = "")
                #colnames(Table1)[2] <-paste(names(dataframes_list)[j],"Test_model_", i, sep = "")
                top_table<-cbind(top_table, Table1) #add to table
                
              }
              
              #plot using plotly
              nround_plot<-plot_ly(data = top_table, x = ~iterations)
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[2])), name = colnames(top_table)[2], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[3])), name = colnames(top_table)[3], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[4])), name = colnames(top_table)[4], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[5])), name = colnames(top_table)[5], mode = 'lines')
              nround_plot <-  nround_plot %>% add_trace(y = as.formula(paste0('~', top_table[6])), name = colnames(top_table)[6], mode = 'lines')
              nround_plot <- nround_plot %>% layout(title = paste(names(dataframes_list)[j], "Plot", sep =" "), xaxis = list(title = 'Number of Rounds'),
                                                    yaxis = list(title = 'Test Error Rate'))
              
              #makes download quality of plot better
              nround_plot <- nround_plot%>% config(toImageButtonOptions = list(format = "jpeg", width = 1500, height = 750))
              
              assign(paste("nround_", names(dataframes_list)[j], sep = ""), nround_plot)
              
              #importance plot
              #build model for the top model 
              params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[Model_Num,10], gamma= Top_Acc[Model_Num,12], 
                             max_depth= Top_Acc[Model_Num,11], min_child_weight= Top_Acc[Model_Num, 13], subsample= 1, colsample_bytree= 1)
              xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = input$XGB_nrounds, eval.metric = "error", early.stop.rounds=10,  silent = 0)
              
              
              XGB_Importance = xgb.importance(model = xgb_Mod1) #important matrix
              
              Fatty_Acid <- 1:length(XGB_Importance$Feature)
              cbind(XGB_Importance, Fatty_Acid)
              
              #fatty acid points 
              for (i in 1:length(XGB_Importance$Feature)){
                XGB_Importance$Fatty_Acid[i]<-SPMs_FA[XGB_Importance$Feature[i]]
              }
              
              # #change the names to be in the proper format
              XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)] <- paste(gsub("^X", "'", XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)]), "'", sep = "")
              XGB_Importance$Feature <- gsub("4$", "[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 4
              XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature ) #B4 converted to subscript 4
              XGB_Importance$Feature  <- gsub("\\.", "-", XGB_Importance$Feature ) # Replace "." with "-" when required
              XGB_Importance$Feature <- gsub("n.3.DPA", "[n-3~DPA]", XGB_Importance$Feature ) # Transform n-3 DPA as a subscript:
              
              
              # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
              # (this with the purpose of get the LM in decreasing order in the figure):
              
              XGB_Importance <- XGB_Importance[order(XGB_Importance$Gain), ]
              
              XGB_Importance$Feature<- factor(XGB_Importance$Feature, 
                                              levels = XGB_Importance$Feature[
                                                order(XGB_Importance$Gain)])
              lipids <- as.character(XGB_Importance$Feature )
              
              #add colors to x axis labels
              lm_classes <- as.factor(XGB_Importance$Fatty_Acid)
              names(lm_classes) <- XGB_Importance$Feature
              
              dha_index<-which((lm_classes=="DHA")== TRUE)
              n_three_index<-which((lm_classes=="n3DPA")== TRUE)
              epa_index<-which((lm_classes=="EPA")== TRUE)
              aa_index<-which((lm_classes=="AA")== TRUE)
              
              lm_colors <- NULL
              lm_colors[dha_index] <- "blue"
              lm_colors[n_three_index] <- "brown"
              lm_colors[epa_index] <- "darkgoldenrod1"
              lm_colors[aa_index] <- "darkslategray"
              
              XGB_VIP_Plot <- ggplot(data = XGB_Importance, mapping = aes(x = Feature, y = Gain, color = Fatty_Acid)) + geom_point(size = 3) +
                scale_y_continuous(name = "Gain") +
                labs(x = "Features", title = paste(names(dataframes_list)[j], "Plot", sep =" ")) +
                scale_x_discrete(labels = parse(text = lipids)) +
                scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                     "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
                coord_flip() +
                theme(axis.title = element_text(size = 20),
                      axis.text.x  = element_text(size = 15, hjust = 0.5),
                      axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                      legend.position = "top",
                      panel.background = element_rect(fill = "white"),
                      panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
              
              #save xgboost VIP plot
              assign(paste("XGB_VIP_", names(dataframes_list)[j], sep = ""), XGB_VIP_Plot)
              #Save Xgboost model 
              assign(paste("xgb_Mod_", names(dataframes_list)[j], sep = ""), xgb_Mod1)
            }
            #---> MACHINE LEARNING (Classyfire R):
            
            # Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
            # and creates a novelty detection for the creation of the model.
            
            # The idea is to create several models and see which one fits the best. The models will be based on the whole
            # lipid profiles and the different groups based on substrates.
            
            # "cfBuild" to create the SVM:
            # Clasyfire requieres matrix:
            if (input$SVM_Build == TRUE){
              #convert dataframe to matrix
              lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])
              
              #building SVM model
              support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses,
                                                  bootNum = input$SVM_BootNum, ensNum = input$SVM_Ensemble, cpus = 4)
              #obtaining confusion matrix
              conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale))
              
              # SVM table:
              Support_vector_table <- data.frame(machine_learning = "SVM",
                                                 groups = names(dataframes_list)[j],
                                                 percentage_accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                                 sensitivity = conf_matrix[4, 3]/100,
                                                 specificity = conf_matrix[1, 3]/100,
                                                 TP_per = conf_matrix[4, 3],
                                                 FP_per = conf_matrix[2, 3],
                                                 TN_per = conf_matrix[1, 3],
                                                 FN_per = conf_matrix[3, 3],
                                                 stringsAsFactors = FALSE)
              final_table <- rbind(final_table, Support_vector_table)
              
              # Ensemble Plot:
              ensAcc   <- getAcc(support_lmprofiles_scale)$Test
              meanVal  <- ensAcc[1]
              for (i in 2:length(ensAcc)) {
                meanVal <- c(meanVal, mean(ensAcc[1:i]))
              }
              ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc),
                                          AvgAcc = meanVal)
              ensemble_plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) +
                geom_point(aes(colour = AvgAcc), size = 5) +
                geom_line(linetype = "dotted", size = 1) +
                ggtitle(names(dataframes_list)[j]) +
                scale_x_continuous(name = "Ensemble interaction") +
                scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
                theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                      axis.title.y = element_text(size = 25, colour = "black"),
                      axis.title.x = element_text(size = 25, colour = "black"),
                      axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                      axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep.
                      legend.position = ("none"),
                      panel.background = element_rect(fill = "white"),
                      axis.ticks.length = unit(0.4, "cm"))
              assign(paste("ensemble_", names(dataframes_list)[j], sep = ""), ensemble_plot)
              
              #Save SVM model 
              assign(paste("svm_Mod_", names(dataframes_list)[j], sep = ""), support_lmprofiles_scale)
            }
            
            #---> ELASTIC NET REGRESSION (caret R):
            if (input$LA_Build == TRUE){
              # Get the explanatory variables as a matrix:
              explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
              
              # LASSO Analysis:
              model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet",
                                 trControl = trainControl("boot", number = input$LA_BootNum))
              # Get the confusion matrix of the model:
              conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
              
              # Final Elastic net model table:
              net_table <- data.frame(machine_learning = "GLMNET",
                                      groups = names(dataframes_list)[j],
                                      percentage_accuracy = max(model_net$results$Accuracy)*100,
                                      sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                                      specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                                      TP_per = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                      FP_per = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                      TN_per = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                      FN_per = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                      stringsAsFactors = FALSE)
              final_table <- rbind(final_table, net_table)
              
              # Parameter tuning figure:
              scaleFUN <- function(x) sprintf("%.2f", x)
              
              boot_net_plot <- ggplot(model_net, highlight = TRUE) +
                scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
                scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
                ggtitle(names(dataframes_list)[j]) +
                scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
                scale_shape_manual(values=c(16, 16, 16)) +
                labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
                guides(shape = FALSE) +
                theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                      axis.title.y = element_text(size = 25, colour = "black"),
                      axis.title.x = element_text(size = 25, colour = "black"),
                      axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                      axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
                      axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep.
                      legend.title = element_text(size = 25),
                      legend.text  = element_text(size = 20),
                      legend.key = element_rect(fill = NA),
                      legend.key.size = unit(1.3, "cm"),
                      panel.background = element_rect(fill = "white"),
                      axis.ticks.length = unit(0.4, "cm"))
              assign(paste("net_", names(dataframes_list)[j], sep = ""), boot_net_plot)
              
              #rename and save elastic net model 
              assign(paste("la_Mod_", names(dataframes_list)[j], sep = ""), model_net)
            }
            
            #---> BAYESIAN MODEL (Caret R):
            
            if (input$BC_Build == TRUE){
              explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
              
              bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm",
                                trControl = trainControl("boot", number = input$BC_BootNum))
              conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)
              
              # Final Bayesian model table:
              bay_table <- data.frame(machine_learning = "BAYES",
                                      groups = names(dataframes_list)[j],
                                      percentage_accuracy = max(bayesian$results$Accuracy)*100,
                                      sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                                      specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                                      TP_per = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                      FP_per = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                      TN_per = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                      FN_per = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                      stringsAsFactors = FALSE)
              final_table <- rbind(final_table, bay_table)
              
              #rename and save bayesian model 
              assign(paste("bc_Mod_", names(dataframes_list)[j], sep = ""), bayesian)
            }
            
            incProgress(1/length(dataframes_list), detail = paste(names(dataframes_list)[j]))
          }
          
          final_table <- final_table[-1, ]
          table <- data.frame(Model = factor(final_table$groups, levels = unique(final_table$groups)),
                              methodology = final_table$machine_learning,
                              accuracy = round(final_table$percentage_accuracy, 0))
          accuracy <- ggplot(data = table, aes(x = Model, y = accuracy, fill = methodology)) +
            geom_bar(stat = "identity",  position = position_dodge2(preserve = 'single'), colour = "black", size = 2) +
            scale_fill_manual(values=c("dodgerblue2","firebrick2",'goldenrod1', 'lightslategray', "purple")) + 
            geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
                      aes(label = paste(accuracy, "%", sep = "")), size = 8) +
            scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20),
                               expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
            coord_cartesian(ylim = c(1, 120)) +
            theme(axis.title = element_text(size = 40),
                  axis.title.x = element_blank(),
                  axis.text.x  =  element_text(size = 40, hjust = 0.5, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 40, hjust = 1, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.0), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.0), # Color and thickness of every axis sep. 
                  panel.background = element_rect(fill = "white"),
                  legend.title = element_blank(),
                  legend.position = "top",
                  legend.key.size = unit(1.3, "cm"), 
                  legend.text  = element_text(size = 25),
                  legend.spacing.x = unit(1, "cm"),
                  axis.ticks.length = unit(0.4, "cm"))
          
          #Put figures together: 
          if (input$RF_Build == TRUE){
            rf_plot_list <- list("ALL LM" = `ntree_ALL LM`)
            rf_VIPplot_list <- list("ALL LM" = `rf_VIP_ALL LM`)
            rf_Mod_list<-list("ALL LM" = `rf_Mod_ALL LM`)
          }
          if (input$XGB_Build == TRUE){
            XGB_VIPplot_list <- list("ALL LM" = `XGB_VIP_ALL LM`)
            xgb_Mod_list<-list("ALL LM" = `xgb_Mod_ALL LM`)}
          if (input$SVM_Build == TRUE){
            svm_plot_list <- list("ALL LM" = `ensemble_ALL LM`)
            svm_Mod_list<-list("ALL LM" = `svm_Mod_ALL LM`)}
          if (input$LA_Build == TRUE){
            net_plot_list <- list("ALL LM" = `net_ALL LM`)
            la_Mod_list<-list("ALL LM" = `la_Mod_ALL LM`)}
          if (input$BC_Build == TRUE){bc_Mod_list<-list("ALL LM" = `bc_Mod_ALL LM`)}
          
          
          column <- 2
          if (exists('DHA') == TRUE) {
            if (input$RF_Build == TRUE){
              rf_plot_list[["DHA"]] <- ntree_DHA
              rf_VIPplot_list[["DHA"]] <- rf_VIP_DHA
              rf_Mod_list[["DHA"]] <- rf_Mod_DHA
            }
            if (input$XGB_Build == TRUE){
              XGB_VIPplot_list[["DHA"]] <- XGB_VIP_DHA
              xgb_Mod_list[["DHA"]] <- xgb_Mod_DHA}
            if (input$SVM_Build == TRUE){
              svm_plot_list[["DHA"]] <- ensemble_DHA
              svm_Mod_list[["DHA"]] <- svm_Mod_DHA}
            if (input$LA_Build == TRUE){
              net_plot_list[["DHA"]] <- net_DHA
              la_Mod_list[["DHA"]] <- la_Mod_DHA}
            if (input$BC_Build == TRUE){bc_Mod_list[["DHA"]] <- bc_Mod_DHA}
            
            column <- column + 1}
          
          if (exists('n3DPA') == TRUE) {
            if (input$RF_Build == TRUE){
              rf_plot_list[["n3DPA"]] <- ntree_n3DPA
              rf_VIPplot_list[["n3DPA"]] <- rf_VIP_n3DPA
              rf_Mod_list[["n3DPA"]] <- rf_Mod_n3DPA
            }
            if (input$XGB_Build == TRUE){
              XGB_VIPplot_list[["n3DPA"]] <- XGB_VIP_n3DPA
              xgb_Mod_list[["n3DPA"]] <- xgb_Mod_n3DPA}
            if (input$SVM_Build == TRUE){
              svm_plot_list[["n3DPA"]] <- ensemble_n3DPA
              svm_Mod_list[["n3DPA"]] <- svm_Mod_n3DPA}
            if (input$LA_Build == TRUE){
              net_plot_list[["n3DPA"]] <- net_n3DPA
              la_Mod_list[["n3DPA"]] <- la_Mod_n3DPA}
            if (input$BC_Build == TRUE){bc_Mod_list[["n3DPA"]] <- bc_Mod_n3DPA}
            
            column <- column + 1}
          
          
          if (exists('EPA') == TRUE) {
            if (input$RF_Build == TRUE){
              rf_plot_list[["EPA"]] <- ntree_EPA
              rf_VIPplot_list[["EPA"]] <- rf_VIP_EPA
              rf_Mod_list[["EPA"]] <- rf_Mod_EPA
            }
            if (input$XGB_Build == TRUE){
              XGB_VIPplot_list[["AA"]] <- XGB_VIP_AA
              xgb_Mod_list[["EPA"]] <- xgb_Mod_EPA}
            if (input$SVM_Build == TRUE){
              svm_plot_list[["EPA"]] <- ensemble_EPA
              svm_Mod_list[["EPA"]] <- svm_Mod_EPA}
            if (input$LA_Build == TRUE){
              net_plot_list[["EPA"]] <- net_EPA
              la_Mod_list[["EPA"]] <- la_Mod_EPA}
            if (input$BC_Build == TRUE){bc_Mod_list[["EPA"]] <- bc_Mod_EPA}
            column <- column + 1}
          
          if (exists('AA') == TRUE) {
            if (input$RF_Build == TRUE){
              rf_plot_list[["AA"]] <- ntree_AA
              rf_VIPplot_list[["AA"]] <- rf_VIP_AA
              rf_Mod_list[["AA"]] <- rf_Mod_AA}
            if (input$XGB_Build == TRUE){
              XGB_VIPplot_list[["AA"]] <- XGB_VIP_AA
              xgb_Mod_list[["AA"]] <- xgb_Mod_AA}
            if (input$SVM_Build == TRUE){
              svm_plot_list[["AA"]] <- ensemble_AA
              svm_Mod_list[["AA"]] <- svm_Mod_AA}
            if (input$LA_Build == TRUE){
              net_plot_list[["AA"]] <- net_AA
              la_Mod_list[["AA"]] <- la_Mod_AA}
            if (input$BC_Build == TRUE){bc_Mod_list[["AA"]] <- bc_Mod_AA}
            
            column <- column + 1}
          if (column >= 6) {column <- 4}
          
          #if (RF_ML){grid.arrange(grobs = rf_plot_list, ncol = column/2)}
          else if (input$SVM_Build == TRUE){grid.arrange(grobs = svm_plot_list, ncol = column/2)}
          else if (input$LA_Build == TRUE){grid.arrange(grobs = net_plot_list, ncol = column/2)}
          
          # Final table  outputs:
          final_table$`% Accuracy Score` <- round(final_table$percentage_accuracy, 0)
          final_table$Sensitivity <- round(final_table$sensitivity, 2)
          final_table$Specificity <- round(final_table$specificity, 2)
          final_table$TP <- round(final_table$TP_per, 0)
          final_table$FP <- round(final_table$FP_per, 0)
          final_table$TN <- round(final_table$TN_per, 0)
          final_table$FN <- round(final_table$FN_per, 0)
          
          final_table <- final_table[, c(1, 2, 10:16)]
          colnames(final_table)[1] <- "Machine Learning Methodology"
          colnames(final_table)[2] <- "Model"
          inal_table<-as.data.frame(final_table)
          
          #outputs for accuracy tab
          output$Build_AccPlot_Title <- renderUI({req(input$Build_ML); h2("% Accuracy Score Figure for different ML models", align = "center") })
          output$Build_AccTable_Title <- renderUI({req(input$Build_ML); h2("Model Table Summary", align = "center") })
          
          output$Build_Accuracy_ML <- renderPlot({return(accuracy)}) #accuracy plot output
          output$downloadBuild_Acc_ML_Plot <- downloadHandler(filename = function(){paste("Accuracy_Plot",'.png',sep='')},
                                                              content = function(file){
                                                                ggsave(file,plot= accuracy, width = 15, height = 10)})
          
          output$Build_ML_Table <- renderDataTable({return(final_table)}) #accuracy table output 
          output$downloadBuild_Acc_ML_Table <- downloadHandler(filename = function(){"Accuracy_Table.csv"}, 
                                                               content = function(fname){
                                                                 write.csv(final_table, fname)})
          
          #outputs if Random forest model is run 
          if (input$RF_Build == TRUE){
            
            #ui outputs for the page nemaes 
            output$Build_RF_Plot_Title <- renderUI({req(input$Build_ML); h2("Optimal Parameters RandomForest", align = "center") })
            output$Build_RF_VIPPlot_Title <- renderUI({req(input$Build_ML); h2("Importance of Variance Plot", align = "center") })
            
            #output for Random forest optimal parametesd plot and download 
            output$Build_RF_Plot <- renderPlot({return(grid.arrange(grobs = rf_plot_list, ncol = 1))})
            rf_plot_grid<-grid.arrange(grobs = rf_plot_list, ncol = 1)
            output$downloadBuild_RF_Plot <- downloadHandler(filename = function(){paste("Random_Forest_Plot",'.png',sep='')},
                                                            content = function(file){
                                                              ggsave(file, plot= rf_plot_grid, width = 12, height = 20)})
            
            #output and download for random forest VIP plot 
            output$Build_RF_VIPPlot <- renderPlot({return(grid.arrange(grobs = rf_VIPplot_list, ncol = 2))})
            output$downloadBuild_RF_VIPPlot <- downloadHandler(filename = function(){paste("Random_Forest_VIP_Plot",'.png',sep='')},
                                                               content = function(file){
                                                                 ggsave(file,plot= rf_VIPplot_list[[match(input$Build_RF_Mod_Num, name_List)]], width = 12, height = 12)})
            
            #downloading model information for paramters tab and importance tab
            output$downloadBuild_RF_Mod <- downloadHandler(filename = function(){paste("Random_Forest_Model_",input$Build_RF_Mod_Num,'.rds',sep='')},
                                                           content = function(fname){
                                                             saveRDS(rf_Mod_list[[match(input$Build_RF_Mod_Num, name_List)]], fname)})
            
            output$downloadBuild_RF_Mod1 <- downloadHandler(filename = function(){paste("Random_Forest_Model_",input$Build_RF_Mod_Num1,'.rds',sep='')},
                                                            content = function(fname){
                                                              saveRDS(rf_Mod_list[[match(input$Build_RF_Mod_Num1, name_List)]], compress = TRUE, fname)})
            
          }
          
          #outputs for xgbooost model 
          if (input$XGB_Build == TRUE){
            
            #ui outputs for plot titles 
            output$Build_XGB_Plot_Title <- renderUI({req(input$Build_ML); h2("Top 5 models for XGBoost test error rate", align = "center") })
            output$Build_XGB_VIPPlot_Title <- renderUI({req(input$Build_ML); h2("Importance of Variance Plot", align = "center") })
            
            #output for optimal parameters plot  and table 
            output$Build_XGB_Plot1 <- renderPlotly({return(`nround_ALL LM`)})
            if (exists('DHA') == TRUE){output$Build_XGB_Plot2 <- renderPlotly({return(nround_DHA)})}
            if (exists('n3DPA') == TRUE) {output$Build_XGB_Plot3 <- renderPlotly({return(nround_n3DPA)})}
            if (exists('EPA') == TRUE) {output$Build_XGB_Plot4 <- renderPlotly({return(nround_EPA)})}
            if (exists('AA') == TRUE) {output$Build_XGB_Plot5 <- renderPlotly({return(nround_AA)})}
            
            XGB_Model_Table<-XGB_Model_Table[-1,]
            XGB_Model_Table[,-1] <-round(XGB_Model_Table[,-1],2)
            output$XGB_Build_Models<-renderDataTable(XGB_Model_Table, rownames = FALSE,
                                                     options = list(paging = TRUE, pagelength = 20, scrollX = TRUE, scrollY = TRUE, autoWidth = TRUE))
            
            output$downloadBuild_XGB_Table <- downloadHandler(filename = function(){"XGB_Model_Table.csv"}, 
                                                              content = function(fname){
                                                                write.csv(XGB_Model_Table, fname)})
            
            #XGboost VIP plot outputs 
            output$Build_XGB_VIPPlot <- renderPlot({return(grid.arrange(grobs = XGB_VIPplot_list, ncol = 2))})
            output$downloadBuild_XGB_VIPPlot <- downloadHandler(filename = function(){paste("XGBoost_VIP_Plot",'.png',sep='')},
                                                                content = function(file){
                                                                  ggsave(file,plot= XGB_VIPplot_list[[match(input$Build_XGB_Mod_Num, name_List)]], width = 12, height = 12)})
            
            #downloading model information for paramters tab and importance tab for cgboost
            output$downloadBuild_XGB_Mod <- downloadHandler(filename = function(){paste("XGBoost_Model_",input$Build_XGB_Mod_Num,'.rds',sep='')},
                                                            content = function(fname){
                                                              saveRDS(xgb_Mod_list[[match(input$Build_XGB_Mod_Num, name_List)]], fname)})
            
            output$downloadBuild_XGB_Mod1 <- downloadHandler(filename = function(){paste("XGBoost_Model_",input$Build_XGB_Mod_Num1,'.rds',sep='')},
                                                             content = function(fname){
                                                               saveRDS(xgb_Mod_list[[match(input$Build_XGB_Mod_Num1, name_List)]], fname)})
            
            
          }
          
          
          #outputs for SVM page
          if (input$SVM_Build == TRUE){
            
            #ui outpute for title 
            output$Build_SVM_Plot_Title <- renderUI({req(input$Build_ML); h2("Optimal Parameters SVM", align = "center") })
            
            #optimal paraters plot and download 
            output$Build_SVM_Plot <- renderPlot({return(grid.arrange(grobs = svm_plot_list, ncol = 1))})
            svm_plot_grid<-grid.arrange(grobs = svm_plot_list, ncol = 1)
            output$downloadBuild_SVM_Plot <- downloadHandler(filename = function(){paste("SVM_Ensemble_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= svm_plot_grid, width = 12, height = 20)})
            #download button for svm model
            output$downloadBuild_SVM_Mod <- downloadHandler(filename = function(){paste("SVM_Model_",input$Build_SVM_Mod_Num,'.rds',sep='')},
                                                            content = function(fname){
                                                              saveRDS(svm_Mod_list[[match(input$Build_SVM_Mod_Num, name_List)]], fname)})
            
          }
          
          
          
          if (input$LA_Build == TRUE){
            
            #ui title output for Elastic Net page 
            output$Build_LA_Plot_Title <- renderUI({req(input$Build_ML); h2("Optimal Parameters Elastic Net", align = "center") })
            
            #plot output for Elastic Net page 
            output$Build_LA_Plot <- renderPlot({return(grid.arrange(grobs = net_plot_list, ncol = 1))})
            net_plot_grid<-grid.arrange(grobs = net_plot_list, ncol = 1)
            output$downloadBuild_LA_Plot <- downloadHandler(filename = function(){paste("LASSO_Analysis_Plot",'.png',sep='')},
                                                            content = function(file){
                                                              ggsave(file,plot= net_plot_grid, width = 12, height = 20)})
            
            #download button for LA model
            output$downloadBuild_LA_Mod <- downloadHandler(filename = function(){paste("ElasticNet_Model_",input$Build_LA_Mod_Num,'.rds',sep='')},
                                                           content = function(fname){
                                                             saveRDS(la_Mod_list[[match(input$Build_LA_Mod_Num, name_List)]], fname)})
          }
          
          #output for bayesian model 
          if (input$BC_Build == TRUE){
            
            #download button for BC model
            output$downloadBuild_BC_Mod <- downloadHandler(filename = function(){paste("BayesianGLM_Model_",input$Build_BC_Mod_Num,'.rds',sep='')},
                                                           content = function(fname){
                                                             saveRDS(bc_Mod_list[[match(input$Build_BC_Mod_Num, name_List)]], fname)})
          }
          
        })
        
        output$Build_Accuracy_ML_Error<- output$Build_RF_Plot_Error<- output$Build_RF_VIPPlot_Error<-output$Build_XGB_Models_Error<-
          output$Build_XGB_VIPPlot_Error<-output$Build_SVM_Plot_Error<-output$Build_LA_Plot_Error<-output$Build_BC_Mod_Error<-
          renderUI(return())
          
      
    },
    error = function(e) {
      #notification box with specific error
      showNotification(paste0(e), type = 'err')
      #error messages for each tab explaining what could be wrong
      output$Build_Accuracy_ML_Error<- output$Build_RF_Plot_Error<- output$Build_RF_VIPPlot_Error<-output$Build_XGB_Models_Error<-
        output$Build_XGB_VIPPlot_Error<-output$Build_SVM_Plot_Error<-output$Build_LA_Plot_Error<-output$Build_BC_Mod_Error<-
        renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                             " - File Formats and seperator is correct and All inputs filled.", 
                             " - The correct Group column is selected and only 2 groups",
                             " - The correct sample location is selected",
                             " - The first row is the fatty acid family",
                             " - For XGBoost, the checkboxes are clicked for all the fatty acid families in the datasets",
                             sep = "<br/>"))})
      #make table/plots outputs empty
      output$Build_XGB_Plot1<-output$Build_XGB_Plot2<-output$Build_XGB_Plot3<-output$Build_XGB_Plot4<-output$Build_XGB_Plot5<-renderPlotly(return())
      output$Build_Accuracy_ML<-output$Build_RF_Plot<-output$Build_RF_VIPPlot<-output$Build_SVM_Plot<-output$Build_LA_Plot<-
        renderPlot(return())
      output$Build_ML_Table<-output$XGB_Build_Models<-renderDataTable(return())
      output$Build_AccPlot_Title<-output$Build_AccTable_Title<-output$Build_RF_Plot_Title<-output$Build_RF_VIPPlot_Title<-output$Build_XGB_Plot_Title<-
        output$Build_XGB_VIPPlot_Title<-output$Build_SVM_Plot_Title<-output$Build_LA_Plot_Title<-
        renderUI(return())
      
    }, silent=FALSE)
               
    
  })
  
  
  observe({
    
    req(input$LM_ML_File2_5)
    
    #error handling for the code 
    tryCatch(
      
      #try the code to see if it works 
      {
        lm_profile = read.table(input$LM_ML_File2_5$datapath, 
                                sep=input$sep_LM_ML2_5,
                                header = TRUE,
                                row.names = 1,
                                stringsAsFactors = FALSE)
        
        if (input$samp_ML2_5 == "col_ML2_5"){
          
          lm_profile <- data.frame(t(lm_profile))
        }
        
        Met_Rename<-c()
        #remove spaces from column names 
        colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
        
        
        for (x in 1:ncol(lm_profile)){
          if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
            Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
          }else{
            Met_Rename[x] <-colnames(lm_profile)[x]
          }
        }
        colnames(lm_profile)<-Met_Rename
        
        #get the column names and update 
        Col_NamesML<-unique(as.character(unname(unlist(colnames(lm_profile)))))
        #if error in file and column names are empty
        if (Col_NamesML == ""){
          updateSelectInput(session,"Group_ML2_2", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"MetName_ML", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        }
        
        updateSelectInput(session, "Group_ML2_5", choices = Col_NamesML, selected = Col_NamesML[[1]])
        updateSelectInput(session, "MetName_ML", choices = Col_NamesML, selected = Col_NamesML[[1]])
        
      },
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Group_ML2_2", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"MetName_ML", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
    
  })
  
  observeEvent(input$Build_ML2_5, {
    
    req(input$LM_ML_File2_5)
    
    tryCatch(
      {
        # Call the lipid mediator profiling file variant:
        inFile <- input$LM_ML_File2_5
        
        print(inFile)
        
        if(is.null(inFile))
          
          return(NULL)
        
        
        set.seed(415) # To get same results even with the random part.
        options(digits = 3) # To get only part of the decimals.
        
        # Open the file:
        
        
        lm_profile = read.table(inFile$datapath,
                                sep=input$sep_LM_ML2_5,
                                header = TRUE,
                                row.names = 1,
                                stringsAsFactors = FALSE)
        
        if (input$samp_ML2_5 == "col_ML2_5"){
          
          lm_profile <- data.frame(t(lm_profile))
          
          #remove spaces from column names 
          colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
          
          for (x in 1:ncol(lm_profile)){
            if (substr(colnames(lm_profile)[x],1,1) %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
              colnames(lm_profile)[x]= paste("X", colnames(lm_profile)[x], sep = '')
            }
            else{colnames(lm_profile)[x] = colnames(lm_profile)[x]}
          }
          
        }
        
        Met_Rename<-c()
        #remove spaces from column names 
        colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
        
        for (x in 1:ncol(lm_profile)){
          if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
            Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
          }else{
            Met_Rename[x] <-colnames(lm_profile)[x]
          }
        }
        colnames(lm_profile)<-Met_Rename
        
        
        
        #remove fatty acid info row 
        lm_profile<-lm_profile[-1,]
        
        #obtain the responses info column
        ResponseName<-as.character(input$Group_ML2_5)
        Responses<-lm_profile[,ResponseName]
        
        #obtain columns of interest 
        lm_profileMet<-lm_profile[,(names(lm_profile) %in% input$MetName_ML)]
        
        # Creates the data frame that is going to collect all the info regarding the models:
        
        final_table <- data.frame(machine_learning = "del",
                                  Accuracy = 1,
                                  Sensitivity = 1,
                                  Specificity = 1,
                                  TP = 1,
                                  FP = 1,
                                  TN = 1,
                                  FN = 1,
                                  stringsAsFactors = FALSE)
        
        #progress bar
        #withProgress(message = 'Building Models', style = style, value = 0.1, {
        # Sys.sleep(0.25)})
        
        
        #withProgress(message = 'Building Models', value = 0, {
        
        # Save all the values as numeric:
        lm_profile_number <- sapply(lm_profileMet, function(x) as.numeric(x))
        row.names(lm_profile_number) <- row.names(lm_profileMet)
        
        # Scale the data because is better when you are running Machine Learning models:
        lm_profiles_scale <- as.data.frame(scale(lm_profile_number, center = FALSE, scale = TRUE))
        
        # If all the values from one column are the same, the scalation will give you NA. For those cases, to avoid errors
        # replace the NA for zeros.
        lm_profiles_scale[is.na(lm_profiles_scale)] <- 0
        #--------------------
        ### change 0 to 1/5 of lowest value
        lm_profile_zero<-lm_profiles_scale[-1,]
        
        #get vector the same length as number of columns
        cols<- 1:(ncol(lm_profile_zero))
        
        #replace zeros for each column to 1/5 the smallest value for each column
        lm_profile_zero[cols] <- lapply(lm_profile_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
        
        #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
        lm_profile_zero<-as.matrix(lm_profile_zero)
        lm_profile_zero[is.infinite(lm_profile_zero)] <- NA
        min_index <- which.min(lm_profile_zero)
        zero_replace <- (lm_profile_zero[min_index]/5)
        lm_profile_zero <- as.data.frame(lm_profile_zero)
        lm_profile_zero[is.na(lm_profile_zero)] <- zero_replace
        
        #add first row to dataframe
        lm_profiles_scale<-rbind(lm_profiles_scale[1,] ,lm_profile_zero)
        
        #remove columns that have all the same value 
        #lm_profiles_scale<-Filter(var, lm_profiles_scale[-1,])
        
        #---------------------
        # Add the classification variable to the data frame (Responder and non responder):
        
        lm_profiles_scale$responses <- factor(Responses)
        
        # Make sure that column names do not represent a problem to randomForest making them a valid name to R.
        names(lm_profiles_scale) <- make.names(names(lm_profiles_scale))
        
        #---> MACHINE LEARNING (randomForest R):
        if (input$RF_BuildMet == TRUE){
          
          # In this case we defined the lipid mediators to create a loop to define which mtry is the best one for our model.
          oob_error <- double(ncol(lm_profiles_scale) - 1) #Define number of variable. -1 is because the last column is responses.
          for (mtry in 1:(ncol(lm_profiles_scale) - 1)) {
            rf_lm_profiles_scales <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = mtry,
                                                  importance = TRUE, ntree = 10000)
            oob_error[mtry] <- 100 - ((rf_lm_profiles_scales$err.rate[10000])*100)
          }
          # Define the best mtry according to the best prediction value.
          final_mtry <- which.max(oob_error)
          
          # Run the model again with the right mtry value.
          RF_Model <- randomForest(responses ~ ., data = lm_profiles_scale, mtry = final_mtry,
                                   importance = TRUE, ntree = 10000)
          # Get the confusion matrix of the model, sensitivity and specificity:
          RF_ConfMat <- as.data.frame(RF_Model$confusion)
          RF_ConfMat[is.na(RF_ConfMat)] <- 0
          
          # Calculates sensitivity, specificity and AUC.
          Sensi_RF <- RF_ConfMat[2, 2]/(RF_ConfMat[2, 2] + RF_ConfMat[2, 1])
          Speci_RF <- RF_ConfMat[1, 1]/(RF_ConfMat[1, 1] + RF_ConfMat[1, 2])
          
          # Final table for random forest:
          RF_Table <- data.frame(machine_learning = "RF",
                                 Accuracy = 100 - ((RF_Model$err.rate[10000])*100),
                                 Sensitivity = Sensi_RF,
                                 Specificity = Speci_RF,
                                 TP = (1 - RF_ConfMat[2, 3])*100,
                                 FP = RF_ConfMat[1, 3]*100,
                                 TN = (1 - RF_ConfMat[1, 3])*100,
                                 FN = RF_ConfMat[2, 3]*100,
                                 stringsAsFactors = FALSE)
          
          final_table <- rbind(final_table, RF_Table)
          # Number of trees plot:
          tree_table <- data.frame(Ensemble = c(1:10000),
                                   OBB = RF_Model$err.rate[, 1],
                                   AvgAcc = 100 - ((RF_Model$err.rate[, 1])*100))
          
          RF_Plot <- ggplot(data = tree_table, aes(x = Ensemble, y = AvgAcc)) +
            geom_line(linetype = "solid", size = 1.2, color = "navyblue") +
            ggtitle("Random Forest nTrees Plot") +
            scale_x_continuous(name = "Trees") +
            scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
            theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 25, colour = "black"),
                  axis.title.x = element_text(size = 25, colour = "black"),
                  axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep.
                  legend.position = ("none"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(0.4, "cm"))
          #get importance dataframe
          rf_Importance<-as.data.frame(importance(RF_Model, type = 1))
          rf_Importance$lipid_mediators <- rownames(rf_Importance)
          
          Fatty_Acid <- 1:length(rf_Importance$lipid_mediators)
          cbind(rf_Importance, Fatty_Acid)
          
          #fatty acid points
          for (i in 1:length(rf_Importance$lipid_mediators)){
            rf_Importance$Fatty_Acid[i]<-SPMs_FA[rf_Importance$lipid_mediators[i]]
          }
          #change the names to be in the proper format
          rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)] <- paste(gsub("^X", "'", rf_Importance$lipid_mediators[grep("^X", rf_Importance$lipid_mediators)]), "'", sep = "")
          rf_Importance$lipid_mediators <- gsub("4$", "[4]", rf_Importance$lipid_mediators) #converts names ending with 4 to subscript 5
          rf_Importance$lipid_mediators <- gsub("B4'", "B'[4]", rf_Importance$lipid_mediators) #B4 converted to subscript 4
          rf_Importance$lipid_mediators<- gsub("\\.", "-", rf_Importance$lipid_mediators) # Replace "." with "-" when required
          rf_Importance$lipid_mediators <- gsub("n.3.DPA", "[n-3~DPA]", rf_Importance$lipid_mediators) # Transform n-3 DPA as a subscript:
          
          # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
          # (this with the purpose of get the LM in decreasing order in the figure):
          rf_Importance <- rf_Importance[order(rf_Importance$MeanDecreaseAccuracy), ]
          rf_Importance$lipid_mediators <- factor( rf_Importance$lipid_mediators,
                                                   levels = rf_Importance$lipid_mediators[order(rf_Importance$MeanDecreaseAccuracy)])
          
          lipids <- as.character(rf_Importance$lipid_mediators)
          
          #add colors to SPM names
          lm_classes <- as.factor(rf_Importance$Fatty_Acid)
          names(lm_classes) <- rf_Importance$lipid_mediators
          
          dha_index<-which((lm_classes=="DHA")== TRUE)
          n_three_index<-which((lm_classes=="n3DPA")== TRUE)
          epa_index<-which((lm_classes=="EPA")== TRUE)
          aa_index<-which((lm_classes=="AA")== TRUE)
          
          lm_colors <- NULL
          lm_colors[dha_index] <- "blue"
          lm_colors[n_three_index] <- "brown"
          lm_colors[epa_index] <- "darkgoldenrod1"
          lm_colors[aa_index] <- "darkslategray"
          
          #rf_Importance1<- subset(rf_Importance, MeanDecreaseAccuracy >1)
          RF_VIP_Plot <- ggplot(data = rf_Importance, mapping = aes(x = lipid_mediators, y = MeanDecreaseAccuracy, color = Fatty_Acid)) + geom_point(size = 3) +
            scale_y_continuous(name = "Mean Decrease Accuracy") +
            labs(x = "Lipid Mediators", title = "Radnom Forest VIP Plot" )+
            scale_x_discrete(labels = parse(text = lipids)) +
            scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                 "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
            coord_flip() +
            theme(axis.title = element_text(size = 20),
                  axis.text.x  = element_text(size = 15, hjust = 0.5),
                  axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                  legend.position = "top",
                  aspect.ratio = 2/1,
                  legend.title = element_text(size = 20),
                  legend.text  = element_text(size = 15),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
          
          #outputs for random forest 
          output$BuildMet_RF_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("Optimal paramters for Random Forest", align = "center") })
          output$BuildMet_RF_VIPPlot_Title <- renderUI({req(input$Build_ML2_5); h2("Importance of Variance Plot", align = "center") })
          
          output$BuildMet_RF_Plot <- renderPlot(RF_Plot)
          output$BuildMet_RF_VIPPlot <- renderPlot(RF_VIP_Plot)
          
          
          
          output$downloadBuildMet_RF_Plot <- downloadHandler(filename = function(){paste("RF_ntrees_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= RF_Plot, width = 12, height = 8)})
          output$downloadBuildMet_RF_VIPPlot <- downloadHandler(filename = function(){paste("RF_VIP_Plot",'.png',sep='')},
                                                                content = function(file){
                                                                  ggsave(file,plot= RF_VIP_Plot, width = 12, height = 8)})
          output$downloadBuildMet_RF_Mod <- downloadHandler(filename = function(){"RandomForest_Model.rds"},
                                                            content = function(fname){
                                                              saveRDS(RF_Model, fname)})
          output$downloadBuildMet_RF_Mod1 <- downloadHandler(filename = function(){"RandomForest_Model.rds"},
                                                             content = function(fname){
                                                               saveRDS(RF_Model, fname)})
        }
        
        #---> MACHINE LEARNING (EXTREME GRADIENT BOOSTING)
        if (input$XGB_BuildMet == TRUE){
          df<-lm_profiles_scale
          
          #need to convert responses to integer of 0 and 1
          df$responses<-as.factor(df$responses)
          responses<-df$responses
          label<-as.integer(df$responses) -1
          
          #remove labels from dataframe
          df$responses<-NULL
          df<-sapply(df, function(x) as.numeric(x))
          
          #split it into training and testing dataset with 75/20
          n = nrow(df) #get number of rows for dataframe
          train.index = sample(n,floor(0.7*n)) #randomize rows to get dataset
          train.data = as.matrix(df[train.index,])
          train.label = label[train.index]
          test.data = as.matrix(df[-train.index,])
          test.label = label[-train.index]
          
          #transform dataset to xgb matrix
          xgb.train = xgb.DMatrix(data=train.data,label=train.label)
          xgb.test = xgb.DMatrix(data=test.data, label=test.label)
          
          # #create grid with all the options for xgboost
          # searchGridSubCol <- expand.grid(subsample = c(0.5, 0.75, 1),
          #                                 colsample_bytree = c(0.5, 0.75, 1),
          #                                 max_depth = c(2, 4, 6, 8, 10),
          #                                 gamma_val = c(0, 1, 5, 10),
          #                                 eta = c(0.01, 0.1, 0.2, 0.3),
          #                                 min_child = seq(1))
          
          searchGridSubCol <- expand.grid(subsample = c(0.6, 0.8, 1), 
                                          colsample_bytree = c(0.6, 0.8, 1),
                                          max_depth = c(3, 4, 5),
                                          gamma_val = c(0.5, 1, 1.5, 2, 5), 
                                          eta = c(0.1, 0.2, 0.3),
                                          min_child = c(1, 5, 10)
          )
          
          #run through each combination of the grid
          system.time(
            ErrorsHyperparameters <- apply(searchGridSubCol, 1, function(parameterList){
              #Extract Parameters to test
              currentSubsampleRate <- parameterList[["subsample"]]
              currentColsampleRate <- parameterList[["colsample_bytree"]]
              currentDepth <- parameterList[["max_depth"]]
              currentEta <- parameterList[["eta"]]
              currentGamma <- parameterList[["gamma_val"]]
              currentMinChild <- parameterList[["min_child"]]
              
              #run selected parameter through the model
              xgboostModelCV <- xgb.cv(data =  xgb.train, nrounds = 10000, nfold = 5, showsd = TRUE,
                                       metrics = "error", verbose = TRUE, "eval_metric" = "error",
                                       "objective" = "binary:logistic", "max.depth" = currentDepth, "eta" = currentEta,
                                       "subsample" = currentSubsampleRate, "colsample_bytree" = currentColsampleRate, 
                                       print_every_n = 10, booster = "gbtree",
                                       early_stopping_rounds = 10, "gamma" = currentGamma, "min_child_weight" = currentMinChild)
              
              #have error evaluation score as dataframe
              xvalidationScores <- as.data.frame(xgboostModelCV$evaluation_log)
              test_error <- tail(xvalidationScores$test_error_mean, 1)
              train_error <- tail(xvalidationScores$train_error_mean,1)
              output <- return(c(test_error, train_error, currentSubsampleRate, currentColsampleRate, currentDepth, currentEta, currentGamma, currentMinChild))})) #add to final table called output
          
          output_Error <- as.data.frame(t(ErrorsHyperparameters))
          varnames <- c("TestError", "TrainError", "SubSampRate", "ColSampRate", "Depth", "eta", "gamma", "currentMinChild")
          names(output_Error) <- varnames #change the variable names
          top_output<-output_Error[order(output_Error$TestError, output_Error$TrainError), ] #order it by the lowest error rate
          
          #empty table to put the top 10 model information into
          Final_Acc <- data.frame(machine_learning = "test",
                                  Accuracy = 0,
                                  Specificity = 0,
                                  Sensitivity = 0,
                                  TP = 0,
                                  FP = 0,
                                  TN = 0,
                                  FN = 0,
                                  nrounds = 0,
                                  eta= 0,
                                  max_depth= 0,
                                  gamma= 0,
                                  min_child_weight=0,
                                  subsample=0,
                                  colsample_bytree= 0)
          top<-1:50
          #loop to build models for the parameters with the lowest error rate
          for (i in top){
            
            #parameter list with the lowest output
            params <- list(booster = "gbtree", objective = "binary:logistic", eta= top_output[i,6], gamma= top_output[i,7],
                           max_depth= top_output[i,5], min_child_weight=top_output[i,8], subsample=top_output[i,3],
                           colsample_bytree= top_output[i,4])
            
            #Build model based on the training data
            xgb<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
            
            #run model on test data = see how good it is at predicting
            pred<-predict(xgb, test.data)
            
            #obtain confusion matrix
            pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1
            pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
            pred<-as.numeric(pred) #makes it numeric
            test.label_num <-as.numeric(test.label) #makes it numeric
            Responses_List<-unique(Responses)
            pred_y <- cut(pred, 2, label = c(Responses_List[1], Responses_List[2])) #converts to factor of responder or non responder
            test_y <-cut(test.label_num, 2, label = c(Responses_List[1], Responses_List[2])) #convers to factor of responder or non responder
            ConMatrix <- confusionMatrix(reference = test_y, data = pred_y) #build confusion matrix
            TP = as.numeric(ConMatrix$table[2,2])
            FP = as.numeric(ConMatrix$table[2,1])
            TN = as.numeric(ConMatrix$table[1,1])
            FN = as.numeric(ConMatrix$table[1,2])
            
            #convert build table with accuracy for each model
            Table <- data.frame(machine_learning = "XGBoost",
                                Accuracy = as.numeric(ConMatrix$overall[1])*100,
                                Specificity = as.numeric(ConMatrix$byClass[2]),
                                Sensitivity = as.numeric(ConMatrix$byClass[1]),
                                TP = (TP/(TP+FN)*100),
                                FP = (FP/(FP+TN)*100),
                                TN = (TN/(FP+TN)*100),
                                FN = (FN/(TP+FN)*100),
                                nrounds = 10000,
                                eta= top_output[i,6],
                                max_depth= top_output[i,5],
                                gamma= top_output[i,7],
                                min_child_weight=top_output[i,8],
                                subsample=top_output[i,3],
                                colsample_bytree= top_output[i,4])
            Final_Acc<-rbind(Final_Acc, Table)
          }
          
          #remove the first row of random values
          Final_Acc<-Final_Acc[-1,]
          
          #order the top 10 models based on highest accuracy, sensitivity, and specificity
          Top_Acc<-Final_Acc[with(Final_Acc, order(-Accuracy, -Sensitivity, -Specificity)), ]
          
          
          #preparing nrounds plots for top 5 models
          iterations<-1:10000
          top_table<- data.frame(iterations)
          
          params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[1,10], gamma= Top_Acc[1,12],
                         max_depth= Top_Acc[1,11], min_child_weight= Top_Acc[1,13], subsample= Top_Acc[1,14])
          
          #run the cv to get the error rate for each nround
          xgbcv <- xgb.cv( params = params, data = xgb.train, nrounds = 10000, nfold = 5, showsd = T,
                           stratified = T, early_stop_round = 20, maximize = F, metrics = "error", verbose = FALSE)
          # 
          Table1<-data.frame( test = as.numeric(xgbcv$evaluation_log$test_error_mean))
          colnames(Table1)[1] <-"XGBoost_Model"
          top_table<-cbind(top_table, Table1) #add to table
          
          #plot using plotly
          nround_plot<-plot_ly(data = top_table, x = ~iterations)
          nround_plot <-  nround_plot %>% add_trace(y = ~XGBoost_Model, name = colnames(top_table)[2], mode = 'lines')
          nround_plot <- nround_plot %>% layout(title =  "XGBoost Plot", xaxis = list(title = 'Number of Rounds'),
                                                yaxis = list(title = 'Test Error Rate'))
          #makes download quality of plot better
          nround_plot <- nround_plot%>% config(toImageButtonOptions = list(format = "jpeg", width = 1500, height = 750))
          
          
          
          #importance plot
          #build model for the top model
          params <- list(booster = "gbtree", objective = "binary:logistic", eta= Top_Acc[1,10], gamma= Top_Acc[1,12],
                         max_depth= Top_Acc[1,11], min_child_weight= Top_Acc[1,13], subsample= Top_Acc[1,14], colsample_bytree= Top_Acc[1,15])
          xgb_Mod1<-xgb.train(data = xgb.train, params = params, nrounds = 10000, eval.metric = "error", early.stop.rounds=10,  silent = 0)
          
          
          XGB_Importance = xgb.importance(model = xgb_Mod1) #important matrix
          
          Fatty_Acid <- 1:length(XGB_Importance$Feature)
          cbind(XGB_Importance, Fatty_Acid)
          
          #fatty acid points
          for (i in 1:length(XGB_Importance$Feature)){
            XGB_Importance$Fatty_Acid[i]<-SPMs_FA[XGB_Importance$Feature[i]]
          }
          
          # #change the names to be in the proper format
          XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)] <- paste(gsub("^X", "'", XGB_Importance$Feature [grep("^X", XGB_Importance$Feature)]), "'", sep = "")
          XGB_Importance$Feature <- gsub("4$", "[4]", XGB_Importance$Feature) #converts names ending with 4 to subscript 4
          XGB_Importance$Feature <- gsub("B4'", "B'[4]", XGB_Importance$Feature ) #B4 converted to subscript 4
          XGB_Importance$Feature  <- gsub("\\.", "-", XGB_Importance$Feature ) # Replace "." with "-" when required
          XGB_Importance$Feature <- gsub("n.3.DPA", "[n-3~DPA]", XGB_Importance$Feature ) # Transform n-3 DPA as a subscript:
          
          
          # Order the data in increasing order based on the Mean Decrease Accuracy and create a factor column
          # (this with the purpose of get the LM in decreasing order in the figure):
          
          XGB_Importance <- XGB_Importance[order(XGB_Importance$Gain), ]
          
          XGB_Importance$Feature<- factor(XGB_Importance$Feature,
                                          levels = XGB_Importance$Feature[
                                            order(XGB_Importance$Gain)])
          lipids <- as.character(XGB_Importance$Feature )
          
          #add colors to x axis labels
          lm_classes <- as.factor(XGB_Importance$Fatty_Acid)
          names(lm_classes) <- XGB_Importance$Feature
          
          dha_index<-which((lm_classes=="DHA")== TRUE)
          n_three_index<-which((lm_classes=="n3DPA")== TRUE)
          epa_index<-which((lm_classes=="EPA")== TRUE)
          aa_index<-which((lm_classes=="AA")== TRUE)
          
          lm_colors <- NULL
          lm_colors[dha_index] <- "blue"
          lm_colors[n_three_index] <- "brown"
          lm_colors[epa_index] <- "darkgoldenrod1"
          lm_colors[aa_index] <- "darkslategray"
          
          XGB_VIP_Plot <- ggplot(data = XGB_Importance, mapping = aes(x = Feature, y = Gain, color = Fatty_Acid)) + geom_point(size = 3) +
            scale_y_continuous(name = "Gain") +
            labs(x = "Features", title = "XGBoost VIP Plot") +
            scale_x_discrete(labels = parse(text = lipids)) +
            scale_color_manual(name = "Lipid Mediators Metabolomes:", values = c("DHA"="blue","n3DPA"="brown",
                                                                                 "EPA" ="darkgoldenrod1", "AA" = "darkslategray")) +
            coord_flip() +
            theme(axis.title = element_text(size = 20),
                  axis.text.x  = element_text(size = 15, hjust = 0.5),
                  axis.text.y  = element_text(size = 15, hjust = 1, color = lm_colors),
                  legend.position = "top",
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major.y = element_line(size = 0.5, linetype = 2, colour = "black"))
          XGBoost_Table<-Top_Acc[1, 1:8]
          final_table<-rbind(final_table, XGBoost_Table)
          
          #outputs plots and downloads 
          output$BuildMet_XGB_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("XGBoost nrounds vs Test Error Rate Plot", align = "center") })
          output$BuildMet_XGB_VIPPlot_Title <- renderUI({req(input$Build_ML2_5); h2("Importance of Variance Plot", align = "center") })
          
          output$BuildMet_XGB_Plot <- renderPlotly(nround_plot)
          output$BuildMet_XGB_VIPPlot <- renderPlot(XGB_VIP_Plot)
          output$downloadBuildMet_XGB_VIPPlot <- downloadHandler(filename = function(){paste("XGBoost_VIP_Plot",'.png',sep='')},
                                                                 content = function(file){
                                                                   ggsave(file,plot= XGB_VIP_Plot, width = 12, height = 8)})
          output$downloadBuildMet_XGB_Mod <- downloadHandler(filename = function(){"XGBoost_Model.rds"},
                                                             content = function(fname){
                                                               saveRDS(xgb_Mod1, fname)})
          output$downloadBuildMet_XGB_Mod1 <- downloadHandler(filename = function(){"XGBoost_Model.rds"},
                                                              content = function(fname){
                                                                saveRDS(xgb_Mod1, fname)})
        }
        
        
        #---> MACHINE LEARNING (Classyfire R):
        
        # Classyfire uses Support Vector Machine to create the Machine Learning Model. SVM classifies, makes a reggresion,
        # and creates a novelty detection for the creation of the model.
        
        # The idea is to create several models and see which one fits the best. The models will be based on the whole
        # lipid profiles and the different groups based on substrates.
        
        # "cfBuild" to create the SVM:
        # Clasyfire requieres matrix:
        if (input$SVM_BuildMet == TRUE){
          #convert dataframe to matrix
          lm_profiles_scale_matrix <- as.matrix(lm_profiles_scale[, -(ncol(lm_profiles_scale))])
          
          #building SVM model
          support_lmprofiles_scale <- cfBuild(lm_profiles_scale_matrix, lm_profiles_scale$responses,
                                              bootNum = 70, ensNum = 70, cpus = 4)
          #obtaining confusion matrix
          conf_matrix <- as.data.frame(getConfMatr(support_lmprofiles_scale))
          
          # SVM table:
          Support_vector_table <- data.frame(machine_learning = "SVM",
                                             Accuracy = getAvgAcc(support_lmprofiles_scale)$Test,
                                             Sensitivity = conf_matrix[4, 3]/100,
                                             Specificity = conf_matrix[1, 3]/100,
                                             TP = conf_matrix[4, 3],
                                             FP = conf_matrix[2, 3],
                                             TN = conf_matrix[1, 3],
                                             FN = conf_matrix[3, 3],
                                             stringsAsFactors = FALSE)
          final_table <- rbind(final_table, Support_vector_table)
          
          # Ensemble Plot:
          ensAcc   <- getAcc(support_lmprofiles_scale)$Test
          meanVal  <- ensAcc[1]
          for (i in 2:length(ensAcc)) {
            meanVal <- c(meanVal, mean(ensAcc[1:i]))
          }
          ensembl_table <- data.frame(Ensemble = 1:length(support_lmprofiles_scale$testAcc),
                                      AvgAcc = meanVal)
          SVM_Plot <- ggplot(data = ensembl_table, aes(x = Ensemble, y = AvgAcc)) +
            geom_point(aes(colour = AvgAcc), size = 5) +
            geom_line(linetype = "dotted", size = 1) +
            ggtitle("SVM Ensemble Plot") +
            scale_x_continuous(name = "Ensemble interaction") +
            scale_y_continuous(name = "Average Percent Accuracy", limits=c(50, 100)) +
            theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 25, colour = "black"),
                  axis.title.x = element_text(size = 25, colour = "black"),
                  axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1.5), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1.5), # Color and thickness of every axis sep.
                  legend.position = ("none"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(0.4, "cm"))
          
          output$BuildMet_SVM_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("Optimal Parameters SVM", align = "center") })
          output$BuildMet_SVM_Plot <- renderPlot(SVM_Plot)
          output$downloadBuildMet_SVM_Plot <- downloadHandler(filename = function(){paste("SVM_Ensemble_Plot",'.png',sep='')},
                                                              content = function(file){
                                                                ggsave(file,plot= SVM_plot, width = 12, height = 8)})
          output$downloadBuildMet_SVM_Mod <- downloadHandler(filename = function(){"SVM_Model.rds"},
                                                             content = function(fname){
                                                               saveRDS(support_lmprofiles_scale, fname)})
          
        }
        
        #---> ELASTIC NET REGRESSION (caret R):
        if (input$LA_BuildMet == TRUE){
          # Get the explanatory variables as a matrix:
          explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
          
          # LASSO Analysis:
          model_net <- train(explanatory, lm_profiles_scale$responses, method = "glmnet",
                             trControl = trainControl("boot", number = 70))
          # Get the confusion matrix of the model:
          conf_net <- as.data.frame(confusionMatrix(model_net, "none")$table)
          
          # Final Elastic net model table:
          net_table <- data.frame(machine_learning = "ElasticNet",
                                  Accuracy = max(model_net$results$Accuracy)*100,
                                  Sensitivity = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3])),
                                  Specificity = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3])),
                                  TP = (conf_net[4, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                  FP = (conf_net[2, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                  TN = (conf_net[1, 3]/(conf_net[1, 3] + conf_net[2, 3]))*100,
                                  FN = (conf_net[3, 3]/(conf_net[4, 3] + conf_net[3, 3]))*100,
                                  stringsAsFactors = FALSE)
          final_table <- rbind(final_table, net_table)
          
          # Parameter tuning figure:
          scaleFUN <- function(x) sprintf("%.2f", x)
          
          Elastic_Net_Plot <- ggplot(model_net, highlight = TRUE) +
            scale_x_continuous(name = expression(paste("Alpha (", alpha, ")", sep = ""))) +
            scale_y_continuous(name = "Average Accuracy", labels= scaleFUN) +
            ggtitle("Elastic Net Parameters Plot") +
            scale_color_manual(values = c("darkorchid3", "orangered1", "chartreuse3")) +
            scale_shape_manual(values=c(16, 16, 16)) +
            labs(color = expression(paste("Lambda (", lambda, ")", sep = ""))) +
            guides(shape = FALSE) +
            theme(plot.title = element_text(size = 30, colour = "black", hjust = 0.5),
                  axis.title.y = element_text(size = 25, colour = "black"),
                  axis.title.x = element_text(size = 25, colour = "black"),
                  axis.text.x  =  element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.text.y  = element_text(size = 25, colour = "black"), # Put color to the labels
                  axis.line = element_line(colour = 'black', size = 1), # Color and thickness of axis
                  axis.ticks = element_line(colour = "black", size = 1), # Color and thickness of every axis sep.
                  legend.title = element_text(size = 25),
                  legend.text  = element_text(size = 20),
                  legend.key = element_rect(fill = NA),
                  legend.key.size = unit(1.3, "cm"),
                  panel.background = element_rect(fill = "white"),
                  axis.ticks.length = unit(0.4, "cm"))
          
          output$BuildMet_LA_Plot_Title <- renderUI({req(input$Build_ML2_5); h2("Optimal Parameters Elastic Net Regression", align = "center") })
          output$BuildMet_LA_Plot <- renderPlot(Elastic_Net_Plot)
          output$downloadBuildMet_LA_Plot <- downloadHandler(filename = function(){paste("ElasticNet_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= Elastic_Net_Plot, width = 12, height = 8)})
          output$downloadBuildMet_LA_Mod <- downloadHandler(filename = function(){"ElasticNet_Model.rds"},
                                                            content = function(fname){
                                                              saveRDS(model_net, fname)})
          
        }
        
        #---> BAYESIAN MODEL (Caret R):
        
        if (input$BC_BuildMet == TRUE){
          explanatory <- data.matrix(lm_profiles_scale)[, -(ncol(lm_profiles_scale))]
          
          bayesian <- train(explanatory, lm_profiles_scale$responses, method = "bayesglm",
                            trControl = trainControl("boot", number = 70))
          conf_bay <- as.data.frame(confusionMatrix(bayesian, "none")$table)
          
          # Final Bayesian model table:
          bay_table <- data.frame(machine_learning = "Bayes GLM",
                                  Accuracy = max(bayesian$results$Accuracy)*100,
                                  Sensitivity = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3])),
                                  Specificity = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3])),
                                  TP = (conf_bay[4, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                  FP = (conf_bay[2, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                  TN = (conf_bay[1, 3]/(conf_bay[1, 3] + conf_bay[2, 3]))*100,
                                  FN = (conf_bay[3, 3]/(conf_bay[4, 3] + conf_bay[3, 3]))*100,
                                  stringsAsFactors = FALSE)
          final_table <- rbind(final_table, bay_table)
          
          output$downloadBuildMet_BC_Mod <- downloadHandler(filename = function(){"BayesGLM_Model.rds"},
                                                            content = function(fname){
                                                              saveRDS(bayesian, fname)})
        }
        
        final_table <- final_table[-1, ]
        final_table$Accuracy <- round(final_table$Accuracy, 0)
        final_table$Sensitivity <- round(final_table$Sensitivity, 2)
        final_table$Specificity <- round(final_table$Specificity, 2)
        final_table$TP <- round(final_table$TP, 0)
        final_table$FP <- round(final_table$FP, 0)
        final_table$TN <- round(final_table$TN, 0)
        final_table$FN <- round(final_table$FN, 0)
        
        accuracy <- ggplot(data = final_table, aes(x = machine_learning, y = Accuracy, fill = "#00C19A")) +
          geom_bar(stat = "identity") +
          geom_text(position = position_dodge(w = 0.9), vjust = -0.7,
                    aes(label = paste(Accuracy, "%", sep = "")), size = 8) +
          scale_y_continuous(name = "% Accuracy Score", breaks = seq(from = 0, to = 100, by = 20),
                             expand = c(0, 1)) + # Expand to delete the space between zero and the x-axis. 
          coord_cartesian(ylim = c(1, 120)) +
          theme(axis.title = element_text(size = 20),
                axis.title.x = element_blank(),
                axis.text.x  =  element_text(size = 20, hjust = 0.5, colour = "black"), # Put color to the labels
                axis.text.y  = element_text(size = 20, hjust = 1, colour = "black"), # Put color to the labels
                axis.line = element_line(colour = 'black', size = 1.0), # Color and thickness of axis
                axis.ticks = element_line(colour = "black", size = 1.0), # Color and thickness of every axis sep. 
                legend.key = element_blank(),
                panel.background = element_rect(fill = "white"),
                axis.ticks.length = unit(0.4, "cm"))
        
        
        # Final table cute:
        
        colnames(final_table)[1] <- "Machine Learning Methodology"
        
        output$BuildMet_AccPlot_Title <- renderUI({req(input$Build_ML2_5); h2("% Accuracy Score Figure for different ML models", align = "center") })
        output$BuildMet_AccPlot <- renderPlot(accuracy)
        output$downloadBuildMet_AccPlot <- downloadHandler(filename = function(){paste("Accuracy_Plot",'.png',sep='')},
                                                           content = function(file){
                                                             ggsave(file, plot= accuracy, width = 12, height = 8)})
        
        output$BuildMet_ML_Table <- renderDataTable({return(final_table)})
        output$downloadBuildMet_Acc_ML_Table <- downloadHandler(filename = function(){"Accuracy_Table.csv"},
                                                                content = function(fname){
                                                                  write.csv(final_table, fname)})
        
        output$BuildMet_Accuracy_ML_Error<- output$BuildMet_RF_Plot_Error<- output$BuildMet_RF_VIPPlot_Error<-output$BuildMet_XGB_Models_Error<-
          output$BuildMet_XGB_VIPPlot_Error<-output$BuildMet_SVM_Plot_Error<-output$BuildMet_LA_Plot_Error<-output$BuildMet_BC_Mod_Error<-
          renderUI(return())
      },
      error = function(e) {
        #notification box with specific error
        showNotification(paste0(e), type = 'err')
        #error messages for each tab explaining what could be wrong
        output$BuildMet_Accuracy_ML_Error<- output$BuildMet_RF_Plot_Error<- output$BuildMet_RF_VIPPlot_Error<-output$BuildMet_XGB_Models_Error<-
          output$BuildMet_XGB_VIPPlot_Error<-output$BuildMet_SVM_Plot_Error<-output$BuildMet_LA_Plot_Error<-output$BuildMet_BC_Mod_Error<-
          renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                               " - File Formats and seperator is correct and All inputs filled.", 
                               " - The correct Group column is selected and only 2 groups",
                               " - The correct sample location is selected",
                               " - The first row is the fatty acid family",
                               sep = "<br/>"))})
        #make table/plots outputs empty
        output$BuildMet_XGB_Plot<-renderPlotly(return())
        output$BuildMet_Accuracy_ML<-output$BuildMet_RF_Plot<-output$Build_RFMet_VIPPlot<-output$BuildMet_SVM_Plot<-output$BuildMet_LA_Plot<-
          renderPlot(return())
        output$BuildMet_ML_Table<-renderDataTable(return())
        output$BuildMet_AccPlot_Title<-output$BuildMet_RF_Plot_Title<-output$BuildMet_RF_VIPPlot_Title<-output$BuildMet_XGB_Plot_Title<-
          output$Build_XGBMet_VIPPlot_Title<-output$BuildMet_SVM_Plot_Title<-output$BuildMet_LA_Plot_Title<-
          renderUI(return())
        
      }, silent=FALSE)


    
      



    })
  
  observe({
    
    req(input$LM_ML_File3)
    
    
    tryCatch(
      
      #try the code to see if it works 
      {
        lm_profile = read.table(input$LM_ML_File3$datapath, 
                                sep=input$sep_LM_ML3,
                                header = TRUE,
                                row.names = 1,
                                stringsAsFactors = FALSE)
        
        if (input$samp_ML3 == "col_ML3"){
          
          lm_profile <- data.frame(t(lm_profile))
          
          #remove spaces from column names 
          colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile)) 
          
          for (x in 1:ncol(lm_profile)){
            if (substr(colnames(lm_profile)[x],1,1) %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
              colnames(lm_profile)[x]= paste("X", colnames(lm_profile)[x], sep = '')
            }
            else{colnames(lm_profile)[x] = colnames(lm_profile)[x]}
          }
          
        }
        
        Met_Rename<-c()
        #remove spaces from column names 
        colnames(lm_profile)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profile))
        
        for (x in 1:ncol(lm_profile)){
          if (colnames(lm_profile)[x] %in% Rename_Met_Table$Other_Name){
            Met_Rename[x] <-Rename_Met[[colnames(lm_profile)[x]]]
          }else{
            Met_Rename[x] <-colnames(lm_profile)[x]
          }
        }
        colnames(lm_profile)<-Met_Rename
        
        #get the column names and update 
        Col_NamesML<-unique(as.character(unname(unlist(colnames(lm_profile)))))
        
        #if error in file and column names are empty
        if (Col_NamesML == ""){
          updateSelectInput(session,"Group_ML3", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
          updateSelectInput(session,"MetName_ML3", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        }
        updateSelectInput(session, "Group_ML3", choices = Col_NamesML, selected = Col_NamesML[[1]])
        updateSelectInput(session, "MetName_ML3", choices = Col_NamesML, selected = Col_NamesML[[1]])
        
      },
      #if there is an error in the above code, the error message will display 
      error = function(e) {
        showNotification('Error! Please check inputs are correct','',type = "error")
        updateSelectInput(session,"Group_ML3", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        updateSelectInput(session,"MetName_ML3", choices = "Error!!Please check file and separator", selected = "Error!!Please check file and Separator")
        
      }, silent= FALSE)
    
    
    
    
    
  })


  observeEvent(input$Run_ML,{
    
    req(input$LM_ML_File3)
    tryCatch(
      {
        # Call the lipid mediator profiling file variant:
        inFile <- input$LM_ML_File3
        set.seed(415) # To get same results even with the random part.
        options(digits = 3) # To get only part of the decimals.
        
        # Open the file:
        lm_profiles = read.csv(input$LM_ML_File3$datapath, 
                               sep=input$sep_LM_ML3,
                               header = TRUE,
                               row.names = 1,
                               stringsAsFactors = FALSE)
        
        if (input$samp_ML3 == "col_ML3"){
          
          lm_profiles <- data.frame(t(lm_profiles))
          
          #remove spaces from column names 
          colnames(lm_profiles)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profiles)) 
          
          for (x in 1:ncol(lm_profiles)){
            if (substr(colnames(lm_profiles)[x],1,1) %in% c("0","1", "2", "3", "4", "5", "6", "7", "8", "9")){
              colnames(lm_profiles)[x]= paste("X", colnames(lm_profiles)[x], sep = '')
            }
            else{colnames(lm_profiles)[x] = colnames(lm_profiles)[x]}
          }
          
        }
        Met_Rename<-c()
        #remove spaces from column names 
        colnames(lm_profiles)<-gsub("[^a-zA-Z0-9]","", colnames(lm_profiles))
        
        for (x in 1:ncol(lm_profiles)){
          if (colnames(lm_profiles)[x] %in% Rename_Met_Table$Other_Name){
            Met_Rename[x] <-Rename_Met[[colnames(lm_profiles)[x]]]
          }else{
            Met_Rename[x] <-colnames(lm_profiles)[x]
          }
        }
        colnames(lm_profiles)<-Met_Rename
        
        
        #obtain the responses info column
        ResponseName<-as.character(input$Group_ML3)
        
        #Zero handing and scaling data 
        label <- lm_profiles[,ResponseName]
        if (input$Choose_MetName3 == TRUE){
          Val_Data_Raw<- lm_profiles[,(names(lm_profiles) %in% input$MetName_ML3)]
        } else{
          Val_Data_Raw<- lm_profiles[,!(names(lm_profiles) %in% ResponseName)]
        }
        
        Val_Data_Raw<-sapply(Val_Data_Raw, function(x) as.numeric(x)) #convert to numeric
        row.names(Val_Data_Raw) <- row.names(label)
        Val_Data_Scale<-as.data.frame(scale(Val_Data_Raw, center = FALSE, scale = TRUE)) #scale data
        Val_Data_Scale[is.na(Val_Data_Scale)] <- 0 #convert NA to 0
        
        #zero handling to make 0 to 1/5 lowest value
        
        Val_Data_zero<-Val_Data_Scale
        
        
        #get vector the same length as number of columns
        cols<- 1:(ncol(Val_Data_zero))
        
        
        #replace zeros for each column to 1/5 the smallest value for each column
        
        Val_Data_zero[cols] <- lapply(Val_Data_zero[cols], function(x) replace(x, x == 0, min(x[x>0], na.rm = TRUE)/5))
        
        #if column only has zero, it will produce inf therefore need to replace it with 1/5 lowest value in df
        Val_Data_zero<-as.matrix(Val_Data_zero)
        Val_Data_zero[is.infinite(Val_Data_zero)] <- NA
        min_index <- which.min(Val_Data_zero)
        zero_replace <- (Val_Data_zero[min_index]/5)
        Val_Data_zero <- as.data.frame(Val_Data_zero)
        Val_Data_zero[is.na(Val_Data_zero)] <- zero_replace
        
        Val_Data<-Val_Data_zero
        
        #remove columns with all the same values 
        #Val_Data<-Filter(var, Val_Data_zero)
        
        #make sure column names do not cause issue in random forest
        names(label) <- make.names(names(label))
        names(Val_Data) <- make.names(names(Val_Data))
        
        # #No data preprocessing 
        # label <- lm_profiles[,1]
        # Val_Data<- lm_profiles[,-1]
        # Val_Data<-sapply(Val_Data, function(x) as.numeric(x)) #convert to numeric
        
        #output$Run_ValData<-renderTable(Val_Data)
        ROC_List<- list()
        ROC_Names<-list()
        pos = 1
        
        #final ROC Table
        ROC_Table <- data.frame(Model = "del",
                                AUC = 1,
                                Accuracy = 1,
                                Sensitivity = 1,
                                Specificity = 1,
                                TP_per = 1,
                                FP_per = 1,
                                TN_per = 1,
                                FN_per = 1,
                                stringsAsFactors = FALSE)
        
        if (input$RF_Run == TRUE){
          for (x in 1:input$numRF_Model){
            if (x == 1){RF_Model <- readRDS(input$RF_Mod_File1$datapath)}
            if (x == 2){RF_Model <- readRDS(input$RF_Mod_File2$datapath)}
            if (x == 3){RF_Model <- readRDS(input$RF_Mod_File3$datapath)}
            if (x == 4){RF_Model <- readRDS(input$RF_Mod_File4$datapath)}
            if (x == 5){RF_Model <- readRDS(input$RF_Mod_File5$datapath)}
            
            pred_RF<-as.data.frame(predict(RF_Model, Val_Data, type = 'prob'))
            RF_ROC<-pROC::roc(label,pred_RF$Non_Responder)
            
            RF_ROCPlot<- ggroc(RF_ROC, legacy.axes = TRUE) +
              geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
              theme_bw()
            
            output$Run_RF_ROCPlot<-renderPlot({return(RF_ROCPlot)})
            
            output$downloadRun_RF_ROCPlot <- downloadHandler(filename = function(){paste("RandomForest_ROC_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= RF_ROCPlot, width = 15, height = 10)})
            
            # ROC_List<-append(ROC_List, RF_ROC)
            # ROC_Names<-append(ROC_Names, "Random Forest")
            
            ROC_List[[pos]] <- RF_ROC
            ROC_Names[[pos]] <- paste("RF_Mod", x, sep = '')
            pos = pos + 1
            
            #get table information---confusion matrix 
            # Get the confusion matrix of the model, sensitivity and specificity: 
            RF_ConfMat<- as.data.frame(RF_Model$confusion)
            RF_ConfMat[is.na(RF_ConfMat)] <- 0
            
            # Calculates sensitivity, specificity and AUC.
            RF_Sensi <- RF_ConfMat[2, 2]/(RF_ConfMat[2, 2] + RF_ConfMat[2, 1])
            RF_Speci <- RF_ConfMat[1, 1]/(RF_ConfMat[1, 1] + RF_ConfMat[1, 2])
            
            # Final table for random forest:
            RF_Table <- data.frame(Model = paste("Random Forest Model", x, sep = ''),
                                   AUC = auc(RF_ROC),
                                   Accuracy = 100 - ((RF_Model$err.rate[10000])*100),
                                   Sensitivity = RF_Sensi,
                                   Specificity = RF_Speci,
                                   TP_per = (1 - RF_ConfMat[2, 3])*100,
                                   FP_per = RF_ConfMat[1, 3]*100,
                                   TN_per = (1 - RF_ConfMat[1, 3])*100,
                                   FN_per = RF_ConfMat[2, 3]*100,
                                   stringsAsFactors = FALSE)
            ROC_Table <- rbind(ROC_Table, RF_Table)
          }
          
        }
        
        if (input$XGB_Run == TRUE){
          for (x in 1:input$numXGB_Model){
            if (x == 1){XGB_Model <- readRDS(input$XGB_Mod_File1$datapath)}
            if (x == 2){XGB_Model <- readRDS(input$XGB_Mod_File2$datapath)}
            if (x == 3){XGB_Model <- readRDS(input$XGB_Mod_File3$datapath)}
            if (x == 4){XGB_Model <- readRDS(input$XGB_Mod_File4$datapath)}
            if (x == 5){XGB_Model <- readRDS(input$XGB_Mod_File5$datapath)}
            
            pred_XGB<-as.data.frame(predict(XGB_Model, as.matrix(Val_Data),  type = 'prob'))
            XGB_ROC<-pROC::roc(label, pred_XGB[,1])
            
            XGB_ROCPlot<- ggroc(XGB_ROC, legacy.axes = TRUE) +
              geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
              theme_bw()
            
            #XGB_ROCPlot<-plot(XGB_ROC)
            
            output$Run_XGB_ROCPlot<-renderPlot({return(XGB_ROCPlot)})
            
            output$downloadRun_XGB_ROCPlot <- downloadHandler(filename = function(){paste("XGBoost_ROC_Plot",'.png',sep='')},
                                                              content = function(file){
                                                                ggsave(file,plot= XGB_ROCPlot, width = 15, height = 10)})
            
            # ROC_List<-append(ROC_List, XGB_ROC)
            # ROC_Names<-append(ROC_Names, "XGBoost")
            ROC_List[[pos]] <- XGB_ROC
            ROC_Names[[pos]] <- paste("XGB_Mod", x, sep = '')
            pos = pos + 1
            
            #get table information for Confusion matrix
            pred<-pred_XGB[,1]
            pred[(pred>0.5)] = 1 #any value more than 0.5 is rounded up to 1 
            pred[(pred<0.5)] = 0 #any value less than 0.5 is rounded down to 0
            pred<-as.numeric(pred) #makes it numeric
            label_xgb<-as.factor(label)
            label_xgb<-as.integer(label_xgb) -1
            label_xgb_num<-as.numeric(label_xgb)
            label_list<-unique(label) #list of unique values 
            pred_y <- cut(pred, 2, label = c(label_list[1], label_list[2])) #converts to factor of responder or non responder 
            label_y <-cut(label_xgb_num, 2, label = c(label_list[1], label_list[2])) #convers to factor of responder or non responder 
            XGB_ConfMat <- confusionMatrix(reference = label_y, data = pred_y) #build confusion matrix
            
            TP = as.numeric(XGB_ConfMat$table[2,2])
            FP = as.numeric(XGB_ConfMat$table[2,1])
            TN = as.numeric(XGB_ConfMat$table[1,1])
            FN = as.numeric(XGB_ConfMat$table[1,2])
            
            #convert build table with accuracy for each model
            XGB_Table <- data.frame(Model = paste("Extreme Gradient Boosting", x, sep = ''),
                                    AUC = auc(XGB_ROC),
                                    Accuracy = as.numeric(XGB_ConfMat$overall[1])*100,
                                    Specificity = as.numeric(XGB_ConfMat$byClass[2]),
                                    Sensitivity = as.numeric(XGB_ConfMat$byClass[1]),
                                    TP_per = (TP/(TP+FN)*100),
                                    FP_per = (FP/(FP+TN)*100),
                                    TN_per = (TN/(FP+TN)*100),
                                    FN_per = (FN/(TP+FN)*100),
                                    stringsAsFactors = FALSE)
            
            ROC_Table <- rbind(ROC_Table, XGB_Table)
          }
          
        }
        
        if (input$SVM_Run == TRUE){
          
          for (x in 1:input$numSVM_Model){
            if (x == 1){SVM_Model <- readRDS(input$SVM_Mod_File1$datapath)}
            if (x == 2){SVM_Model <- readRDS(input$SVM_Mod_File2$datapath)}
            if (x == 3){SVM_Model <- readRDS(input$SVM_Mod_File3$datapath)}
            if (x == 4){SVM_Model <- readRDS(input$SVM_Mod_File4$datapath)}
            if (x == 5){SVM_Model <- readRDS(input$SVM_Mod_File5$datapath)}
            
            pred_SVM<-cfPredict(SVM_Model, Val_Data)
            names(pred_SVM)<-c("prediction", "Coef Score")
            label_list<-unique(label) #list of unique values (the group names )
            
            pred_SVM$Responder[pred_SVM$prediction == label_list[1]] <- pred_SVM$`Coef Score`[pred_SVM$prediction == label_list[1]]/100 
            pred_SVM$Responder[pred_SVM$prediction == label_list[2]] <- 1 - (pred_SVM$`Coef Score`[pred_SVM$prediction == label_list[2]]/100)
            pred_SVM$Non_Responder <- 1 - pred_SVM$Responder
            
            SVM_ROC<-pROC::roc(label,pred_SVM$Non_Responder)
            
            SVM_ROCPlot<- ggroc(SVM_ROC, legacy.axes = TRUE) +
              geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
              theme_bw()
            
            output$Run_SVM_ROCPlot<-renderPlot({return(SVM_ROCPlot)})
            
            output$downloadRun_SVM_ROCPlot <- downloadHandler(filename = function(){paste("SVM_ROC_Plot",'.png',sep='')},
                                                              content = function(file){
                                                                ggsave(file,plot= SVM_ROCPlot, width = 15, height = 10)})
            
            
            # ROC_List<-append(ROC_List, SVM_ROC)
            # ROC_Names<-append(ROC_Names, "SVM")
            
            ROC_List[[pos]] <- SVM_ROC
            ROC_Names[[pos]] <- paste("SVM_Mod", x, sep = '')
            pos = pos + 1
            
            #get  confusion matrix informaiton
            #obtaining confusion matrix 
            SVM_ConfMat <- as.data.frame(getConfMatr(SVM_Model))
            
            # SVM table:
            SVM_Table <- data.frame(Model = paste("Support Vector Machine", x, sep = ''),
                                    AUC = auc(SVM_ROC),
                                    Accuracy = getAvgAcc(SVM_Model)$Test,
                                    Sensitivity = SVM_ConfMat[4, 3]/100,
                                    Specificity = SVM_ConfMat[1, 3]/100,
                                    TP_per = SVM_ConfMat[4, 3],
                                    FP_per = SVM_ConfMat[2, 3],
                                    TN_per = SVM_ConfMat[1, 3],
                                    FN_per = SVM_ConfMat[3, 3],
                                    stringsAsFactors = FALSE)
            ROC_Table <- rbind(ROC_Table, SVM_Table)
            
          }
          
        }
        
        
        if (input$LA_Run == TRUE){
          
          for (x in 1:input$numLA_Model){
            if (x == 1){LA_Model <- readRDS(input$LA_Mod_File1$datapath)}
            if (x == 2){LA_Model <- readRDS(input$LA_Mod_File2$datapath)}
            if (x == 3){LA_Model <- readRDS(input$LA_Mod_File3$datapath)}
            if (x == 4){LA_Model <- readRDS(input$LA_Mod_File4$datapath)}
            if (x == 5){LA_Model <- readRDS(input$LA_Mod_File5$datapath)}
            
            
            pred_LA<-as.data.frame(predict(LA_Model, Val_Data, type = 'prob'))
            LA_ROC<-pROC::roc(label,pred_LA$Non_Responder)
            
            LA_ROCPlot<- ggroc(LA_ROC, legacy.axes = TRUE) +
              geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
              theme_bw()
            
            output$Run_LA_ROCPlot<-renderPlot({return(LA_ROCPlot)})
            
            output$downloadRun_LA_ROCPlot <- downloadHandler(filename = function(){paste("ElasticNet_ROC_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= LA_ROCPlot, width = 15, height = 10)})
            
            
            # ROC_List<-append(ROC_List, LA_ROC)
            # ROC_Names<-append(ROC_Names, "Elastic Net")
            
            ROC_List[[pos]] <- LA_ROC
            ROC_Names[[pos]] <- paste("ENR_Mod", x, sep = '')
            pos = pos + 1
            
            #confusion matrxi
            LA_ConfMat <- as.data.frame(confusionMatrix(LA_Model, "none")$table)
            
            # Final Elastic net model table:
            LA_Table <- data.frame(Model = paste("Elastic Net Liner Regression Model", x, sep = ''),
                                   AUC = auc(LA_ROC),
                                   Accuracy = max(LA_Model$results$Accuracy)*100,
                                   Sensitivity = (LA_ConfMat[4, 3]/(LA_ConfMat[4, 3] + LA_ConfMat[3, 3])),
                                   Specificity = (LA_ConfMat[1, 3]/(LA_ConfMat[1, 3] + LA_ConfMat[2, 3])),
                                   TP_per = (LA_ConfMat[4, 3]/(LA_ConfMat[4, 3] + LA_ConfMat[3, 3]))*100,
                                   FP_per = (LA_ConfMat[2, 3]/(LA_ConfMat[1, 3] + LA_ConfMat[2, 3]))*100,
                                   TN_per = (LA_ConfMat[1, 3]/(LA_ConfMat[1, 3] + LA_ConfMat[2, 3]))*100,
                                   FN_per = (LA_ConfMat[3, 3]/(LA_ConfMat[4, 3] + LA_ConfMat[3, 3]))*100,
                                   stringsAsFactors = FALSE)
            
            ROC_Table <- rbind(ROC_Table, LA_Table)
            
          }
          
          
          
        }
        
        if (input$BC_Run == TRUE){
          
          for (x in 1:input$numBC_Model){
            if (x == 1){BC_Model <- readRDS(input$BC_Mod_File1$datapath)}
            if (x == 2){BC_Model <- readRDS(input$BC_Mod_File2$datapath)}
            if (x == 3){BC_Model <- readRDS(input$BC_Mod_File3$datapath)}
            if (x == 4){BC_Model <- readRDS(input$BC_Mod_File4$datapath)}
            if (x == 5){BC_Model <- readRDS(input$BC_Mod_File5$datapath)}
            
            pred_BC<-as.data.frame(predict(BC_Model, Val_Data, type = 'prob'))
            BC_ROC<-pROC::roc(label,pred_BC$Non_Responder)
            
            BC_ROCPlot<- ggroc(BC_ROC, legacy.axes = TRUE) +
              geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
              theme_bw()
            
            output$Run_BC_ROCPlot<-renderPlot({return(BC_ROCPlot)})
            
            output$downloadRun_BC_ROCPlot <- downloadHandler(filename = function(){paste("BayesGLM_ROC_Plot",'.png',sep='')},
                                                             content = function(file){
                                                               ggsave(file,plot= BC_ROCPlot, width = 15, height = 10)})
            
            
            # ROC_List<-append(ROC_List, BC_ROC)
            # ROC_Names<-append(ROC_Names, "Bayes")
            # 
            ROC_List[[pos]] <- BC_ROC
            ROC_Names[[pos]] <- paste("BC_Mod", x, sep = '')
            pos = pos + 1
            
            #confusion matrix
            BC_ConfMat <- as.data.frame(confusionMatrix(BC_Model, "none")$table)
            
            # Final Elastic net model table:
            BC_Table <- data.frame(Model = paste("Bayes Linear Regression Model", x, sep =''),
                                   AUC = auc(BC_ROC),
                                   Accuracy = max(BC_Model$results$Accuracy)*100,
                                   Sensitivity = (BC_ConfMat[4, 3]/(BC_ConfMat[4, 3] + BC_ConfMat[3, 3])),
                                   Specificity = (BC_ConfMat[1, 3]/(BC_ConfMat[1, 3] + BC_ConfMat[2, 3])),
                                   TP_per = (BC_ConfMat[4, 3]/(BC_ConfMat[4, 3] + BC_ConfMat[3, 3]))*100,
                                   FP_per = (BC_ConfMat[2, 3]/(BC_ConfMat[1, 3] + BC_ConfMat[2, 3]))*100,
                                   TN_per = (BC_ConfMat[1, 3]/(BC_ConfMat[1, 3] + BC_ConfMat[2, 3]))*100,
                                   FN_per = (BC_ConfMat[3, 3]/(BC_ConfMat[4, 3] + BC_ConfMat[3, 3]))*100,
                                   stringsAsFactors = FALSE)
            
            ROC_Table <- rbind(ROC_Table, BC_Table)
            
          }
          
        }
        
        
        Final_ROCPlot<-ggroc(ROC_List, legacy.axes = TRUE) +
          geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey", linetype="dashed", size = 0.5) +
          scale_color_discrete(name = "ML_Models", labels = ROC_Names) +
          theme_bw()
        
        output$ROC_Plot<-renderPlot({return(Final_ROCPlot)})
        
        output$downloadRun_ROCPlot <- downloadHandler(filename = function(){paste("ROC_Plot",'.png',sep='')},
                                                      content = function(file){
                                                        ggsave(file,plot= Final_ROCPlot, width = 15, height = 10)})
        
        
        #round table and remove first column
        ROC_Table<-ROC_Table[-1,]
        ROC_Table$Accuracy <- round(as.numeric(ROC_Table$Accuracy), 0)
        ROC_Table$AUC <- round(as.numeric(ROC_Table$AUC), 2)
        ROC_Table$Sensitivity <- round(ROC_Table$Sensitivity, 2)
        ROC_Table$Specificity <- round(ROC_Table$Specificity, 2)
        ROC_Table$TP_per <- round(as.numeric(ROC_Table$TP_per), 0)
        ROC_Table$FP_per <- round(as.numeric(ROC_Table$FP_per), 0)
        ROC_Table$TN_per <- round(as.numeric(ROC_Table$TN_per), 0)
        ROC_Table$FN_per <- round(as.numeric(ROC_Table$FN_per), 0)
        
        output$ROC_Table<-renderDataTable(ROC_Table, rownames = FALSE)
        
        output$downloadRun_ROCTable <- downloadHandler(filename = function(){"ROC_Accuracy_Table.csv"}, 
                                                       content = function(fname){
                                                         write.csv(ROC_Table, fname)})
        
        output$ROC_Plot_Error<- output$Run_RF_ROCPlot_Error<- output$Run_XGB_ROCPlot_Error<-output$Run_SVM_ROCPlot_Error<-
          output$Run_LA_ROCPlot_Error<-output$Run_BC_ROCPlot_Error<-
          renderUI(return())
      },
      error = function(e) {
        #notification box with specific error
        showNotification(paste0(e), type = 'err')
        #error messages for each tab explaining what could be wrong
        output$ROC_Plot_Error<- output$Run_RF_ROCPlot_Error<- output$Run_XGB_ROCPlot_Error<-output$Run_SVM_ROCPlot_Error<-
          output$Run_LA_ROCPlot_Error<-output$Run_BC_ROCPlot_Error<-
          renderUI({HTML(paste("Error!!Please confirm if the following is correct:",
                               " - File Formats and seperator is correct and All inputs filled.", 
                               " - The correct Group column is selected and only 2 groups",
                               " - The correct sample location is selected",
                               " - The first row is the fatty acid family",
                               " - The model is uploaded as RDS file",
                               sep = "<br/>"))})
        #make table/plots outputs empty
        output$ROC_Plot<-output$Run_RF_ROCPlot<-output$Run_XGB_ROCPlot<-output$Run_SVM_ROCPlot<-
          output$Run_LA_ROCPlot<-output$Run_BC_ROCPlot<-
          renderPlot(return())
        output$ROC_Table<-renderDataTable(return())
        
        
      }, silent=FALSE)
    
    
    
  })
  
  
  
  

  
  
}

shinyApp(ui = ui, server = server)

#glitter cursor and cats 
