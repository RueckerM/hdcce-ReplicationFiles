"Application_Code_Paper.R" contains the replication code for the application.

For the code "Application_Code_Paper.R" to succesfully run the input data must be in the form of a dataframe with the following characteristics:

In terms of content:

(1) The first column has name "permno" and contains the stock identifier. (In our case these are the CRSP permno)
(2) The second column has name "date" and contains the date in format "%Y-%m-%d"
(3) The third column has name "ret" and contains the excess returns
(4) All the remaining columns contain the stock characteristics. 
 
The data should be sorted by company and by data within company. The file "data_structure_application.csv" shows the structure of the data employed in the application.
All the data values have been replaced with NAs due to issues with sharing the data as it is not publicly available.
 
