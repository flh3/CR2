The SPSS extension to accompany Huang, F. & Li, X. (2021). Using cluster robust standard errors when analyzing group randomized trials with few clusters. Behavior Research Methods.

For computing cluster robust standard errors with degrees of freedom adjustments in SPSS.

## To install: 

To download, just click on the *.spd file and then choose the download option on the following screen. 
To install, **DO NOT** double click on the *.spd file.
To install, in SPSS, go to Extensions -> Utilities -> Install Custom Dialog (Compatibility Mode)... select the downloaded file and click OK to install.

## Usage: 

A new item is now added under Analyze -> Regression -> Cluster Robust... [for now, can only be used for linear regression]

If you have a large number of clusters (e.g., > 50), just use the CR0 (Liang & Zeger, 1986) or CR1 option with the HLM degrees of freedom (this makes it run faster).

If you have few clusters (e.g., < 50), use the CR2 option together with the Bell & McCaffrey dof adjustment. 

## To uninstall:

In SPSS, choose Extensions -> Utilities -> Custom Dialogue Builder (Compatibility mode)-> File -> Uninstall (select the CLUSTER... file)

## Updates

1.01 Can have L1 or L2 left blank (no predictors at that level).

