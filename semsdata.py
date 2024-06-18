#works with Brechtel datafiles from SEMS v5.3.0

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class SEMSdata():
    def __init__(self,path):
        self.path = path
        #import raw dataframe and stored as self.df
        self.df = pd.read_csv(path, skiprows = 53, delimiter = "\t") 
        
        #create dataframe to store processed data
        self.df_proc = pd.DataFrame()
        self.df_proc['StartDateTime'] = pd.to_datetime(self.df["#StartDate"].astype('str') + 
                                                       " " + self.df["StartTime"], format='%y%m%d %H:%M:%S')
             
        #bin diameter centers
        diameters = self.df.loc[:,'Bin_Dia1':'Bin_Dia60']
        self.diameters = diameters.values[0]
        
        #calculate bin diameter widths
        bin_dia_limits = np.zeros(len(self.diameters)+1)
        logbinsize = np.zeros(len(self.diameters))
        #calculate ends of bins
        for i in np.arange(1,len(bin_dia_limits)-1):
            bin_dia_limits[i] = (self.diameters[i-1]*self.diameters[i])**(1/2)
        bin_dia_limits[0] = self.diameters[0]**2/bin_dia_limits[1]
        bin_dia_limits[-1] = self.diameters[-1]**2/bin_dia_limits[-2]
        for i in np.arange(0,len(logbinsize)):
            logbinsize[i] = np.log10(bin_dia_limits[i + 1]) - np.log10(bin_dia_limits[i])
        #store values for future retrieval
        self.logbinsize = logbinsize
      
        #lists to hold concentrations and calculated parameters for each scan
        concentrations = []
        num_concs = []
        mass_concs = []
        total_num_concs = []
        total_mass_concs = []
        ave_num_diams = []
        ave_vol_diams = []
        
        #go through the dataframe, calculate, and store concentrations and params for each scan
        for scan in range(len(self.df.iloc[:,0])):
            concentration1 = self.df.loc[scan,'Bin_Conc1':'Bin_Conc60']
            concentrations.append(concentration1)
            num_conc = concentration1.values * self.logbinsize
            num_concs.append(num_conc)
            
            total_num_conc = np.sum(num_conc)
            total_num_concs.append(total_num_conc)
            
            mass_conc = num_conc*4/3*3.1415*((self.diameters/2)**3)*1.2*1e6*1e6*1e-21 #assumes density 1.2
            mass_concs.append(mass_conc)
            total_mass_conc = np.sum(mass_conc)
            total_mass_concs.append(total_mass_conc)
            
            ave_num_diam = np.sum(num_conc*self.diameters)/np.sum(num_conc)
            ave_num_diams.append(ave_num_diam)
            ave_vol_diam = np.sum((num_conc*self.diameters**3*3.14159/6)*self.diameters)/np.sum(num_conc*self.diameters**3*3.14159/6)
            ave_vol_diams.append(ave_vol_diam)
        
        #store lists in df_proc dataframe
        self.df_proc['Diams'] = pd.Series([self.diameters] * len(self.df_proc)) 
        self.df_proc['Num_Conc'] = num_concs
        self.df_proc['Mass_Conc'] = mass_concs
        self.df_proc['Total_Num_Conc'] = total_num_concs
        self.df_proc['Total_Mass_Conc'] = total_mass_concs
        self.df_proc['Ave_Num_Diam'] = ave_num_diams
        self.df_proc['Ave_Vol_Diam'] = ave_vol_diams
        #print basic information and the first five rows
        print("DataFrame with ",len(self.df_proc.index), " rows imported.")
        print(self.df_proc.head())
    
        
    def plotAll(self):
        for scan in range(len(self.df_proc.iloc[:,0])):
            plt.plot(self.diameters, self.df_proc['Num_Conc'][scan], '.')
            plt.xlabel('Diameter [nm]')
            plt.ylabel('Concentration [#/cm^3]')
            plt.title('Num. Conc. vs Diam. of Aerosol Particles of Scan '+ 
                      str(scan+1) +' at Time '+ self.df_proc['StartDateTime'][scan].strftime('%y%m%d %H:%M:%S'))
            plt.show()

    def plotNumConc(self, limits = None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            df = self.df_proc.loc[((self.df_proc['StartDateTime'] > start) & (self.df_proc['StartDateTime'] < end))]
        else:
            df = self.df_proc

        plt.figure()
        plt.plot(df['StartDateTime'],df['Total_Num_Conc'])
        plt.xlabel('Time')
        plt.ylabel('Num Conc [# cm^-3]')

    def plotMassConc(self, limits = None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            df = self.df_proc.loc[((self.df_proc['StartDateTime'] > start) & (self.df_proc['StartDateTime'] < end))]
        else:
            df = self.df_proc

        plt.figure()
        plt.plot(df['StartDateTime'],df['Total_Mass_Conc'])
        plt.xlabel('Time')
        plt.ylabel('Mass Conc [ug m^-3]')


    def calcNumConcMean(self,limits = None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            df = self.df_proc.loc[((self.df_proc['StartDateTime'] > start) & (self.df_proc['StartDateTime'] < end))]
        else:
            df = self.df_proc
            
        mean = np.round(df['Total_Num_Conc'].mean(),2)
        std = np.round(df['Total_Num_Conc'].std(),2)
        print(str(mean) + " +- " +  str(std) + " particles/cc")

    
    def calcMassConcMean(self, limits = None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            df = self.df_proc.loc[((self.df_proc['StartDateTime'] > start) & (self.df_proc['StartDateTime'] < end))]
        else:
            df = self.df_proc

        mean = np.round(df['Total_Mass_Conc'].mean(),2)
        std = np.round(df['Total_Mass_Conc'].std(),2)

        print(str(mean) + " +- " + str(std) + " ug/m^3")

    def createHeatmap(self,limits = None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            dfsliced = self.df_proc.loc[((self.df_proc['StartDateTime'] > start) & (self.df_proc['StartDateTime'] < end))]
        else:
            dfsliced = self.df_proc

        df = dfsliced['Num_Conc']
        dfindices = self.df_proc['Diams'][0]
        df = pd.DataFrame(df.tolist())
        df.columns = dfindices
        df['StartDateTime'] = dfsliced['StartDateTime']
        df = df.set_index('StartDateTime')
        df_transposed = df.transpose()

        # Create the heatmap
        #plt.figure(figsize=(10, 8))  # Adjust the figure size
        plt.figure()  # Adjust the figure size
        sns.heatmap(df_transposed, cmap='coolwarm')  # Create heatmap with annotations
        plt.gca().invert_yaxis()
        plt.title('Heatmap of DataFrame')  # Add title
        plt.show()  # Display the heatmap
    
