import pandas as pd
import matplotlib.pyplot as plt

class RHdata():
    def __init__(self,path):
        self.path = path
        #import raw dataframe and stored as self.df
        self.df = pd.read_csv(path, skiprows = 3, delimiter = "\t").dropna(axis=1)
        self.df['time'] = pd.to_datetime(self.df['time'])
    
        
    def plotT(self,limits=None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            df = self.df.loc[((self.df['time'] > start) & (self.df['time'] < end))]
        else:
            df = self.df

        plt.figure()
        plt.plot(df['time'],df['s0_T'])
        plt.xlabel('Time')
        plt.ylabel('T [deg C]')
        
        
    def plotRH(self,limits=None):
        if limits is not None:
            start = pd.to_datetime(limits[0])
            end = pd.to_datetime(limits[1])
            df = self.df.loc[((self.df['time'] > start) & (self.df['time'] < end))]
        else:
            df = self.df

        plt.figure()
        plt.plot(df['time'],df['s0_RH'])
        plt.xlabel('Time')
        plt.ylabel('RH [%]')