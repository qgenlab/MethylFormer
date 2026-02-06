from abc import ABC, abstractmethod
import pandas as pd
import numpy as np

class MethylationAdapter(ABC):
    def __init__(self, df: pd.DataFrame, suffixe:str = "0"):
        self.df = df
        self.suffixe = suffixe
    @abstractmethod
    def get_name(self):
        pass
    @abstractmethod
    def get_chrom(self) -> pd.Series:
        pass
    @abstractmethod
    def get_start(self) -> pd.Series:
        pass
    @abstractmethod
    def get_end(self) -> pd.Series:
        pass
    @abstractmethod
    def get_meth_diff(self) -> pd.Series:
        pass
    @abstractmethod
    def get_significance(self) -> pd.Series:
        pass
    def get_bed_format(self) -> pd.DataFrame:
        return pd.DataFrame({
            'chrom': self.get_chrom(),
            'start': self.get_start(),
            'end': self.get_end(),
            'score': self.get_significance(),
            'diff': self.get_meth_diff()
        })

