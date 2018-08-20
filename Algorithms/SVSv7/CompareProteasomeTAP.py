#/usr/bin/python

import scipy
from scipy import stats

import pandas as pd

class CompareProteasomeTAP():       
    """
    Class: determines the percentile score
    ASSUMPTION: the higher the score of a peptide -> the higher the score of the percentile -> the better the peptide performs in a given process (e.g. better at being cleaved by proteasomal cleavage, better at being transported by TAP)
    """

    def __init__( self, path_prot_tap ):
        """
        Args:
            -path_prot_tap = string that is path to file that contains all proteasome 
                -the file should be tab-delimited
                -file should have columns "peptide", "proteasome_score", "tap_score"
                -NOTE: As of now, the file name is 
        """
        self.df_prot_tap = pd.read_csv( path_prot_tap, sep = '\t' )

    @staticmethod
    def calculate_percentileofscore( list_score, curr_score ):
        """
        Args:
            -list_score = an array (list or numpy array) of integers or floats. This should be a list of scores for either proteasomal cleavage or TAP transport
            -curr_score = float value. This should be a score for a prediction algorithm for either proteasomal cleavage or TAP transport
        Function: determine the 
        """
        #Method 1
        return scipy.stats.percentileofscore( list_score, curr_score, 'strict' )

    def calc_percentile_prot_tap( self, curr_score, pick_process = 1 ):
        """
        Calculates the percentile of the "curr_score" for either process proteasomal cleavage or TAP transport
        Args:
            -curr_score = float value. This should be a score for a prediction algorithm for either proteasomal cleavage or TAP transport
            -pick_process = integer that picks a process (proteasomal cleavage or TAP transport)
                -1 = chooses proteasomal cleavage
                -2 = chooses TAP transport
        """
        if pick_process == 1:
            return CompareProteasomeTAP.calculate_percentileofscore( self.df_prot_tap["proteasome_score"], curr_score )
        else:
            return CompareProteasomeTAP.calculate_percentileofscore( self.df_prot_tap["tap_score"], curr_score )

    def calc_percentile_dataframe( self, df_mhc, calc_perc = 3 ):
        """
        Calculate the percentile score for either proteasome score, TAP score, or both proteasome and TAP scores
        Args:
            -df_mhc = a pandas dataframe coming from class MHC_IEDB.retrieve_epitope_analysis(), where this dataframe contains a column titled 'peptide' that contains peptide sequence of interest
            -calc_perc = integer that determines which score to calculate a percentile score for
                -1 = only proteasome
                -2 = only TAP
                -3 = both proteasome & TAP
        NOTE: the score for proteasome & TAP are calculated with the score -log(IC50), therefore the higher the Proteasome/TAP score, the higher the cleavage likelihood & binding affinity (and therefore TAP transport), respectively.
        """
        def eval_percentile( row ):
            col_prot = "proteasome_score"
            col_tap = "tap_score"
            #see if proteasomal percentile should be calculated
            if calc_perc == 1 or calc_perc == 3:
                row["proteasome_percentile"] = self.calc_percentile_prot_tap( row[col_prot], 1 )
            #see if TAP transport percentile should be calculated
            if calc_perc == 2 or calc_perc == 3:
                row["tap_percentile"] = self.calc_percentile_prot_tap( row[col_tap], 2 )

            return row

        df_mhc = df_mhc.apply( eval_percentile, axis = 1 )

        return df_mhc