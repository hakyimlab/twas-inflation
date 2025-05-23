#! /usr/bin/env python3

# modules
import os
import sqlite3
import pandas as pd
import numpy as np
from scipy.stats import chi2
from scipy.stats import norm

# This is a simple function to correct inflation in TWAS studies
# We assume here you used the latest predictDB model whereby each gene has the inflation factor
# You should also have the sample size of the trait and the heritability.
# In a nut shell the function takes the phi, N and h2 to correct for inflation in the results
# This function has also been implemented in SPrediXcan and you just need to provide --gwas_h2 and --gwas_N arguments
# when runnning SprediXcan with the latest models and you will obtain calibrated results


#calibration functions
def get_phi(db):
    # get phi from SQLite database in PredictDB format
    try:
        conn = sqlite3.connect(db)
        extras = pd.read_sql_query("SELECT gene, phi FROM extra", conn)
        extras = extras.dropna(subset=['phi'])
        return extras
    except (sqlite3.DatabaseError, sqlite3.OperationalError) as e:
        # Log any database-related errors
        logging.error(f"An error occurred while accessing the database: {e}")
    finally:
        # Ensure the connection is closed
        if conn:
            conn.close()

def correct_inf_phi(xcan_df, predict_db, N, h2):
    # correct for inflation in Spredicxan  results
    #xcan_df: These are TWAS results and should contain the gene, pvalue and zscore columns
    extras = get_phi(predict_db)
    xcan_df = xcan_df.merge(extras, on='gene', how='inner')

    xcan_df['uncalibrated_pvalue'] = xcan_df['pvalue']
    xcan_df['uncalibrated_zscore'] = xcan_df['zscore']

    #QC: Replace negative phi values with 0
    xcan_df['phi'] = np.where(xcan_df['phi'] < 0, 0, xcan_df['phi'])
    
    # Calibration value
    denominator = 1 + (xcan_df['phi'] * N * h2)
    
    # calibrated z-score and pvalue
    xcan_df['zscore'] = xcan_df['zscore'] / np.sqrt(denominator)
    xcan_df['pvalue'] = 2 * norm.sf(abs(xcan_df['zscore']))
    
    logging.info("The pvalue and zscore have been calibrated successfully")
    return xcan_df
