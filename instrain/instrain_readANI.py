#!/usr/bin/python


#################################################################
###                  Krys Kibler 2024-07-22                   ###
###                       PYTHON SCRIPT                       ###
### Purpose: obtain readani in figure 3 of instrain           ###
### https://github.com/MrOlm/inStrain/issues/131              ###
#################################################################

import argparse


# arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Script that extracts readANI from fig 3 of instrain"
    )
    parser.add_argument("--instrain_out", required=True, type=str)
    parser.add_argument("--csv", required=True, type=str)

    args = parser.parse_args()

    instrain_out = args.instrain_out
    csv = args.csv


# import instrain things
import inStrain
import inStrain.SNVprofile
import inStrain.plotting
import inStrain.plotting.mapping_plots


# code
IS = inStrain.SNVprofile.SNVprofile(f"{instrain_out}")
Mdb = inStrain.plotting.mapping_plots.prepare_read_ani_dist_plot(IS)
Mdb.to_csv(f"{csv}", index=False)
