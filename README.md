# Challenge_Epigene_Labs

The following project aims at preprocessing data from a clinical study in order to make them comparable to other data form other studies. The project contains the followings :

An input directories which contains clinical data to preprocess for 2 samples from one experiment. The input directory needs to be unzipped before performing the script:

    a metadata file 
    a barcodes file for each sample : 
            contains information about each cells
    a features file for each sample :
            contains information about each genes
    a matrix file for each sample :
            a counting matrix of reads for each genes of each cell

An Outpout_Ex directory that contains the file we are supposed to get after the preprocessing

An Output directory that contains the file we get after the preprocessing

A Figures directory that contains figures comparing the expected output to the one we get

A Preprocessing.py file
    Perform the preprocessing of the clinical data and export it as an AnData object in a .h5ad file

A Figure.py file
    Plot all the figures to compare the expected output to the one we have

