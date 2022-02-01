#!\bin\bash

# --------------------------------------------------
# Download data from David Koes website.
# --------------------------------------------------

wget http://bits.csb.pitt.edu/files/gnina1.0_paper/crossdocked_all_data.tar.gz
tar -xvf crossdocked_all_data.tar.gz

# Message about steps to take before runing the next script
echo "# --------------------------------------------------"
echo "  >>> Protein FKB1A/1F40_PRO.pdb need to be manually"
echo "  >>> amended to have only one model in the PDB file"
echo "# --------------------------------------------------"
