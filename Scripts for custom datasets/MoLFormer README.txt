FOR USING NEW SCRIPTS ON CUSTOM DATASETS:


You can run the appropriate MoLFormer procedures from the corresponding shell file. In order to do this, you need everything in the MoLFormer github https://github.com/IBM/molformer/tree/main, and you need to put your dataset in the ‘data’ former of molformer. 


It MUST be a folder containing three csv files, named train.csv, test.csv, and valid.csv. Python functions that split datasets from the papers we looked at can be found in the csvPrep_dataAug.py file of this repository https://github.com/egilman2/RIPSRelay.


Then, fill out the arguments of the shell accordingly, and run the shell.


Results and checkpoints will appear in the corresponding checkpoint folder – every line will give: the validation error of that epoch, the test error of the epoch, the overall minimum validation error, and the test error of the epoch with the overall minimum validation error.

To run the embeddings/predictions bash file, replace the name of the file at the top of the regression file with the name of the embeddings python file.