Hi, Thanks for downloading this code. 
This is the readme file for the model associated with the paper:
Naze S., Bernard C. Jirsa V. K. (2015). Computational modeling of seizure dynamics using coupled neuronal networks: Factors shaping epileptoform activity.

All the files were supplied by the authors. Copyright Sebastien Naze, 2014. email: sebastien.naze@gmail.com

A class diagram from the model is enclosed in the file Class_diagram_global.pdf.

To run a simulation:
The script to be executed is popEpi_varParam.py.
The parameters taken by default are at the beginning of the script (parse.(...)) together with their explanation.
The nomenclature may differs slightly as compared to the maths equations published in the article (or some constant may be hard-coded) but their significations are written in the header.

Remarks:
- The code is far from optimized. This object-oriented architecture has been improvised in order to start computational modeling of neuron-like equation without previous knowledge in computational neuroscience and related software bundles (background in computer science and software engineering). A simulation of 40 neurons per population (80 in total) for 100s takes about 6h in a 3 GHz machine. 
- Another version of the code using vectorized equations is presented in the file populations_epilepton_HMR_ML_euler.py. 

For further information please contact:

Sebastien Naze
Institut de Neurosciences des Systemes - UMR 1106
27 Bd Jean Moulin
13005 Marseille, FRANCE

http://ins.medecine.univmed.fr/
email: sebastien.naze@gmail.com
