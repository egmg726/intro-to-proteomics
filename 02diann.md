---
title: 'Processing data with DIA-NN'
teaching: 10
exercises: 5
---


:::::::::::::::::::::::::::::::::::::: questions 

- How can I process my output files from a mass spectrometry-based proteomics experiment? 

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Learn how to process proteomics data using DIA-NN.

::::::::::::::::::::::::::::::::::::::::::::::::

## Processing proteomics data

There are many software options for processing mass spectrometry-based proteomics output. These software match mass spectra against an empirical or predicted spectral library to **identify peptides** and **infer the protein** to which they belong.

Some software you may have heard of include Spectronaut, Fragpipe, MaxQuant and DIA-NN. Each software have slightly different functionality and attributes, and are frequently updated. In this tutorial, we are demonstrating the use of **DIA-NN**, a software designed to process DIA proteomics data.

![Proteomics software like DIA-NN match observed spectra against a spectral library to identify peptides and assign them to a protein or protein group. Source: [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/proteomics/tutorials/introduction/slides.html#15)](episodes/fig/02diann_peptideMatching.png)

<br>

### Why DIA-NN?

Key features of DIA-NN software include:

- **Neural network–based scoring:** Improves identification confidence.
- **Interference correction:** Enhances quantification accuracy.
- **Library-free analysis:** Can work without a spectral library (generating one directly from DIA data).
- **High speed and efficiency:** Suitable for large-scale datasets.
- **Cross-run normalization:** Ensures consistent quantification across experiments.
- **Command-line executable:** DIA-NN can be run from a graphical user interface (GUI) or directly from the command line, making it easy to integrate into an automated pipeline.

DIA-NN is free to download. Detailed operating instructions and a link to install the latest version of DIA-NN can be found on the [GitHub page](https://github.com/vdemichev/DiaNN?tab=readme-ov-file).


## DIA-NN workflow

:::::callout
# Important Note

Due to the amount of time it will take to process the data, you do not need to install or operate DIA-NN for the purposes of this tutorial. Pre-processed data is provided, which has been prepared following the procedures described below.

:::::


### Part 1: Creating a spectral library

DIA-NN analysis consists of two steps. The first step is to **generate a predicted spectral library** for your experimental organism.

1. **Download a sequence database** from [UniProt](https://www.uniprot.org/uniprotkb?query=Human&facets=model_organism%3A9606%2Creviewed%3Atrue). 

    a. Go to https://www.uniprot.org/uniprotkb?query=*
    
    b. Select the `relevant organism` for your sample and `Reviewed (Swiss-Prot)` status to filter results.
    
    c. Click `Download`
    
    d. Change the format to `FASTA (canonical and isoform)`
    
    e. Select `No` under Compressed and download the FASTA file
    
:::::callout
# FOR SPECIFIC SEARCHES

If you need to search for a variant protein or peptide that differs from the standard organism proteome, you can manually edit the protein sequences in the UniProt FASTA in a text editor.

:::::
    
2. **Download a contaminant FASTA file** from [here](https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/tree/main). For most studies, the Universal Contaminants.fasta will be appropriate, but you may prefer to use one of the sample-type specific contaminant FASTAs depending on your experiment. ***More information about this resource and contaminants will be discussed in the next lesson.***

3. **Open the DIA-NN** GUI or a command line interface. If there are multiple versions of DIA-NN installed on your device, ensure to update `DIA-NN exe` with the file path to the version you wish to use.

4. **Add the UniProt and Contaminant FASTA files** you just downloaded. 

:::::caution

At the time of writing, DIA-NN 2.2.0 cannot read or write files saved in an external drive. Please ensure all input and output files are saved on a local drive to avoid error.

:::::

5. Examples are included below for the additional **settings we recommend changing from the default**. Depending on your experiment and research goals, you may wish to change additional settings from the default - read more in the [DIA-NN documentation](https://github.com/vdemichev/DiaNN?tab=readme-ov-file#changing-default-settings).

6. **Run**. DIA-NN will output a predicted spectral library in the `.speclib` file type.

![Example of the DIA-NN GUI used to generate a predicted spectral library. Settings which have been altered from the default are highlighted in red.](episodes/fig/02diann_DIANNgui_speclib.png)


``` bash
C:\DIA-NN\2.2.0\diann.exe ` # Replace with file path to DIA-NN on your device
--lib "" `
--out "C:\path\to\output\report.parquet" ` # Replace with file path to output report.parquet
--out-lib "C:\path\to\output\report-lib.parquet" ` # Replace with file path to output report-lib.parquet
--fasta "C:\path\to\UniprotFasta\.fasta" ` # Replace with file path to Uniprot FASTA
--fasta "C:\path\to\ContaminantsFasta\.fasta" ` # Replace with file path to contaminants FASTA
--threads 30 ` # Replace with the number of threads to use for this analysis
--verbose 1 `
--qvalue 0.01 `
--matrices  `
--gen-spec-lib `
--predictor `
--reannotate `
--fasta-search `
--min-fr-mz 200 `
--max-fr-mz 1800 `
--met-excision `
--min-pep-len 7 `
--max-pep-len 30 `
--min-pr-mz 300 `
--max-pr-mz 1800 `
--min-pr-charge 1 `
--max-pr-charge 4 `
--cut K*,R* `
--missed-cleavages 1 `
--unimod4 `
--rt-profiling 
```

This spectral library can be re-used for all analyses of samples from the same organism. It is recommended to re-download the UniProt FASTA and re-run this step approximately once a month to ensure your library is up-to-date.

<br>

### Part 2: Analysis

The second step is to **search your raw files (mass spectrometry output files)** against the spectral library.

1. **Open the DIA-NN** GUI or a command line interface. If there are multiple versions of DIA-NN installed on your device, ensure to update `DIA-NN exe` with the file path to the version you wish to use.

2. **Upload your raw files**.

    - If using the GUI, under `Input` select `Raw` then select all of the files you wish to search.
    - If using the command line, you can use `--f` multiple times to list each file individually, or (recommended) use `--dir [folder]` to instruct DIA-NN to process all raw files saved in a particular folder.


3. **Upload your spectral library**.

:::::callout
# Empirical vs predicted spectral library

You may wish to generate an empirical spectral library from your data, and use this for future analyses. An empirical library will be smaller than a predicted library and therefore increase the speed of analysis; however, may result in fewer identifications. We only recommend using an empirical library if generated from high quality data.

You can generate an empirical library from your samples by selecting **Generate spectral library** and specifying the location of its output.

:::::

4. **Adjust settings** according to your experiment and research goals. Again, examples are included below for the additional settings we recommend changing from the default. You may wish to change additional settings from the default - read more in the [DIA-NN documentation](https://github.com/vdemichev/DiaNN?tab=readme-ov-file#changing-default-settings).

    - In this example, the **protease** Trypsin was used to generate peptides during sample preparation.
    - Leaving **Mass accuracy** and **MS1 accuracy** on the default 0 means DIA-NN will optimise these parameters automatically for the first run in the experiment and then reuse the optimised settings for other runs. The DIA-NN documentation recommends different settings depending on the  type of mass spectrometer used to generate your data. In this example, data were generated on a **TripleTOF 6600**, so we set Mass accuracy to 20.0 and MS1 accuracy to 12.0
    - Set the number of **Threads** to be no more than the number of logical cores on your computer.
    - In this experiment, we have set the **maximum number of variable modifications** to 0; however, it is common for this setting to be relaxed in the range of 1 - 3.
    - DIA-NN's default setting is to activate **match-between-runs**. From the documentation: *In MBR mode, DIA-NN does two passes over the data. During the first pass, it creates an empirical spectral library from the data. During the second pass, it reanalyses the experiment with this empirical library, which may result in much improved identification numbers, data completeness and quantification accuracy.*

![Example of the DIA-NN GUI used to analyse the data described in this tutorial. Settings which have been altered from the default are highlighted in red.](episodes/fig/02diann_DIANNgui_analysis.png)


``` bash
C:\DIA-NN\2.2.0\diann.exe ` # Replace with file path to DIA-NN on your device
--dir "C:\path\to\folder\of\raw\files" ` # Replace with file path to folder containing all raw files to be analysed
--lib "C:\path\to\SpectralLibrary\report-lib.predicted.speclib" ` # Replace with file path to spectral library
--out "C:\path\to\output\report.parquet" ` # Replace with file path to output report.parquet
--threads 30 ` # Replace with the number of threads to use for this analysis
--verbose 1 `
--qvalue 0.01 `
--matrices  `
--unimod4 `
--mass-acc 20 `
--mass-acc-ms1 12.0 `
--reanalyse `
--rt-profiling 
```

:::::challenge

If you want to run DIA-NN from the command line but are unsure how to write the code, you can update your settings in the GUI, click 'Run', and copy the commands output in the log. That's how we got the example code above! **Can you match each line of code to the corresponding setting in the GUI?**

:::solution

Check the description of available commands in the [DIA-NN documentation](https://github.com/vdemichev/DiaNN?tab=readme-ov-file#command-line-reference)

:::
:::::

## DIA-NN output

DIA-NN output includes several files, including:

| Output file        | Description                                                                         |
| -----              | ------------------                                                                        |
| **log.txt**        | A log of what DIA-NN has run. You should always use `ctrl + F` to check for any 'ERROR' or 'WARNING' messages issued by DIA-NN during your run. |
| **.parquet**       | This is what we will be using for subsequent analyses in this tutorial. |
| **pg_matrix.tsv**  | Matrix of protein values generated by DIA-NN. |
| **pr_matrix.tsv**  | Matrix of peptide values generated by DIA-NN. |
| **trends.pdf**     | A PDF report listing the number of peptides, proteins, and other quality control measures detected for each sample. |

::::::::::::::challenge

If DIA-NN is run through the command line rather than the GUI, the trends.pdf and runs.pdf reports are not automatically generated. Using the documentation provided in the [DIA-NN GitHub page](https://github.com/vdemichev/DiaNN?tab=readme-ov-file#changing-default-settings), can you work out how to generate this report via the command line?

:::solution

Run the `diann-stats.py` Python script provided when you download DIA-NN.


``` bash
python C:\DIA-NN\2.2.0\diann-stats.py C:\path\to\report.parquet
```

:::
::::::::::::::

We have now processed our raw files, identifying the peptides present in our samples and inferring the proteins to which they belong.

:::::callout
# Learn more about DIANN

For more background information, see the original publication:  
*Demichev, V., et al. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. Nature Methods.*  [https://doi.org/10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x).

For further details about the different DIA-NN settings, see the [documentation](https://github.com/vdemichev/DiaNN?tab=readme-ov-file).

If you come across any issues, the creators are fairly active at responding to [issues](https://github.com/vdemichev/DiaNN/issues) and [discussions](https://github.com/vdemichev/DiaNN/discussions) on their GitHub.
:::::

::::::::::::::::::::::::::::::::::::: keypoints 

- Raw files from a mass-spectrometry based DIA proteomics experiment can be processed using DIA-NN (peptide identification and protein inference).
- The DIA-NN workflow has two key steps: generated a predicted spectral library, and processing raw files.
- The [DIA-NN documentation](https://github.com/vdemichev/DiaNN?tab=readme-ov-file) contains detailed information about changing default settings as appropriate for your experiment.

::::::::::::::::::::::::::::::::::::::::::::::::



