

Commands
--------

* ``tslGetUcscRnaAligns`` - Load PSL alignments from UCSC all_mrna or all_est tables into an SQLite database. (C program)
* ``tslGetEnsemblRnaAligns`` - Fetch cDNA alignments from Ensembl and load as PSLs into an SQLite3 database.
* ``tslGbffGetProblemCases`` - Scan genbank flat files looking for known problem libraries.
* ``tslGenbankProblemCasesLoad`` - Load problem case tab file generate ``gbffGetProblemCases`` by into an SQLite3 databases.
* ``tslLoadGenbankEvid`` - Build and load all GenBank evidence.
* ``tslCollectSupport`` - Collect support for GENCODE annotations.  Normally run in a cluster job.
* ``tslCollectSupportMkJobs`` - Generate parasol jobs to collect TSLs for GENCODE.
* ``tslCollectSupportJob`` - Job wrapper to run ``rslGencodeCollectSupport``.
* ``tslCollectSupportFinishJobs`` - Combine ``tslCollectSupport`` job results and store in an SQLite3 table.

Data flow
---------

#. GENCODE data
   - ``gencodeDb`` - database of annotations in Sqlite, built from Ensembl or UCSC databases, or maybe directly from Ensembl DB

#. Evidence sets as alignments
   - ``evidenceMeta`` - evidence metadata database (SqlLite or MySql)
   - ``evidenceData`` - alignments for each evidence data set, in tabix format
   
#. Support evaluation
   - ``supportEval`` - results of evidence evaluation for transcripts
